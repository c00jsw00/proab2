#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import numpy as np
from datetime import datetime
from sys import stdout

try:
    # Try new import style first (OpenMM 7.6+)
    import openmm as mm
    from openmm import unit
    from openmm.app import PDBFile, PDBReporter, DCDReporter, StateDataReporter
    from openmm.app import Modeller, ForceField, Simulation, NoCutoff, HBonds
    USING_NEW_OPENMM = True
    
    # Check if we can get the GB models from the new-style import
    try:
        from openmm.app import HCT, OBC1, OBC2
    except ImportError:
        # Define GB models using the app.internal namespace
        from openmm.app.internal.customgbforces import GBSAOBCForce
        HCT = GBSAOBCForce.HCT
        OBC1 = GBSAOBCForce.OBC1
        OBC2 = GBSAOBCForce.OBC2
    
except ImportError:
    # Fall back to old import style (pre-7.6)
    from simtk import unit
    import simtk.openmm as mm
    from simtk.openmm.app import PDBFile, PDBReporter, DCDReporter, StateDataReporter
    from simtk.openmm.app import Modeller, ForceField, HCT, OBC1, OBC2
    from simtk.openmm.app import Simulation, NoCutoff, HBonds
    USING_NEW_OPENMM = False

# Try to import PDBFixer if available
try:
    from pdbfixer import PDBFixer
    HAVE_PDBFIXER = True
except ImportError:
    HAVE_PDBFIXER = False
    print("PDBFixer not found. Some structure preparation features will be limited.")
    print("Consider installing PDBFixer: pip install pdbfixer")

def fix_structure(input_file, output_file=None):
    """
    Use PDBFixer to fix structural issues in a PDB or CIF file
    """
    if not HAVE_PDBFIXER:
        print("PDBFixer not available. Skipping structure fixing.")
        return input_file
        
    if output_file is None:
        base, ext = os.path.splitext(input_file)
        output_file = f"{base}_fixed.pdb"
    
    print(f"Fixing structure with PDBFixer: {input_file}")
    
    # Initialize fixer with the input file
    fixer = PDBFixer(filename=input_file)
    
    # Find missing residues
    print("Finding missing residues...")
    fixer.findMissingResidues()
    
    # Find missing atoms
    print("Finding missing heavy atoms...")
    fixer.findMissingAtoms()
    
    # Add missing atoms
    print("Adding missing atoms...")
    fixer.addMissingAtoms()
    
    # Add missing hydrogens
    print("Adding missing hydrogens...")
    fixer.addMissingHydrogens(7.0)  # pH 7.0
    
    # Write the fixed structure
    print(f"Writing fixed structure to {output_file}")
    with open(output_file, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    return output_file

def run_gbsa_simulation(input_file, output_prefix, simulation_steps=500000, 
                       report_interval=1000, gb_model='OBC2'):
    """
    Run GBSA implicit solvent simulation on a protein structure
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # First, fix the structure using PDBFixer if available
    fixed_pdb = f"{output_prefix}_fixed.pdb"
    fixed_pdb = fix_structure(input_file, fixed_pdb)
    
    print(f"Starting GBSA simulation with structure: {fixed_pdb}")
    start_time = datetime.now()
    
    # Load the structure
    try:
        pdb = PDBFile(fixed_pdb)
    except Exception as e:
        print(f"Error loading fixed PDB: {e}")
        print("Trying to load original file...")
        try:
            # Try loading original file
            if input_file.lower().endswith('.cif'):
                try:
                    from openmm.app import PDBxFile
                    pdb = PDBxFile(input_file)
                except ImportError:
                    from simtk.openmm.app import PDBxFile
                    pdb = PDBxFile(input_file)
            else:
                pdb = PDBFile(input_file)
                
            # Create modeller to add hydrogens
            modeller = Modeller(pdb.topology, pdb.positions)
            modeller.addHydrogens()
            pdb = modeller  # Use modeller as our "pdb" object
        except Exception as e2:
            print(f"Error loading original file: {e2}")
            sys.exit(1)
    
    # Set up force field
    print("Setting up force field...")
    
    # Try different force field combinations
    success = False
    force_field_options = [
        # Try these combinations in order
        ['amber14-all.xml', 'implicit/gbn2.xml'],
        ['amber14-all.xml', 'implicit/gbsaobc2.xml'],
        ['amber14-all.xml', 'implicit/obc.xml'],
        ['amber14-all.xml']  # Last resort - just protein force field
    ]
    
    for ff_files in force_field_options:
        try:
            print(f"Trying force field: {', '.join(ff_files)}")
            forcefield = ForceField(*ff_files)
            success = True
            break
        except Exception as e:
            print(f"Error with force field {ff_files}: {e}")
    
    if not success:
        print("All force field options failed. Exiting.")
        sys.exit(1)
    
    # Define GB models
    gbsa_methods = {
        'HCT': HCT,
        'OBC1': OBC1,
        'OBC2': OBC2
    }
    
    if gb_model not in gbsa_methods:
        print(f"Invalid GB model: {gb_model}. Using OBC2 as default.")
        gb_model = 'OBC2'
    
    gb_method = gbsa_methods[gb_model]
    print(f"Creating system with {gb_model} implicit solvent model...")
    
    # Create system - handle different OpenMM versions
    system = None
    try:
        # Different parameter naming for different OpenMM versions
        if USING_NEW_OPENMM:
            # Try new version parameter name
            print("Using new OpenMM parameter format...")
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=NoCutoff,
                constraints=HBonds,
                implicitSolventSaltConc=0.0*unit.moles/unit.liter,
                implicit_solvent=gb_method  # New parameter name in newer OpenMM
            )
        else:
            # Use old version parameter name
            print("Using legacy OpenMM parameter format...")
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=NoCutoff,
                constraints=HBonds,
                implicitSolvent=gb_method  # Old parameter name
            )
    except (ValueError, TypeError) as e:
        print(f"Error creating system with explicit implicit solvent: {e}")
        print("Trying without specifying implicit solvent explicitly...")
        try:
            # Last resort - create system without specifying implicit solvent
            system = forcefield.createSystem(
                pdb.topology,
                nonbondedMethod=NoCutoff,
                constraints=HBonds
            )
            print("Success - created system without explicit implicit solvent parameter")
            print("Warning: May not be using the intended solvation model")
        except Exception as e2:
            print(f"Final system creation attempt failed: {e2}")
            sys.exit(1)
    
    # Set up integrator
    temp = 300 * unit.kelvin
    friction = 1.0 / unit.picosecond
    timestep = 2.0 * unit.femtosecond
    integrator = mm.LangevinIntegrator(temp, friction, timestep)
    
    # Set up simulation
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    
    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # Set up reporters
    pdb_reporter = PDBReporter(f"{output_prefix}_trajectory.pdb", report_interval)
    dcd_reporter = DCDReporter(f"{output_prefix}_trajectory.dcd", report_interval)
    data_reporter = StateDataReporter(
        f"{output_prefix}_data.csv", report_interval,
        step=True, time=True, potentialEnergy=True, temperature=True, 
        progress=True, totalSteps=simulation_steps, 
        remainingTime=True, speed=True
    )
    console_reporter = StateDataReporter(
        stdout, report_interval*10,
        step=True, time=True, progress=True, potentialEnergy=True,
        temperature=True, remainingTime=True, speed=True,
        totalSteps=simulation_steps
    )
    
    simulation.reporters.append(pdb_reporter)
    simulation.reporters.append(dcd_reporter)
    simulation.reporters.append(data_reporter)
    simulation.reporters.append(console_reporter)
    
    # Save minimized structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(f"{output_prefix}_minimized.pdb", 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    
    # Run simulation
    print(f"Running {simulation_steps} steps of GBSA simulation...")
    simulation.step(simulation_steps)
    
    # Save final structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(f"{output_prefix}_final.pdb", 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    
    end_time = datetime.now()
    elapsed = end_time - start_time
    print(f"Simulation completed in {elapsed}")
    print(f"Output files saved with prefix: {output_prefix}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run OpenMM GBSA simulation on a protein structure')
    parser.add_argument('--cif', required=True, help='Input CIF or PDB structure file')
    parser.add_argument('--output', default='output/sim', help='Output prefix for simulation files')
    parser.add_argument('--steps', type=int, default=500000, help='Number of simulation steps')
    parser.add_argument('--report', type=int, default=1000, help='Reporting interval')
    parser.add_argument('--gbmodel', default='OBC2', choices=['HCT', 'OBC1', 'OBC2'], 
                      help='GB model to use')
    
    args = parser.parse_args()
    
    run_gbsa_simulation(args.cif, args.output, args.steps, args.report, args.gbmodel)
