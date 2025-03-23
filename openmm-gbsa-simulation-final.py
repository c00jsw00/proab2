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

def find_disulfide_bonds(topology, positions, distance_cutoff=0.3):
    """
    Automatically detect disulfide bonds based on distance between sulfur atoms
    
    Parameters:
    -----------
    topology : openmm.app.Topology
        Topology object of the protein
    positions : list or array
        Positions of all atoms
    distance_cutoff : float
        Maximum distance between sulfur atoms to be considered a disulfide bond (in nm)
        
    Returns:
    --------
    list of tuples
        Each tuple contains (chain1, resid1, chain2, resid2) for a disulfide bond
    """
    print("Scanning for disulfide bonds...")
    disulfides = []
    
    # Collect all CYS SG atoms
    cys_sg_atoms = []
    for chain in topology.chains():
        for res in chain.residues():
            if res.name == 'CYS':
                for atom in res.atoms():
                    if atom.name == 'SG':
                        cys_sg_atoms.append((chain.id, res.id, atom, res))
    
    # Check distances between all SG atoms
    num_cys = len(cys_sg_atoms)
    print(f"Found {num_cys} cysteine residues")
    
    # Convert distance_cutoff to nm if positions have units
    try:
        # Check if positions have units
        has_units = hasattr(positions[0][0], 'unit') or hasattr(positions[0][0], 'value_in_unit')
    except (IndexError, TypeError):
        has_units = False
        
    # Convert distance cutoff to proper units if needed
    if has_units:
        try:
            from openmm import unit as openmm_unit
            distance_cutoff = distance_cutoff * openmm_unit.nanometer
        except ImportError:
            try:
                from simtk import unit as simtk_unit
                distance_cutoff = distance_cutoff * simtk_unit.nanometer
            except ImportError:
                print("WARNING: Could not import unit module, units might be inconsistent")
    
    for i in range(num_cys):
        chain_i, resid_i, atom_i, res_i = cys_sg_atoms[i]
        pos_i = positions[atom_i.index]
        
        for j in range(i+1, num_cys):
            chain_j, resid_j, atom_j, res_j = cys_sg_atoms[j]
            pos_j = positions[atom_j.index]
            
            # Calculate distance - handle units
            try:
                # Try with units
                dx = pos_i[0] - pos_j[0]
                dy = pos_i[1] - pos_j[1]
                dz = pos_i[2] - pos_j[2]
                distance = (dx*dx + dy*dy + dz*dz)**0.5
                
                # Compare with units
                if distance < distance_cutoff:
                    disulfides.append((chain_i, resid_i, chain_j, resid_j))
                    print(f"  Detected disulfide bond: Chain {chain_i} CYS {resid_i} - Chain {chain_j} CYS {resid_j} (distance: {distance})")
            except (TypeError, AttributeError):
                # If units cause problems, try to strip units
                try:
                    if has_units:
                        # Try to get values in nanometers
                        from openmm import unit as openmm_unit
                        pos_i_nm = [v.value_in_unit(openmm_unit.nanometer) for v in pos_i]
                        pos_j_nm = [v.value_in_unit(openmm_unit.nanometer) for v in pos_j]
                    else:
                        # Assume values are already in nanometers
                        pos_i_nm = [float(v) for v in pos_i]
                        pos_j_nm = [float(v) for v in pos_j]
                        
                    dx = pos_i_nm[0] - pos_j_nm[0]
                    dy = pos_i_nm[1] - pos_j_nm[1]
                    dz = pos_i_nm[2] - pos_j_nm[2]
                    distance_nm = (dx*dx + dy*dy + dz*dz)**0.5
                    
                    if distance_nm < float(distance_cutoff):
                        disulfides.append((chain_i, resid_i, chain_j, resid_j))
                        print(f"  Detected disulfide bond: Chain {chain_i} CYS {resid_i} - Chain {chain_j} CYS {resid_j} (distance: {distance_nm:.3f} nm)")
                except Exception as e:
                    print(f"  WARNING: Error calculating distance between CYS {resid_i} and CYS {resid_j}: {e}")
                    continue
    
    return disulfides

def add_disulfide_bonds(system, topology, disulfide_bonds):
    """
    Add disulfide bonds to the system
    
    Parameters:
    -----------
    system : openmm.System
        The OpenMM system object
    topology : openmm.app.Topology
        Topology object of the protein
    disulfide_bonds : list of tuples
        List of disulfide bonds as returned by find_disulfide_bonds
    
    Returns:
    --------
    bool
        True if all disulfide bonds were successfully added
    """
    if not disulfide_bonds:
        print("No disulfide bonds to add")
        return True
    
    print(f"Adding {len(disulfide_bonds)} disulfide bonds to the system")
    
    # Function to find the SG atom index of a specific CYS residue
    def find_sg_atom(chain_id, res_id):
        for chain in topology.chains():
            if chain.id == chain_id:
                for res in chain.residues():
                    if res.id == res_id and res.name == 'CYS':
                        for atom in res.atoms():
                            if atom.name == 'SG':
                                return atom.index
        return None
    
    # Add harmonic bond constraints for each disulfide
    force = None
    
    # Try to find an existing HarmonicBondForce
    for i in range(system.getNumForces()):
        if system.getForce(i).__class__.__name__ == 'HarmonicBondForce':
            force = system.getForce(i)
            break
    
    # If no HarmonicBondForce exists, create one
    if force is None:
        try:
            # For newer OpenMM versions
            from openmm import HarmonicBondForce
            force = HarmonicBondForce()
            system.addForce(force)
        except ImportError:
            # For older OpenMM versions
            try:
                from simtk.openmm import HarmonicBondForce
                force = HarmonicBondForce()
                system.addForce(force)
            except ImportError:
                print("ERROR: Could not create HarmonicBondForce - disulfide bonds will not be added")
                return False
    
    # Constants for the disulfide bond
    length = 0.204  # nm, typical S-S bond length
    k = 1000.0      # kJ/mol/nm^2, force constant
    
    success = True
    for bond in disulfide_bonds:
        chain1, resid1, chain2, resid2 = bond
        sg1 = find_sg_atom(chain1, resid1)
        sg2 = find_sg_atom(chain2, resid2)
        
        if sg1 is not None and sg2 is not None:
            force.addBond(sg1, sg2, length, k)
            print(f"  Added disulfide bond between atom indices {sg1} and {sg2}")
        else:
            print(f"  WARNING: Could not find SG atoms for CYS {chain1}:{resid1} - {chain2}:{resid2}")
            success = False
    
    return success

def run_gbsa_simulation(input_file, output_prefix, simulation_steps=500000, 
                       report_interval=1000, gb_model='OBC2', ss_cutoff=0.3, manual_ss_bonds=None):
    """
    Run GBSA implicit solvent simulation on a protein structure
    
    Parameters:
    -----------
    input_file : str
        Path to the input structure file (CIF or PDB)
    output_prefix : str
        Prefix for output files
    simulation_steps : int
        Number of simulation steps to run
    report_interval : int
        Interval for reporting and saving trajectory frames
    gb_model : str
        GB model to use ('OBC1', 'OBC2', or 'HCT')
    ss_cutoff : float
        Distance cutoff for automatic disulfide bond detection (nm)
    manual_ss_bonds : list of tuples, optional
        Manually specified disulfide bonds, each as (chain1, resid1, chain2, resid2)
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
            
    # Detect and add disulfide bonds
    if manual_ss_bonds:
        print("Using manually specified disulfide bonds:")
        for bond in manual_ss_bonds:
            print(f"  Manual disulfide bond: Chain {bond[0]} CYS {bond[1]} - Chain {bond[2]} CYS {bond[3]}")
        disulfide_bonds = manual_ss_bonds
    else:
        # Automatically detect disulfide bonds
        disulfide_bonds = find_disulfide_bonds(pdb.topology, pdb.positions, distance_cutoff=ss_cutoff)
    if disulfide_bonds:
        add_disulfide_bonds(system, pdb.topology, disulfide_bonds)
    
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
    parser.add_argument('--ss-cutoff', type=float, default=0.3, 
                      help='Distance cutoff for automatic disulfide bond detection (nm)')
    
    args = parser.parse_args()
    
    run_gbsa_simulation(args.cif, args.output, args.steps, args.report, args.gbmodel, args.ss_cutoff)
