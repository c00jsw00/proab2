#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime
import shutil
from Bio.PDB import PDBIO, MMCIFParser, PDBParser
import Bio.PDB
import freesasa
import csv
from collections import defaultdict

def parse_range(range_str):
    """Parse a range string like '1-31' into a list of numbers"""
    start, end = map(int, range_str.split('-'))
    return range(start, end + 1)

def read_config(config_file):
    """Read configuration file with input files and analysis parameters"""
    config = {
        'input_files': [],
        'cdr_residues': {},
        'ab_lock_residues': [],
        'substrate_residues': []
    }
    
    current_section = None
    with open(config_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            if line.endswith(':'):
                current_section = line[:-1].lower()
                continue
                
            if current_section == 'pdb_files':
                config['input_files'].append(line)
            elif current_section == 'cdr_residues':
                cdr_name, rest = line.split(':')
                chain, range_str = rest.strip().split()
                residue_range = parse_range(range_str)
                config['cdr_residues'][cdr_name] = [(chain, i) for i in residue_range]
            elif current_section == 'ab_lock_residues':
                chain, range_str = line.split(':')
                chain = chain.strip()
                residue_range = parse_range(range_str.strip())
                config['ab_lock_residues'].extend([(chain, i) for i in residue_range])
            elif current_section == 'substrate_residues':
                chain, range_str = line.split(':')
                chain = chain.strip()
                residue_range = parse_range(range_str.strip())
                config['substrate_residues'].extend([(chain, i) for i in residue_range])
    
    return config

class FrameAnalyzer:
    def __init__(self, pdb_file, model=None):
        """Initialize with PDB file and optional pre-loaded model"""
        if not os.path.exists(pdb_file) and model is None:
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        self.pdb_file = pdb_file
        
        if model is not None:
            # Use pre-loaded model
            self.model = model
        else:
            # Load model from file
            parser = Bio.PDB.PDBParser(QUIET=True)
            structure = parser.get_structure('frame', pdb_file)
            self.model = structure[0]
    
    def calc_blocking_rate(self, cdr_residues, ab_lock_residues):
        """
        Calculate blocking rate based on CDR coverage
        Args:
            cdr_residues: Dictionary of CDR regions and their residues
            ab_lock_residues: List of tuples (chain_id, residue_number)
        Returns:
            Dictionary of blocking rates for each CDR and total
        """
        results = {}
        
        # Calculate blocking rate for each CDR region
        for cdr_name, residues in cdr_residues.items():
            covered = 0
            total = len(residues)
            
            for cdr_res in residues:
                # Get CDR residue atoms
                cdr_atoms = self.get_residue_atoms(cdr_res)
                if not cdr_atoms:
                    continue
                
                # Check if covered by Ab lock
                is_covered = False
                for lock_res in ab_lock_residues:
                    lock_atoms = self.get_residue_atoms(lock_res)
                    if not lock_atoms:
                        continue
                    
                    # Check if any atom pair is within 4Å and above 120 degrees
                    for cdr_atom in cdr_atoms:
                        for lock_atom in lock_atoms:
                            dist = self.calc_distance(cdr_atom, lock_atom)
                            angle = self.calc_angle(cdr_atom, lock_atom)
                            if dist <= 4.0 and angle >= 120:
                                is_covered = True
                                break
                        if is_covered:
                            break
                    if is_covered:
                        break
                        
                if is_covered:
                    covered += 1
            
            # Calculate rate for this CDR
            rate = (covered/total) * 100 if total > 0 else 0
            results[cdr_name] = rate
            
        # Calculate average of all CDR rates
        cdr_rates = [rate for cdr_name, rate in results.items()]
        results['total'] = sum(cdr_rates) / len(cdr_rates) if cdr_rates else 0
        
        return results

    def calc_pasa(self, substrate_residues, probe_radius=40.0):
        """
        Calculate Protease Accessible Surface Area using FreeSASA
        Args:
            substrate_residues: List of tuples (chain_id, residue_number)
            probe_radius: Radius of probe in Å (default 40.0Å for MMP-2)
        Returns:
            Dictionary with PASA values for each chain and total
        """
        # Group residues by chain
        chain_ranges = {}
        for chain_id, res_num in substrate_residues:
            if chain_id not in chain_ranges:
                chain_ranges[chain_id] = []
            chain_ranges[chain_id].append(res_num)

        results = {}
        
        try:
            # Try different methods to set probe radius depending on freesasa version
            try:
                # Method 1: Using Parameters class
                parameters = freesasa.Parameters()
                # Try different attribute names
                if hasattr(parameters, 'probe_radius'):
                    parameters.probe_radius = probe_radius
                elif hasattr(parameters, 'probeRadius'):
                    parameters.probeRadius = probe_radius
                
                # Create structure and calculate with parameters
                structure = freesasa.Structure(self.pdb_file)
                result = freesasa.calc(structure, parameters)
                
            except (AttributeError, TypeError):
                # Method 2: Using default configuration but with classifier for larger probe
                try:
                    # Create classifier with larger probe
                    classifier = freesasa.Classifier()
                    # Try to set the probe radius if possible
                    if hasattr(classifier, 'setProbeRadius'):
                        classifier.setProbeRadius(probe_radius)
                    
                    structure = freesasa.Structure(self.pdb_file, classifier)
                    result = freesasa.calc(structure)
                    
                except (AttributeError, TypeError):
                    # Method 3: Just use default parameters
                    print(f"Warning: Could not set probe radius to {probe_radius}Å. Using default value.")
                    structure = freesasa.Structure(self.pdb_file)
                    result = freesasa.calc(structure)
            
            # Calculate SASA for each chain's residue range
            for chain_id, residues in chain_ranges.items():
                residues.sort()  # Sort residues for proper range selection
                # Create selection expression using Pymol syntax
                selection_str = f"s_{chain_id}, chain {chain_id} and resi "
                selection_str += "+".join([str(r) for r in residues])
                
                # Calculate area for selected residues
                areas = freesasa.selectArea([selection_str], structure, result)
                area = areas.get(f"s_{chain_id}", 0)
                results[chain_id] = area
            
            # Calculate total PASA as sum of chain PASAs
            results['total'] = sum(results.values()) if results else 0
                
        except Exception as e:
            print(f"Error calculating SASA: {e}")
            import traceback
            traceback.print_exc()
            return {'total': 0}
                
        return results

    def get_residue(self, res_id):
        """Get residue by chain ID and residue number"""
        chain_id, res_num = res_id
        try:
            for chain in self.model:
                if chain.id == chain_id:
                    for residue in chain:
                        if residue.id[1] == res_num:
                            return residue
        except Exception as e:
            print(f"Error accessing residue {chain_id}:{res_num}: {e}")
        return None

    def get_residue_atoms(self, res_id):
        """Get atoms for a residue by ID"""
        residue = self.get_residue(res_id)
        return list(residue.get_atoms()) if residue else []

    def calc_distance(self, atom1, atom2):
        """Calculate distance between two atoms"""
        diff_vector = atom1.coord - atom2.coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

    def calc_angle(self, atom1, atom2):
        """Calculate angle between two atoms relative to vertical"""
        vector = atom2.coord - atom1.coord
        vertical = np.array([0, 0, 1])
        cos_angle = np.dot(vector, vertical)/(np.linalg.norm(vector)*np.linalg.norm(vertical))
        return np.arccos(cos_angle) * 180/np.pi

class TrajectoryAnalyzer:
    def __init__(self, trajectory_file, frames_to_analyze=0):
        """
        Initialize with a trajectory file path
        
        Args:
            trajectory_file: Path to PDB trajectory file
            frames_to_analyze: Number of frames to analyze (0 for all frames)
        """
        if not os.path.exists(trajectory_file):
            raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
        
        self.trajectory_file = trajectory_file
        self.parser = Bio.PDB.PDBParser(QUIET=True)
        
        # Load all models/frames from the trajectory
        print(f"Loading trajectory from {trajectory_file}...")
        try:
            self.structure = self.parser.get_structure('trajectory', trajectory_file)
            
            # Count the number of frames
            self.num_frames = len(self.structure)
            print(f"Found {self.num_frames} frames in trajectory")
            
            # Determine which frames to analyze
            self.frames_to_analyze = frames_to_analyze if frames_to_analyze > 0 else self.num_frames
            if self.frames_to_analyze > self.num_frames:
                self.frames_to_analyze = self.num_frames
                print(f"Adjusting to analyze all {self.num_frames} frames")
            else:
                print(f"Will analyze {self.frames_to_analyze} frames")
        except Exception as e:
            print(f"Error loading trajectory: {e}")
            self.structure = None
            self.num_frames = 0
            self.frames_to_analyze = 0
        
    def analyze_frames(self, cdr_residues, ab_lock_residues, substrate_residues, 
                      stride=1, output_csv="trajectory_analysis.csv", probe_radius=40.0):
        """
        Analyze selected frames from the trajectory
        
        Args:
            cdr_residues: Dictionary of CDR regions and their residues
            ab_lock_residues: List of tuples (chain_id, residue_number)
            substrate_residues: List of tuples (chain_id, residue_number)
            stride: Analyze every Nth frame
            output_csv: Path to output CSV file
            probe_radius: Probe radius for PASA calculation in Angstroms
        
        Returns:
            Dictionary with analysis results
        """
        if self.structure is None:
            print("No valid structure to analyze")
            return None
            
        results = []
        frame_indices = range(0, self.frames_to_analyze, stride)
        total_frames = len(frame_indices)
        
        # Create analyzer for individual frames
        for i, frame_idx in enumerate(frame_indices):
            try:
                print(f"\nAnalyzing frame {frame_idx+1}/{self.frames_to_analyze} ({i+1}/{total_frames})...")
                model = self.structure[frame_idx]
                
                # Create a temporary PDB file for this frame
                temp_pdb = f"temp_frame_{frame_idx}.pdb"
                io = PDBIO()
                io.set_structure(model)
                io.save(temp_pdb)
                
                # Create frame analyzer
                frame_analyzer = FrameAnalyzer(temp_pdb, model)
                
                # Calculate metrics
                blocking_rates = frame_analyzer.calc_blocking_rate(cdr_residues, ab_lock_residues)
                pasa_results = frame_analyzer.calc_pasa(substrate_residues, probe_radius=probe_radius)
                
                # Collect results
                frame_result = {
                    'Frame': frame_idx,
                    'Time': frame_idx,  # Adjust if time step is known
                    **{f"BR_{cdr}": rate for cdr, rate in blocking_rates.items()},
                    **{f"PASA_{chain}": area for chain, area in pasa_results.items() if chain != 'total'},
                    'Total_PASA': pasa_results.get('total', 0)
                }
                results.append(frame_result)
                
                # Clean up temporary file
                if os.path.exists(temp_pdb):
                    os.remove(temp_pdb)
                
                # Print current results
                print("\nBlocking Rates:")
                for cdr in sorted([k for k in blocking_rates.keys() if k != 'total']):
                    print(f"{cdr:8}: {blocking_rates[cdr]:6.1f}%")
                print(f"{'Total':8}: {blocking_rates.get('total', 0):6.1f}%")
                
                print("\nPASA Values:")
                for chain in sorted([k for k in pasa_results.keys() if k != 'total']):
                    print(f"Chain {chain}: {pasa_results[chain]:6.2f} Å²")
                print(f"Total     : {pasa_results.get('total', 0):6.2f} Å²")
                
            except Exception as e:
                print(f"Error analyzing frame {frame_idx}: {str(e)}")
                import traceback
                traceback.print_exc()
                continue
        
        # Write results to CSV
        if results:
            headers = list(results[0].keys())
            with open(output_csv, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=headers)
                writer.writeheader()
                writer.writerows(results)
            
            print(f"\nAnalysis complete. Results saved to {output_csv}")
            
            # Calculate statistics
            stats_file = output_csv.replace('.csv', '_stats.csv')
            self.write_statistics(results, stats_file)
            
        else:
            print("No valid results to write to CSV")
            
        return results
    
    def write_statistics(self, results, output_file):
        """Calculate and write statistics to a file"""
        if not results:
            return
        
        # Initialize statistics
        stats = {
            'mean': {},
            'std': {},
            'min': {},
            'max': {}
        }
        
        # Collect all numeric columns
        numeric_columns = []
        for key in results[0].keys():
            if key not in ['Frame', 'Time']:
                numeric_columns.append(key)
        
        # Calculate statistics
        for col in numeric_columns:
            values = [result[col] for result in results]
            stats['mean'][col] = np.mean(values)
            stats['std'][col] = np.std(values)
            stats['min'][col] = np.min(values)
            stats['max'][col] = np.max(values)
        
        # Write to CSV
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow(['Metric'] + numeric_columns)
            
            # Write statistics
            for stat_name, stat_values in stats.items():
                row = [stat_name.capitalize()]
                for col in numeric_columns:
                    if 'BR_' in col or col == 'Total_BR':
                        row.append(f"{stat_values[col]:.1f}")
                    else:
                        row.append(f"{stat_values[col]:.2f}")
                writer.writerow(row)
        
        print(f"Statistics saved to {output_file}")
        
        return stats

def run_gbsa_simulation(input_file, output_prefix, simulation_steps=10000, 
                       report_interval=1000, gb_model='OBC2', ss_cutoff=0.25):
    """
    Run OpenMM GBSA simulation with a CIF file
    
    Args:
        input_file: Path to input CIF file
        output_prefix: Prefix for output files
        simulation_steps: Number of simulation steps
        report_interval: Interval for reporting progress
        gb_model: GB model to use ('OBC1', 'OBC2', or 'HCT')
        ss_cutoff: Distance cutoff for disulfide bond detection

    Returns:
        Path to trajectory PDB file or None if simulation failed
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_prefix)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Set paths for output files
    trajectory_file = f"{output_prefix}_trajectory.pdb"
    
    # Command to run OpenMM GBSA simulation
    cmd = [
        "python", "openmm-gbsa-simulation-final.py",
        "--cif", input_file,
        "--output", output_prefix,
        "--steps", str(simulation_steps),
        "--report", str(report_interval),
        "--gbmodel", gb_model,
        "--ss-cutoff", str(ss_cutoff)
    ]
    
    print(f"Running GBSA simulation for {input_file}...")
    print(f"Command: {' '.join(cmd)}")
    
    # Run the command
    try:
        subprocess.run(cmd, check=True)
        
        # Check if trajectory file was created
        if os.path.exists(trajectory_file):
            print(f"Simulation completed successfully. Trajectory saved to {trajectory_file}")
            return trajectory_file
        else:
            print(f"Simulation completed but trajectory file not found: {trajectory_file}")
            return None
    except subprocess.CalledProcessError as e:
        print(f"Error running simulation: {e}")
        return None
    except Exception as e:
        print(f"Unexpected error: {e}")
        return None

def analyze_trajectory(trajectory_file, config, output_file, frames=0, stride=1, probe_radius=40.0):
    """
    Analyze a trajectory PDB file
    
    Args:
        trajectory_file: Path to trajectory PDB file
        config: Configuration dictionary with CDR, lock, and substrate residues
        output_file: Path to output CSV file
        frames: Number of frames to analyze (0 for all)
        stride: Analyze every Nth frame
        probe_radius: Probe radius for PASA calculation
        
    Returns:
        Path to analysis CSV file or None if analysis failed
    """
    if not os.path.exists(trajectory_file):
        print(f"Trajectory file not found: {trajectory_file}")
        return None
    
    print(f"Analyzing trajectory file: {trajectory_file}")
    analyzer = TrajectoryAnalyzer(trajectory_file, frames)
    
    results = analyzer.analyze_frames(
        config['cdr_residues'],
        config['ab_lock_residues'],
        config['substrate_residues'],
        stride,
        output_file,
        probe_radius
    )
    
    if results:
        print(f"Analysis completed successfully. Results saved to {output_file}")
        return output_file
    else:
        print(f"Analysis failed or no valid frames found in {trajectory_file}")
        return None

def combine_analysis_files(analysis_files, output_file):
    """
    Combine multiple analysis CSV files into one
    
    Args:
        analysis_files: List of tuples (file_label, file_path)
        output_file: Path to combined output CSV file
        
    Returns:
        Path to combined CSV file or None if combination failed
    """
    if not analysis_files:
        print("No analysis files to combine")
        return None
    
    # Initialize empty DataFrame for combined results
    combined_df = pd.DataFrame()
    
    # Read each analysis file and append to combined DataFrame
    for label, file_path in analysis_files:
        if os.path.exists(file_path):
            try:
                # Read CSV file
                df = pd.read_csv(file_path)
                
                # Add a column for the source file
                df['Source'] = label
                
                # Append to combined DataFrame
                if combined_df.empty:
                    combined_df = df
                else:
                    combined_df = pd.concat([combined_df, df], ignore_index=True)
                    
                print(f"Added data from {file_path}")
            except Exception as e:
                print(f"Error reading {file_path}: {e}")
        else:
            print(f"File not found: {file_path}")
    
    # Write combined DataFrame to CSV
    if not combined_df.empty:
        combined_df.to_csv(output_file, index=False)
        print(f"Combined results saved to {output_file}")
        return output_file
    else:
        print("No valid data to combine")
        return None

def combine_statistics(analysis_files, output_file):
    """
    Combine statistics from multiple analyses into one summary file
    
    Args:
        analysis_files: List of tuples (file_label, file_path)
        output_file: Path to combined statistics output CSV file
        
    Returns:
        Path to combined statistics file or None if combination failed
    """
    if not analysis_files:
        print("No analysis files to combine statistics from")
        return None
    
    # Initialize dictionary for combined statistics
    combined_stats = {}
    
    # Read statistics for each analysis file
    for label, file_path in analysis_files:
        # Convert to stats file path
        stats_file = file_path.replace('.csv', '_stats.csv')
        
        if os.path.exists(stats_file):
            try:
                # Read statistics CSV file
                df = pd.read_csv(stats_file)
                
                # Extract mean values and add to combined statistics
                try:
                    mean_row = df[df.iloc[:, 0] == 'Mean']
                    if not mean_row.empty:
                        # Convert mean values to a dictionary
                        mean_dict = {}
                        for col in df.columns[1:]:
                            mean_dict[col] = float(mean_row[col].values[0])
                        
                        # Add to combined statistics
                        combined_stats[label] = mean_dict
                        print(f"Added statistics from {stats_file}")
                except Exception as e:
                    print(f"Error extracting mean values from {stats_file}: {e}")
            except Exception as e:
                print(f"Error reading statistics file {stats_file}: {e}")
        else:
            print(f"Statistics file not found: {stats_file}")
    
    # Convert combined statistics to DataFrame
    if combined_stats:
        # Create a list of all metrics from all files
        all_metrics = set()
        for label, stats in combined_stats.items():
            all_metrics.update(stats.keys())
        
        # Initialize DataFrame with all metrics
        summary_df = pd.DataFrame(index=combined_stats.keys(), columns=sorted(all_metrics))
        
        # Fill in values
        for label, stats in combined_stats.items():
            for metric, value in stats.items():
                summary_df.loc[label, metric] = value
        
        # Format values
        for col in summary_df.columns:
            if 'BR_' in col or col == 'Total_BR':
                summary_df[col] = summary_df[col].map(lambda x: f"{x:.1f}" if pd.notnull(x) else "")
            else:
                summary_df[col] = summary_df[col].map(lambda x: f"{x:.2f}" if pd.notnull(x) else "")
        
        # Reset index and rename index column
        summary_df.reset_index(inplace=True)
        summary_df.rename(columns={'index': 'Structure'}, inplace=True)
        
        # Write to CSV
        summary_df.to_csv(output_file, index=False)
        print(f"Combined statistics summary saved to {output_file}")
        return output_file
    else:
        print("No valid statistics to combine")
        return None

def main():
    """Main workflow for integrated GBSA simulation and analysis"""
    parser = argparse.ArgumentParser(description='Integrated GBSA Simulation and Analysis Tool')
    parser.add_argument('-c', '--config', default='input_config.txt',
                      help='Configuration file')
    parser.add_argument('-o', '--output-dir', default='results',
                      help='Output directory for all results')
    parser.add_argument('-s', '--steps', type=int, default=10000,
                      help='Number of simulation steps (default: 10000)')
    parser.add_argument('-r', '--report', type=int, default=1000,
                      help='Reporting interval for simulation (default: 1000)')
    parser.add_argument('-g', '--gbmodel', default='OBC2', choices=['HCT', 'OBC1', 'OBC2'],
                      help='GB model to use (default: OBC2)')
    parser.add_argument('--ss-cutoff', type=float, default=0.25,
                      help='Distance cutoff for disulfide bond detection (default: 0.25 nm)')
    parser.add_argument('-f', '--frames', type=int, default=0,
                      help='Number of frames to analyze (0 for all) (default: 0)')
    parser.add_argument('--stride', type=int, default=1,
                      help='Analyze every Nth frame (default: 1)')
    parser.add_argument('-p', '--probe', type=float, default=40.0,
                      help='Probe radius for PASA calculation in Angstroms (default: 40.0)')
    parser.add_argument('--skip-sim', action='store_true',
                      help='Skip simulation and only perform analysis')
    parser.add_argument('--skip-analysis', action='store_true',
                      help='Skip analysis and only perform simulation')
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    try:
        # Read configuration
        print(f"Reading configuration from {args.config}...")
        config = read_config(args.config)
        
        if not config['input_files']:
            raise ValueError("No input files specified in config")
        
        # Process each input file
        trajectory_files = []
        analysis_files = []
        
        for i, input_file in enumerate(config['input_files']):
            if not os.path.exists(input_file):
                print(f"Input file not found: {input_file}")
                continue
            
            # Get base name for output files
            basename = os.path.splitext(os.path.basename(input_file))[0]
            sim_output_prefix = os.path.join(args.output_dir, basename)
            trajectory_file = f"{sim_output_prefix}_trajectory.pdb"
            analysis_file = f"{sim_output_prefix}_analysis.csv"
            
            print(f"\n{'='*80}\nProcessing {i+1}/{len(config['input_files'])}: {input_file}\n{'='*80}")
            
            # Run GBSA simulation if not skipped
            if not args.skip_sim:
                # Check if trajectory already exists
                if os.path.exists(trajectory_file):
                    print(f"Trajectory file already exists: {trajectory_file}")
                    overwrite = input("Do you want to rerun the simulation? (y/n): ").lower().strip() == 'y'
                    if not overwrite:
                        print("Skipping simulation, using existing trajectory file")
                        trajectory_files.append((basename, trajectory_file))
                        continue
                
                # Run simulation
                result_trajectory = run_gbsa_simulation(
                    input_file,
                    sim_output_prefix,
                    args.steps,
                    args.report,
                    args.gbmodel,
                    args.ss_cutoff
                )
                
                if result_trajectory:
                    trajectory_files.append((basename, result_trajectory))
                else:
                    print(f"Simulation failed for {input_file}, skipping analysis")
                    continue
            else:
                # If simulation is skipped, check if trajectory exists
                if os.path.exists(trajectory_file):
                    print(f"Using existing trajectory file: {trajectory_file}")
                    trajectory_files.append((basename, trajectory_file))
                else:
                    print(f"Trajectory file not found: {trajectory_file}")
                    print(f"Cannot analyze {input_file} without trajectory file")
                    continue
            
            # Run trajectory analysis if not skipped
            if not args.skip_analysis:
                # Check if analysis already exists
                if os.path.exists(analysis_file):
                    print(f"Analysis file already exists: {analysis_file}")
                    overwrite = input("Do you want to rerun the analysis? (y/n): ").lower().strip() == 'y'
                    if not overwrite:
                        print("Skipping analysis, using existing analysis file")
                        analysis_files.append((basename, analysis_file))
                        continue
                
                # Run analysis
                result_analysis = analyze_trajectory(
                    trajectory_files[-1][1],  # Use the latest trajectory file
                    config,
                    analysis_file,
                    args.frames,
                    args.stride,
                    args.probe
                )
                
                if result_analysis:
                    analysis_files.append((basename, result_analysis))
            else:
                # If analysis is skipped, check if analysis file exists
                if os.path.exists(analysis_file):
                    print(f"Using existing analysis file: {analysis_file}")
                    analysis_files.append((basename, analysis_file))
        
        # Combine all analysis results
        if analysis_files:
            print(f"\n{'='*80}\nCombining analysis results from {len(analysis_files)} files\n{'='*80}")
            
            # Combine all analysis CSV files
            combined_csv = os.path.join(args.output_dir, "combined_analysis.csv")
            combined_result = combine_analysis_files(analysis_files, combined_csv)
            
            # Combine statistics summary
            combined_stats = os.path.join(args.output_dir, "combined_statistics.csv")
            stats_result = combine_statistics(analysis_files, combined_stats)
            
            if combined_result and stats_result:
                print(f"\nAll processing complete!")
                print(f"Combined analysis results: {combined_csv}")
                print(f"Combined statistics summary: {combined_stats}")
        else:
            print("\nNo analysis files available to combine")
    
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()