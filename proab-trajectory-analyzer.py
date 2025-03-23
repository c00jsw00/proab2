#!/usr/bin/env python
# -*- coding: utf-8 -*-

import Bio.PDB
import numpy as np
import freesasa
import os
import csv
import argparse
from collections import defaultdict
from Bio.PDB import PDBIO, MMCIFParser, PDBParser

def parse_range(range_str):
    """Parse a range string like '1-31' into a list of numbers"""
    start, end = map(int, range_str.split('-'))
    return range(start, end + 1)

def read_config(config_file):
    """Read configuration file with trajectory processing info"""
    config = {
        'trajectory_files': [],
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
                config['trajectory_files'].append(line)
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
            Path to the generated CSV file
        """
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
            self.write_statistics(results, output_csv.replace('.csv', '_stats.csv'))
            
            return output_csv
        else:
            print("No valid results to write to CSV")
            return None
    
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

def main():
    """Main analysis workflow"""
    parser = argparse.ArgumentParser(description='MD Trajectory ProAntibody Analysis Tool')
    parser.add_argument('-c', '--config', default='sim_trajectory_config.txt',
                      help='Configuration file')
    parser.add_argument('-f', '--frames', type=int, default=0,
                      help='Number of frames to analyze (0 for all)')
    parser.add_argument('-s', '--stride', type=int, default=1,
                      help='Analyze every Nth frame')
    parser.add_argument('-o', '--output', default='trajectory_analysis.csv',
                      help='Output CSV file')
    parser.add_argument('-p', '--probe', type=float, default=40.0,
                      help='Probe radius for PASA calculation in Angstroms (default: 40.0)')
    args = parser.parse_args()
    
    try:
        print("Loading configuration...")
        config = read_config(args.config)
        
        if not config['trajectory_files']:
            raise ValueError("No trajectory files specified in config")
        
        # Analyze the first trajectory file
        trajectory_file = config['trajectory_files'][0]
        if not os.path.exists(trajectory_file):
            raise FileNotFoundError(f"Trajectory file not found: {trajectory_file}")
        
        # Create analyzer
        analyzer = TrajectoryAnalyzer(trajectory_file, args.frames)
        
        print(f"Using probe radius of {args.probe}Å for PASA calculations")
        
        # Analyze frames
        analyzer.analyze_frames(
            config['cdr_residues'],
            config['ab_lock_residues'],
            config['substrate_residues'],
            args.stride,
            args.output,
            args.probe
        )
    
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()
