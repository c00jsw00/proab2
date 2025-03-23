import Bio.PDB
import numpy as np
import freesasa
import os
import csv
import argparse
from collections import defaultdict
from Bio.PDB import PDBIO, MMCIFParser

def parse_range(range_str):
    """Parse a range string like '1-31' into a list of numbers"""
    start, end = map(int, range_str.split('-'))
    return range(start, end + 1)

def cif_to_pdb(cif_path, pdb_path):
    """Convert CIF file to PDB format using Biopython"""
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('cif_structure', cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        print(f"Converted {cif_path} to {pdb_path}")
        return True
    except Exception as e:
        print(f"Failed to convert {cif_path}: {str(e)}")
        return False

def process_cif_files(file_list):
    """Process file list and convert CIF to PDB"""
    converted_files = []
    for file_path in file_list:
        if not file_path.endswith('.cif'):
            if os.path.exists(file_path):
                converted_files.append(file_path)
            else:
                print(f"File not found: {file_path}")
            continue
        
        pdb_path = file_path.replace('.cif', '.pdb')
        if os.path.exists(pdb_path):
            print(f"Using existing PDB file: {pdb_path}")
            converted_files.append(pdb_path)
        else:
            if cif_to_pdb(file_path, pdb_path):
                converted_files.append(pdb_path)
            else:
                print(f"Skipping invalid file: {file_path}")
    
    return converted_files

def read_config(config_file):
    """Read configuration file with CIF/PDB processing"""
    config = {
        'original_files': [],
        'pdb_files': [],
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
                config['original_files'].append(line)
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
    
    # Process file conversion
    config['pdb_files'] = process_cif_files(config['original_files'])
    return config

def create_example_config(filename):
    """Create example configuration file"""
    example_config = """# PDB/CIF Files
pdb_files:
input1.cif
input2.pdb

# CDR Residues
cdr_residues:
CDR-H1: H 58-65
CDR-H2: H 83-89
CDR-H3: H 128-140
CDR-L1: L 58-63
CDR-L2: L 81-83
CDR-L3: L 120-128

# Ab Lock Residues
ab_lock_residues:
H: 1-32
L: 1-31

# Substrate Residues
substrate_residues:
H: 26-30
L: 21-26
"""
    with open(filename, 'w') as f:
        f.write(example_config)
    print(f"Example config created: {filename}")

class ProAntibodyAnalyzer:
    def __init__(self, pdb_file):
        """Initialize with PDB file"""
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        self.parser = Bio.PDB.PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('proab', pdb_file)
        self.model = self.structure[0]
        self.pdb_file = pdb_file
        
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
                
                # Check if covered by Ab lock
                is_covered = False
                for lock_res in ab_lock_residues:
                    lock_atoms = self.get_residue_atoms(lock_res)
                    
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
            probe_radius: Radius of probe in Å (default 40Å for MMP-2)
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
            # Create structure from PDB file
            structure = freesasa.Structure(self.pdb_file)
            result = freesasa.calc(structure)
            
            # Calculate SASA for each chain's residue range
            for chain_id, residues in chain_ranges.items():
                residues.sort()  # Sort residues for proper range selection
                # Create selection expression using Pymol syntax
                selection = f"s_{chain_id}, chain {chain_id} and resi {residues[0]}-{residues[-1]}"
                
                # Calculate area for selected residues
                areas = freesasa.selectArea([selection], structure, result)
                area = areas[f"s_{chain_id}"]
                results[chain_id] = area
            
            # Calculate total PASA as minimum of chain PASAs
            results['total'] = min(results.values()) if results else 0
                
        except Exception as e:
            print(f"Error calculating SASA: {e}")
            return {'total': 0}
                
        return results

    def get_residue(self, res_id):
        """Get residue by chain ID and residue number"""
        chain_id, res_num = res_id
        for chain in self.model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] == res_num:
                        return residue
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
    parser = argparse.ArgumentParser(description='ProAntibody Analysis Tool')
    parser.add_argument('-c', '--config', default='input_config.txt',
                       help='Configuration file (default: input_config.txt)')
    parser.add_argument('--create-config', action='store_true',
                       help='Create example configuration file')
    args = parser.parse_args()
    
    if args.create_config:
        create_example_config("example_input_config.txt")
        return
    
    try:
        print("Loading configuration...")
        config = read_config(args.config)
        valid_files = [f for f in config['pdb_files'] if os.path.exists(f)]
        if not valid_files:
            raise ValueError("No valid structure files available")
        
        print(f"Found {len(valid_files)} valid structure files")
        
        # 初始化數據收集
        csv_data = []
        headers = ['PDB', 'CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3', 'Total_BR', 
                  'PASA_H', 'PASA_L', 'Total_PASA']
        
        all_blocking_rates = defaultdict(list)
        all_pasa_values = defaultdict(list)
        
        for pdb_file in valid_files:
            print(f"\nAnalyzing {pdb_file}:")
            try:
                analyzer = ProAntibodyAnalyzer(pdb_file)
                
                # 計算阻斷率
                blocking_rates = analyzer.calc_blocking_rate(config['cdr_residues'], config['ab_lock_residues'])
                
                # 計算PASA
                pasa_results = analyzer.calc_pasa(config['substrate_residues'])
                
                # 收集數據
                row = [os.path.basename(pdb_file)]
                for cdr in ['CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3']:
                    row.append(f"{blocking_rates.get(cdr, 0):.1f}")
                row.append(f"{blocking_rates.get('total', 0):.1f}")
                row.append(f"{pasa_results.get('A', 0):.2f}")
                row.append(f"{pasa_results.get('B', 0):.2f}")
                row.append(f"{pasa_results.get('total', 0):.2f}")
                csv_data.append(row)
                
                # 統計用數據
                for cdr, rate in blocking_rates.items():
                    all_blocking_rates[cdr].append(rate)
                for chain, area in pasa_results.items():
                    all_pasa_values[chain].append(area)
                
                # 打印即時結果
                print("\nBlocking Rates:")
                for cdr in sorted(blocking_rates.keys()):
                    print(f"{cdr:8}: {blocking_rates[cdr]:6.1f}%")
                
                print("\nPASA Values:")
                for chain in sorted(pasa_results.keys()):
                    print(f"Chain {chain}: {pasa_results[chain]:6.2f} Å²")
                    
            except Exception as e:
                print(f"Error processing {pdb_file}: {str(e)}")
                continue
        
        # 計算統計數據
        if csv_data:
            avg_row = ['Average']
            std_row = ['Std Dev']
            
            # 阻斷率統計
            for cdr in headers[1:7] + ['Total_BR']:
                key = cdr.replace('_BR', '') if cdr != 'Total_BR' else 'total'
                values = [float(row[headers.index(cdr)]) for row in csv_data]
                avg_row.append(f"{np.mean(values):.1f}")
                std_row.append(f"{np.std(values):.1f}")
            
            # PASA統計
            for pasacol in headers[8:11]:
                key = pasacol.split('_')[-1]
                values = [float(row[headers.index(pasacol)]) for row in csv_data]
                avg_row.append(f"{np.mean(values):.2f}")
                std_row.append(f"{np.std(values):.2f}")
            
            csv_data.extend([avg_row, std_row])
        
        # 寫入CSV
        output_file = 'proab_analysis_results.csv'
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(csv_data)
        
        print(f"\nAnalysis complete. Results saved to {output_file}")
    
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        if isinstance(e, FileNotFoundError):
            print("Please check:")
            print("1. Input file paths")
            print("2. Required packages (biopython, freesasa, numpy)")
            print("3. CIF file format compliance")

if __name__ == "__main__":
    main()