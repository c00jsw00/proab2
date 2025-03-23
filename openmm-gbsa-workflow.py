def compare_pre_post_md_analysis(pre_results, post_results, output_prefix='output'):
    """
    比较MD模拟前后的结构分析结果，生成差异报告
    
    Parameters
    ----------
    pre_results : tuple
        MD前的阻断率和PASA结果元组(blocking_rates, pasa_results)
    post_results : tuple
        MD后的阻断率和PASA结果元组(blocking_rates, pasa_results)
    output_prefix : str
        输出文件前缀
    """
    if not pre_results or not post_results:
        print("无法生成比较结果：缺少MD前或MD后的分析数据")
        return
    
    pre_br, pre_pasa = pre_results
    post_br, post_pasa = post_results
    
    # 准备比较数据
    headers = ['指标', 'MD前', 'MD后', '变化', '变化百分比(%)']
    rows = []
    
    # 比较阻断率
    for cdr in sorted(pre_br.keys()):
        if cdr in post_br:
            pre_value = pre_br[cdr]
            post_value = post_br[cdr]
            change = post_value - pre_value
            percent_change = (change / pre_value * 100) if pre_value != 0 else float('inf')
            
            rows.append([
                f"阻断率-{cdr}", 
                f"{pre_value:.1f}%", 
                f"{post_value:.1f}%", 
                f"{change:+.1f}%", 
                f"{percent_change:+.1f}" if percent_change != float('inf') else "N/A"
            ])
    
    # 比较PASA值
    for chain in sorted(pre_pasa.keys()):
        if chain in post_pasa:
            pre_value = pre_pasa[chain]
            post_value = post_pasa[chain]
            change = post_value - pre_value
            percent_change = (change / pre_value * 100) if pre_value != 0 else float('inf')
            
            chain_label = "总PASA" if chain == 'total' else f"PASA-{chain}"
            rows.append([
                chain_label, 
                f"{pre_value:.2f} Å²", 
                f"{post_value:.2f} Å²", 
                f"{change:+.2f} Å²", 
                f"{percent_change:+.1f}" if percent_change != float('inf') else "N/A"
            ])
    
    # 写入CSV
    output_file = f'{output_prefix}_md_effect_comparison.csv'
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(rows)
    
    print(f"MD前后比较分析完成。结果已保存到{output_file}")
    
    # 打印主要变化
    print("\nMD模拟导致的主要结构变化:")
    for row in rows:
        # 只打印变化超过5%的指标
        if row[4] != "N/A" and abs(float(row[4].strip('%+'))) > 5:
            print(f"{row[0]}: {row[1]} → {row[2]} ({row[3]}, {row[4]}%)")
    
    return rowsdef 

def analyze_after_md(pdb_file, config, trajectory_file=None, output_prefix='output'):
    """
    在MD模拟后分析蛋白质结构
    
    Parameters
    ----------
    pdb_file : str
        拓扑文件路径（PDB格式）
    config : dict
        包含CDR残基和其他分析参数的配置字典
    trajectory_file : str, optional
        轨迹文件路径（DCD格式），如果提供则分析轨迹中的最后一帧
    output_prefix : str
        输出文件前缀
    """
    print("\n开始MD后结构分析...")
    
    # 如果提供了轨迹文件，提取最后一帧
    if trajectory_file and os.path.exists(trajectory_file):
        print(f"从轨迹提取最后一帧: {trajectory_file}")
        final_structure = f"{output_prefix}_final_frame.pdb"
        
        try:
            # 使用MDTraj加载轨迹并保存最后一帧
            traj = md.load(trajectory_file, top=pdb_file)
            last_frame = traj[-1]
            last_frame.save(final_structure)
            print(f"最后一帧已保存到: {final_structure}")
            
            # 使用最后一帧进行结构分析
            pdb_file = final_structure
        except Exception as e:
            print(f"提取轨迹最后一帧时出错: {str(e)}")
            print("将使用初始结构进行分析")
    
    # 分析结构
    print(f"分析结构: {pdb_file}")
    blocking_rates, pasa_results = analyze_protein_structure(pdb_file, config)
    
    if blocking_rates and pasa_results:
        # 收集数据
        headers = ['PDB', 'CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3', 'HBCDR', 'Total_BR', 
                  'PASA_H', 'PASA_L', 'Total_PASA']
        
        row = [os.path.basename(pdb_file)]
        for cdr in ['CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3', 'HBCDR']:
            row.append(f"{blocking_rates.get(cdr, 0):.1f}")
        row.append(f"{blocking_rates.get('total', 0):.1f}")
        row.append(f"{pasa_results.get('H', 0):.2f}")
        row.append(f"{pasa_results.get('L', 0):.2f}")
        row.append(f"{pasa_results.get('total', 0):.2f}")
        
        # 写入CSV
        output_file = f'{output_prefix}_post_md_analysis.csv'
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerow(row)
        
        print(f"MD后结构分析完成。结果已保存到{output_file}")
        
        # 返回结果供后续使用
        return blocking_rates, pasa_results
    
    return None, None#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
使用OpenMM对含有二硫键(SSBOND)的AlphaFold生成的蛋白质结构进行GBSA隐式溶剂分子动力学模拟
"""

import os
import sys
import numpy as np
import argparse
import csv
from collections import defaultdict
from datetime import datetime
from openmm import unit
from openmm import app
import openmm as mm
from openmm.app import PDBFile, PDBxFile, ForceField, Simulation, StateDataReporter
from openmm.app import Modeller, OBC2, HBonds, AllBonds
from openmm.app import CutoffNonPeriodic, PME
import mdtraj as md
import Bio.PDB
from Bio.PDB import PDBIO, MMCIFParser
import freesasa

def fix_structure_and_add_ssbonds(input_file, output_pdb):
    """
    读取AlphaFold生成的CIF文件，修正结构并添加二硫键连接
    
    Parameters
    ----------
    input_file : str
        输入的AlphaFold生成的CIF文件路径
    output_pdb : str
        输出的修正后PDB文件路径
    
    Returns
    -------
    list
        二硫键残基对列表，每个元素为(chain_id1, resid1, chain_id2, resid2)
    """
    # 使用MDTraj读取CIF文件
    traj = md.load(input_file)
    
    # 查找潜在的二硫键
    # 在蛋白质中，二硫键通常由两个半胱氨酸(CYS)残基之间形成
    # 当两个CYS残基的硫原子距离小于约2.05埃时，它们可能形成二硫键
    
    # 获取所有CYS残基
    cys_indices = [atom.index for atom in traj.topology.atoms if atom.name == 'SG' and atom.residue.name == 'CYS']
    
    # 如果少于2个CYS残基，则无法形成二硫键
    if len(cys_indices) < 2:
        print("警告：结构中没有足够的CYS残基来形成二硫键")
        traj.save(output_pdb)
        return []
    
    # 计算所有CYS残基的SG原子之间的距离
    distances = md.compute_distances(traj, [[i, j] for i in cys_indices for j in cys_indices if i < j])[0]
    
    # 获取可能形成二硫键的CYS残基对
    potential_bonds = []
    pair_indices = [[i, j] for i in cys_indices for j in cys_indices if i < j]
    
    for i, (atom_i, atom_j) in enumerate(pair_indices):
        if distances[i] < 0.25:  # 0.25 nm = 2.5 Å，二硫键的典型距离
            res_i = traj.topology.atom(atom_i).residue
            res_j = traj.topology.atom(atom_j).residue
            print(f"发现潜在的二硫键: {res_i.chain.index+1}:{res_i.name}{res_i.resSeq} - {res_j.chain.index+1}:{res_j.name}{res_j.resSeq}, 距离: {distances[i]*10:.2f}Å")
            potential_bonds.append((res_i.chain.index, res_i.resSeq, res_j.chain.index, res_j.resSeq))
    
    # 保存修正后的结构为PDB格式
    traj.save(output_pdb)
    
    return potential_bonds


def prepare_structure_for_md(input_file, output_pdb, disulfide_bonds=None):
    """
    准备结构用于分子动力学模拟，修复常见问题
    
    Parameters
    ----------
    input_file : str
        输入的结构文件路径（CIF或PDB格式）
    output_pdb : str
        输出的修正后PDB文件路径
    disulfide_bonds : list, optional
        二硫键残基对列表
        
    Returns
    -------
    str
        准备好的结构文件路径
    """
    from Bio.PDB import PDBParser, PDBIO, Select
    from Bio.PDB.Polypeptide import is_aa
    import os
    import tempfile
    import subprocess
    
    print(f"准备结构用于MD模拟: {input_file}")
    
    # 1. 使用BioPython读取并清理结构
    temp_pdb = tempfile.NamedTemporaryFile(suffix='.pdb', delete=False).name
    
    if input_file.endswith('.cif'):
        # 如果是CIF文件，先转换为PDB
        from Bio.PDB import MMCIFParser
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('structure', input_file)
    else:
        # 否则直接读取PDB
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('structure', input_file)
    
    # 创建一个选择器，只保留蛋白质部分
    class ProteinSelect(Select):
        def accept_residue(self, residue):
            # 只接受氨基酸残基
            if is_aa(residue):
                return 1
            else:
                return 0
    
    # 保存清理后的结构
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb, ProteinSelect())
    
    print(f"已清理结构并保存到临时文件: {temp_pdb}")
    
    # 2. 使用pdbfixer处理结构
    # 如果系统中有pdbfixer，则使用它进行进一步清理
    fixed_pdb = output_pdb
    try:
        from pdbfixer import PDBFixer
        from simtk.openmm.app import PDBFile
        
        print("使用PDBFixer修复结构...")
        fixer = PDBFixer(filename=temp_pdb)
        
        # 查找缺失残基
        fixer.findMissingResidues()
        
        # 添加缺失的原子
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        # 添加缺失的氢原子
        fixer.addMissingHydrogens(7.0)  # pH 7.0
        
        # 处理二硫键
        if disulfide_bonds:
            print("添加二硫键...")
            # 清除自动检测的二硫键
            fixer.disulfideBonds = []
            
            # 添加指定的二硫键
            for chain_id1, resid1, chain_id2, resid2 in disulfide_bonds:
                chain1 = str(chain_id1 + 1) if isinstance(chain_id1, int) else chain_id1
                chain2 = str(chain_id2 + 1) if isinstance(chain_id2, int) else chain_id2
                
                res1 = (chain1, resid1)
                res2 = (chain2, resid2)
                
                # 检查残基是否存在于结构中
                found1 = False
                found2 = False
                for chain in fixer.topology.chains():
                    if chain.id == chain1:
                        for residue in chain.residues():
                            if residue.id == str(resid1):
                                found1 = True
                    if chain.id == chain2:
                        for residue in chain.residues():
                            if residue.id == str(resid2):
                                found2 = True
                
                if found1 and found2:
                    fixer.disulfideBonds.append((res1, res2))
                    print(f"添加二硫键: Chain {chain1} Res {resid1} - Chain {chain2} Res {resid2}")
                else:
                    print(f"警告: 未找到用于二硫键的残基 Chain {chain1} Res {resid1} 或 Chain {chain2} Res {resid2}")
        
        # 保存修复后的结构
        with open(fixed_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
        
        print(f"结构已修复并保存到: {fixed_pdb}")
    
    except ImportError:
        print("PDBFixer未安装，尝试使用替代方法...")
        # 如果没有pdbfixer，使用OpenMM的Modeller直接添加氢原子
        try:
            from simtk.openmm.app import PDBFile, Modeller, ForceField
            
            print("使用OpenMM Modeller添加氢原子...")
            pdb = PDBFile(temp_pdb)
            modeller = Modeller(pdb.topology, pdb.positions)
            
            # 尝试添加氢原子，但不使用力场（避免模板匹配问题）
            modeller.addHydrogens()
            
            # 保存修复后的结构
            with open(fixed_pdb, 'w') as f:
                PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
            print(f"已添加氢原子并保存到: {fixed_pdb}")
        
        except Exception as e:
            print(f"使用OpenMM Modeller时出错: {e}")
            print("将使用原始结构进行后续处理...")
            # 如果都失败了，至少保留清理后的结构
            import shutil
            shutil.copy(temp_pdb, fixed_pdb)
    
    # 清理临时文件
    try:
        os.remove(temp_pdb)
    except:
        pass
    
    return fixed_pdb

def run_gbsa_md(pdb_file, disulfide_bonds, output_prefix, temperature=300, sim_time=10, report_interval=10000, 
                dcd_interval=10000, min_steps=5000, equil_steps=100000):
    """
    使用OpenMM进行GBSA隐式溶剂的分子动力学模拟
    
    Parameters
    ----------
    pdb_file : str
        输入的PDB文件路径
    disulfide_bonds : list
        二硫键残基对列表，每个元素为(chain_id1, resid1, chain_id2, resid2)
    output_prefix : str
        输出文件的前缀
    temperature : float, optional
        模拟温度(K)，默认为300K
    sim_time : float, optional
        模拟时间(ns)，默认为10ns
    report_interval : int, optional
        报告间隔(步数)，默认为10000步
    dcd_interval : int, optional
        轨迹保存间隔(步数)，默认为10000步
    min_steps : int, optional
        能量最小化步数，默认为5000步
    equil_steps : int, optional
        平衡步数，默认为100000步
    """
    # 设置单位
    from openmm import unit
    from openmm import app
    import openmm as mm
    
    temperature = temperature * unit.kelvin
    step_size = 2.0 * unit.femtoseconds
    sim_time = sim_time * unit.nanoseconds
    
    # 计算总步数
    total_steps = int(sim_time / step_size)
    
    # 首先准备结构
    print(f"准备输入结构: {pdb_file}")
    prepared_pdb = f"{output_prefix}_prepared.pdb"
    prepared_pdb = prepare_structure_for_md(pdb_file, prepared_pdb, disulfide_bonds)
    
    try:
        # 读取PDB结构
        pdb = app.PDBFile(prepared_pdb)
        
        # 创建力场
        print("创建力场...")
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        
        # 创建模型系统，但不再添加氢原子（已在准备阶段添加）
        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        # 创建系统
        print("创建系统...")
        system = forcefield.createSystem(
            modeller.topology, 
            nonbondedMethod=app.CutoffNonPeriodic,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=app.HBonds,
            implicitSolvent=app.OBC2,
            implicitSolventSaltConc=0.15*unit.mole/unit.liter
        )
        
        # 添加Langevin积分器
        integrator = mm.LangevinMiddleIntegrator(
            temperature, 
            1.0/unit.picosecond, 
            step_size
        )
        
        # 创建模拟对象
        simulation = app.Simulation(modeller.topology, system, integrator)
        simulation.context.setPositions(modeller.positions)
        
        # 设置报告器
        simulation.reporters.append(
            app.StateDataReporter(
                f'{output_prefix}_data.csv', 
                report_interval, 
                step=True, 
                time=True, 
                potentialEnergy=True, 
                kineticEnergy=True, 
                totalEnergy=True, 
                temperature=True
            )
        )
        
        # 设置轨迹报告器
        simulation.reporters.append(
            app.DCDReporter(
                f'{output_prefix}_trajectory.dcd', 
                dcd_interval
            )
        )
        
        # 保存初始结构
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(f'{output_prefix}_initial.pdb', 'w'))
        
        # 能量最小化
        print("执行能量最小化...")
        simulation.minimizeEnergy(maxIterations=min_steps)
        
        # 保存最小化后的结构
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(f'{output_prefix}_minimized.pdb', 'w'))
        
        # 平衡模拟
        print(f"执行平衡模拟（{equil_steps}步）...")
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(equil_steps)
        
        # 保存平衡后的结构
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(f'{output_prefix}_equilibrated.pdb', 'w'))
        
        # 生产模拟
        print(f"执行生产模拟（{total_steps}步，{sim_time.value_in_unit(unit.nanoseconds):.1f} ns）...")
        simulation.step(total_steps)
        
        # 保存最终结构
        positions = simulation.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(simulation.topology, positions, open(f'{output_prefix}_final.pdb', 'w'))
        
        print("模拟完成！")
        return True
    
    except Exception as e:
        print(f"运行模拟时出错: {str(e)}")
        import traceback
        traceback.print_exc()
        return False 

def analyze_trajectory(trajectory_file, topology_file, output_prefix):
    """
    分析MD轨迹，计算RMSD、RMSF等
    
    Parameters
    ----------
    trajectory_file : str
        轨迹文件路径（DCD格式）
    topology_file : str
        拓扑文件路径（PDB格式）
    output_prefix : str
        输出文件前缀
    """
    # 加载轨迹
    traj = md.load(trajectory_file, top=topology_file)
    
    # 保存为XTC格式（可选，更紧凑的格式）
    traj.save(f'{output_prefix}_trajectory.xtc')
    
    # 计算蛋白质骨架的RMSD
    print("计算蛋白质骨架RMSD...")
    backbone = traj.topology.select('backbone')
    rmsd = md.rmsd(traj, traj, 0, atom_indices=backbone)
    
    # 保存RMSD数据
    np.savetxt(f'{output_prefix}_rmsd.dat', rmsd)
    
    # 计算每个残基的RMSF
    print("计算每个残基的RMSF...")
    ca_indices = traj.topology.select('name CA')
    rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_indices)
    
    # 保存RMSF数据
    np.savetxt(f'{output_prefix}_rmsf.dat', rmsf)
    
    # 计算二面角
    print("计算二面角...")
    phi_indices, phi_angles = md.compute_phi(traj)
    psi_indices, psi_angles = md.compute_psi(traj)
    
    # 保存二面角数据
    np.savetxt(f'{output_prefix}_phi.dat', phi_angles.reshape(phi_angles.shape[0], -1))
    np.savetxt(f'{output_prefix}_psi.dat', psi_angles.reshape(psi_angles.shape[0], -1))
    
    print("轨迹分析完成！")

# 添加蛋白质分析功能代码
def parse_residues(residue_str):
    """解析残基字符串，可以是范围(1-31)或单个残基(26, 30, 31)"""
    parts = residue_str.strip().split(',')
    residues = []
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = map(int, part.split('-'))
            residues.extend(range(start, end + 1))
        else:
            residues.append(int(part))
    return residues

def cif_to_pdb(cif_path, pdb_path):
    """使用Biopython将CIF文件转换为PDB格式"""
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('cif_structure', cif_path)
        io = PDBIO()
        io.set_structure(structure)
        io.save(pdb_path)
        print(f"已转换 {cif_path} 为 {pdb_path}")
        return True
    except Exception as e:
        print(f"转换失败 {cif_path}: {str(e)}")
        return False

def read_config(config_file):
    """读取配置文件，处理CIF/PDB文件和CDR残基信息"""
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
                if cdr_name.strip() == 'HBCDR':
                    # 特殊处理HBCDR格式
                    residue_pairs = []
                    parts = rest.strip().split(',')
                    for part in parts:
                        part = part.strip()
                        chain, residue = part.split()
                        residue_pairs.append((chain, int(residue)))
                    config['cdr_residues'][cdr_name.strip()] = residue_pairs
                else:
                    # 常规CDR格式，包含范围
                    chain, range_str = rest.strip().split()
                    residue_range = parse_residues(range_str)
                    config['cdr_residues'][cdr_name.strip()] = [(chain, i) for i in residue_range]
            elif current_section == 'ab_lock_residues':
                chain, range_str = line.split(':')
                chain = chain.strip()
                residue_range = parse_residues(range_str.strip())
                config['ab_lock_residues'].extend([(chain, i) for i in residue_range])
            elif current_section == 'substrate_residues':
                chain, range_str = line.split(':')
                chain = chain.strip()
                residue_range = parse_residues(range_str.strip())
                config['substrate_residues'].extend([(chain, i) for i in residue_range])
    
    # 处理文件转换
    processed_files = []
    for file_path in config['original_files']:
        if not file_path.endswith('.cif'):
            if os.path.exists(file_path):
                processed_files.append(file_path)
            else:
                print(f"未找到文件: {file_path}")
            continue
        
        pdb_path = file_path.replace('.cif', '.pdb')
        if os.path.exists(pdb_path):
            print(f"使用现有PDB文件: {pdb_path}")
            processed_files.append(pdb_path)
        else:
            if cif_to_pdb(file_path, pdb_path):
                processed_files.append(pdb_path)
            else:
                print(f"跳过无效文件: {file_path}")
    
    config['pdb_files'] = processed_files
    return config

class ProAntibodyAnalyzer:
    def __init__(self, pdb_file):
        """使用PDB文件初始化"""
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB文件未找到: {pdb_file}")
        
        self.parser = Bio.PDB.PDBParser(QUIET=True)
        self.structure = self.parser.get_structure('proab', pdb_file)
        self.model = self.structure[0]
        self.pdb_file = pdb_file
        
    def calc_blocking_rate(self, cdr_residues, ab_lock_residues):
        """
        基于CDR覆盖率计算阻断率
        Args:
            cdr_residues: CDR区域及其残基的字典
            ab_lock_residues: 元组(chain_id, residue_number)的列表
        Returns:
            每个CDR的阻断率字典和总阻断率
        """
        results = {}
        
        # 计算每个CDR区域的阻断率
        for cdr_name, residues in cdr_residues.items():
            covered = 0
            total = len(residues)
            
            for cdr_res in residues:
                # 获取CDR残基原子
                cdr_atoms = self.get_residue_atoms(cdr_res)
                
                # 检查是否被Ab lock覆盖
                is_covered = False
                for lock_res in ab_lock_residues:
                    lock_atoms = self.get_residue_atoms(lock_res)
                    
                    # 检查是否有任何原子对在4Å内且角度大于120度
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
            
            # 计算该CDR的阻断率
            rate = (covered/total) * 100 if total > 0 else 0
            results[cdr_name] = rate
            
        # 计算所有CDR阻断率的平均值
        cdr_rates = [rate for cdr_name, rate in results.items()]
        results['total'] = sum(cdr_rates) / len(cdr_rates) if cdr_rates else 0
        
        return results

    def calc_pasa(self, substrate_residues, probe_radius=40.0):
        """
        使用FreeSASA计算蛋白酶可及表面积
        Args:
            substrate_residues: 元组(chain_id, residue_number)的列表
            probe_radius: 探针半径(默认40Å，适用于MMP-2)
        Returns:
            各链和总PASA值的字典
        """
        # 按链分组残基
        chain_ranges = {}
        for chain_id, res_num in substrate_residues:
            if chain_id not in chain_ranges:
                chain_ranges[chain_id] = []
            chain_ranges[chain_id].append(res_num)

        results = {}
        
        try:
            # 从PDB文件创建结构
            structure = freesasa.Structure(self.pdb_file)
            result = freesasa.calc(structure)
            
            # 计算每条链残基范围的SASA
            for chain_id, residues in chain_ranges.items():
                residues.sort()  # 排序残基以便选择适当的范围
                # 使用Pymol语法创建选择表达式
                selection = f"s_{chain_id}, chain {chain_id} and resi {residues[0]}-{residues[-1]}"
                
                # 计算选定残基的面积
                areas = freesasa.selectArea([selection], structure, result)
                area = areas[f"s_{chain_id}"]
                results[chain_id] = area
            
            # 计算总PASA为链PASA的最小值
            results['total'] = min(results.values()) if results else 0
                
        except Exception as e:
            print(f"计算SASA时出错: {e}")
            return {'total': 0}
                
        return results

    def get_residue(self, res_id):
        """通过链ID和残基号获取残基"""
        chain_id, res_num = res_id
        for chain in self.model:
            if chain.id == chain_id:
                for residue in chain:
                    if residue.id[1] == res_num:
                        return residue
        return None

    def get_residue_atoms(self, res_id):
        """通过ID获取残基的原子"""
        residue = self.get_residue(res_id)
        return list(residue.get_atoms()) if residue else []

    def calc_distance(self, atom1, atom2):
        """计算两个原子之间的距离"""
        diff_vector = atom1.coord - atom2.coord
        return np.sqrt(np.sum(diff_vector * diff_vector))

    def calc_angle(self, atom1, atom2):
        """计算两个原子相对于垂直方向的角度"""
        vector = atom2.coord - atom1.coord
        vertical = np.array([0, 0, 1])
        cos_angle = np.dot(vector, vertical)/(np.linalg.norm(vector)*np.linalg.norm(vertical))
        return np.arccos(cos_angle) * 180/np.pi

def analyze_protein_structure(pdb_file, config):
    """分析蛋白质结构，计算阻断率和PASA值"""
    try:
        analyzer = ProAntibodyAnalyzer(pdb_file)
        
        # 计算阻断率
        blocking_rates = analyzer.calc_blocking_rate(config['cdr_residues'], config['ab_lock_residues'])
        
        # 计算PASA
        pasa_results = analyzer.calc_pasa(config['substrate_residues'])
        
        # 打印结果
        print("\n阻断率:")
        for cdr in sorted(blocking_rates.keys()):
            print(f"{cdr:8}: {blocking_rates[cdr]:6.1f}%")
        
        print("\nPASA值:")
        for chain in sorted(pasa_results.keys()):
            if chain != 'total':
                print(f"链 {chain}: {pasa_results[chain]:6.2f} Å²")
            else:
                print(f"总PASA: {pasa_results[chain]:6.2f} Å²")
                
        return blocking_rates, pasa_results
    
    except Exception as e:
        print(f"处理{pdb_file}时出错: {str(e)}")
        return None, None

def main():
    """主函数"""
    parser = argparse.ArgumentParser(description='AlphaFold蛋白质结构分析和GBSA MD模拟')
    parser.add_argument('--input', '-i', help='输入的AlphaFold生成的CIF文件')
    parser.add_argument('--config', '-c', help='配置文件（用于HBCDR分析）')
    parser.add_argument('--output', '-o', default='output', help='输出文件前缀')
    parser.add_argument('--temp', '-t', type=float, default=300.0, help='模拟温度(K)，默认300K')
    parser.add_argument('--time', type=float, default=10.0, help='模拟时间(ns)，默认10ns')
    parser.add_argument('--no-analyze', action='store_true', help='是否跳过轨迹分析')
    parser.add_argument('--only-structure-analysis', action='store_true', help='仅进行结构分析，不做MD模拟')
    parser.add_argument('--analyze-after-md', action='store_true', help='在MD模拟后进行结构分析')
    parser.add_argument('--create-example-config', action='store_true', help='创建示例配置文件并退出')
    args = parser.parse_args()
    
    # 创建示例配置文件
    if args.create_example_config:
        create_example_config()
        return
    
    # 创建输出目录
    os.makedirs(os.path.dirname(args.output) if os.path.dirname(args.output) else '.', exist_ok=True)
    
    # 初始化配置和MD前分析结果变量
    config = None
    pre_md_results = None
    
    # 如果提供了配置文件，则读取
    if args.config:
        print(f"读取配置文件: {args.config}")
        config = read_config(args.config)
        
        # 仅在不是MD后分析模式时执行批量分析
        if not args.analyze_after_md:
            # 验证文件
            valid_files = [f for f in config['pdb_files'] if os.path.exists(f)]
            if not valid_files:
                print("没有找到有效的结构文件")
                return
            
            print(f"找到{len(valid_files)}个有效的结构文件")
            
            # 初始化数据收集
            csv_data = []
            headers = ['PDB', 'CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3', 'HBCDR', 'Total_BR', 
                      'PASA_H', 'PASA_L', 'Total_PASA']
            
            # 分析每个结构
            for pdb_file in valid_files:
                print(f"\n分析{pdb_file}:")
                blocking_rates, pasa_results = analyze_protein_structure(pdb_file, config)
                
                if blocking_rates and pasa_results:
                    # 收集数据
                    row = [os.path.basename(pdb_file)]
                    for cdr in ['CDR-H1', 'CDR-H2', 'CDR-H3', 'CDR-L1', 'CDR-L2', 'CDR-L3', 'HBCDR']:
                        row.append(f"{blocking_rates.get(cdr, 0):.1f}")
                    row.append(f"{blocking_rates.get('total', 0):.1f}")
                    row.append(f"{pasa_results.get('H', 0):.2f}")
                    row.append(f"{pasa_results.get('L', 0):.2f}")
                    row.append(f"{pasa_results.get('total', 0):.2f}")
                    csv_data.append(row)
                    
                    # 如果这是将要用于MD模拟的文件，保存分析结果以便后续比较
                    if args.input and os.path.basename(pdb_file) == os.path.basename(args.input):
                        pre_md_results = (blocking_rates, pasa_results)
            
            # 计算统计数据
            if csv_data:
                avg_row = ['平均值']
                std_row = ['标准差']
                
                # 阻断率统计
                for cdr in headers[1:8]:  # 包括HBCDR
                    values = [float(row[headers.index(cdr)]) for row in csv_data]
                    avg_row.append(f"{np.mean(values):.1f}")
                    std_row.append(f"{np.std(values):.1f}")
                
                # 总阻断率统计
                values = [float(row[headers.index('Total_BR')]) for row in csv_data]
                avg_row.append(f"{np.mean(values):.1f}")
                std_row.append(f"{np.std(values):.1f}")
                
                # PASA统计
                for pasacol in ['PASA_H', 'PASA_L', 'Total_PASA']:
                    values = [float(row[headers.index(pasacol)]) for row in csv_data]
                    avg_row.append(f"{np.mean(values):.2f}")
                    std_row.append(f"{np.std(values):.2f}")
                
                csv_data.extend([avg_row, std_row])
                
                # 写入CSV
                output_file = f'{args.output}_structure_analysis.csv'
                with open(output_file, 'w', newline='') as f:
                    writer = csv.writer(f)
                    writer.writerow(headers)
                    writer.writerows(csv_data)
                
                print(f"\n结构分析完成。结果已保存到{output_file}")
    
    # MD模拟流程
    if args.input and not args.only_structure_analysis:
        # 修正结构并添加二硫键
        print(f"\n处理MD输入文件: {args.input}")
        fixed_pdb = f"{args.output}_fixed.pdb"
        disulfide_bonds = fix_structure_and_add_ssbonds(args.input, fixed_pdb)
        
        # 如果需要进行MD前结构分析且之前没有进行过
        if args.analyze_after_md and config and not pre_md_results:
            print("\n进行MD前结构分析...")
            pre_md_results = analyze_protein_structure(fixed_pdb, config)
        
        # 运行MD模拟
        print(f"\n开始MD模拟，温度: {args.temp}K，时间: {args.time}ns")
        run_gbsa_md(fixed_pdb, disulfide_bonds, args.output, 
                    temperature=args.temp, sim_time=args.time)
        
        # 常规轨迹分析
        if not args.no_analyze:
            print("\n开始轨迹分析...")
            analyze_trajectory(f'{args.output}_trajectory.dcd', 
                              f'{args.output}_initial.pdb', 
                              args.output)
        
        # MD后结构分析
        if args.analyze_after_md and config:
            print("\n开始MD后结构分析...")
            # 使用最终结构文件进行分析
            final_pdb = f"{args.output}_final.pdb"
            if not os.path.exists(final_pdb):
                print(f"警告: 未找到最终结构文件 {final_pdb}，将使用轨迹最后一帧")
                post_md_results = analyze_after_md(f'{args.output}_initial.pdb', config, 
                                                  f'{args.output}_trajectory.dcd', args.output)
            else:
                post_md_results = analyze_after_md(final_pdb, config, output_prefix=args.output)
            
            # 比较MD前后的结构变化
            if pre_md_results and post_md_results:
                print("\n比较MD前后的结构变化...")
                compare_pre_post_md_analysis(pre_md_results, post_md_results, args.output)
    
    # 仅运行MD后分析（用于已有的模拟结果）
    elif args.analyze_after_md and config:
        trajectory_file = f'{args.output}_trajectory.dcd'
        initial_pdb = f'{args.output}_initial.pdb'
        final_pdb = f'{args.output}_final.pdb'
        
        if not os.path.exists(trajectory_file) or not os.path.exists(initial_pdb):
            print(f"错误: 未找到必要的轨迹文件或初始结构文件")
            print(f"请确保轨迹文件({trajectory_file})和初始结构文件({initial_pdb})存在")
            return
        
        # 分析初始结构（MD前）
        print("\n分析MD前结构...")
        pre_md_results = analyze_protein_structure(initial_pdb, config)
        
        # 分析MD后结构
        print("\n分析MD后结构...")
        if os.path.exists(final_pdb):
            post_md_results = analyze_after_md(final_pdb, config, output_prefix=args.output)
        else:
            post_md_results = analyze_after_md(initial_pdb, config, trajectory_file, args.output)
        
        # 比较MD前后的结构变化
        if pre_md_results and post_md_results:
            print("\n比较MD前后的结构变化...")
            compare_pre_post_md_analysis(pre_md_results, post_md_results, args.output)
    
    elif not args.input and not args.analyze_after_md:
        parser.print_help()

def create_example_config(filename="example_config.txt"):
    """创建示例配置文件"""
    example_config = """# PDB/CIF文件
pdb_files:
pred.model_idx_0.cif
pred.model_idx_1.cif
pred.model_idx_2.cif
pred.model_idx_3.cif
pred.model_idx_4.cif

# CDR残基
cdr_residues:
CDR-H1: H 58-65
CDR-H2: H 83-89
CDR-H3: H 128-140
CDR-L1: L 58-63
CDR-L2: L 81-83
CDR-L3: L 120-128
HBCDR: H 86, H 90, H 134, H 135, H 136, L 122, L 123, L 125

# Ab锁残基
ab_lock_residues:
H: 1-32
L: 1-31

# 底物残基
substrate_residues:
H: 26-30
L: 21-26
"""
    with open(filename, 'w') as f:
        f.write(example_config)
    print(f"示例配置文件已创建: {filename}")

if __name__ == "__main__":
    start_time = datetime.now()
    main()
    end_time = datetime.now()
    print(f"总运行时间: {end_time - start_time}")
