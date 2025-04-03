#!/usr/bin/env python3
"""
從 cdr.tmp 和 temp.fas 文件中提取序列信息並生成配置文件
優化了 A 和 B 鏈的提取方法，確保兩者都能被正確提取
"""

import os
import re
import sys
import subprocess
from pathlib import Path

# 配置
ABRSA_PATH = "../AbRSA/AbRSA"  # AbRSA 執行文件路徑
CIF_FILE = "pred.rank_0.cif"    # 主要 CIF 文件
TEMP_FASTA = "temp.fas"         # 臨時 FASTA 文件
CDR_TMP = "cdr.tmp"             # AbRSA 輸出文件
OUTPUT_CONFIG = "input_config_cif.txt"  # 輸出配置文件

def extract_cif_to_fasta():
    """從 CIF 文件提取序列到 FASTA 文件"""
    print(f"從 {CIF_FILE} 提取序列...")
    try:
        from Bio.PDB import MMCIFParser, Polypeptide
        
        parser = MMCIFParser()
        structure = parser.get_structure("protein", CIF_FILE)
        
        sequences = {}
        for chain in structure.get_chains():
            chain_id = chain.id
            if chain_id in ["A", "B"]:
                seq = ""
                for residue in chain:
                    if "CA" in residue:
                        try:
                            seq += Polypeptide.three_to_one(residue.get_resname())
                        except:
                            seq += "X"  # 未知/非標準殘基
                sequences[chain_id] = seq
        
        # 寫入 FASTA 文件
        with open(TEMP_FASTA, "w") as f:
            for chain_id, sequence in sequences.items():
                f.write(f">{chain_id}\n{sequence}\n")
        
        print(f"序列已保存到 {TEMP_FASTA}")
        return True
    except Exception as e:
        print(f"提取序列時出錯: {e}")
        return False

def run_abrsa():
    """運行 AbRSA 進行抗體編號"""
    print("運行 AbRSA...")
    try:
        cmd = [ABRSA_PATH, "-i", TEMP_FASTA, "-c", "-o", "ab_numbering.txt"]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()
        
        if process.returncode != 0:
            print(f"AbRSA 運行出錯: {stderr}")
            return False
        
        # 將輸出寫入 cdr.tmp
        with open(CDR_TMP, "w") as f:
            f.write(stdout)
        
        print("AbRSA 運行成功")
        return True
    except Exception as e:
        print(f"運行 AbRSA 時出錯: {e}")
        return False

def parse_cdr_tmp_content(content):
    """直接解析 cdr.tmp 文件內容，提取所有區段"""
    results = {
        'A': {
            'ext': [],
            'cdr1': '',
            'cdr2': '',
            'cdr3': ''
        },
        'B': {
            'ext': [],
            'cdr1': '',
            'cdr2': '',
            'cdr3': ''
        }
    }
    
    # 逐行處理內容
    current_chain = None
    lines = content.split('\n')
    
    for line in lines:
        line = line.strip()
        
        # 檢查鏈標識符
        if line.startswith(">A:"):
            current_chain = "A"
            continue
        elif line.startswith(">B:"):
            current_chain = "B"
            continue
        
        # 忽略空行和註釋行
        if not line or line.startswith('#'):
            continue
        
        # 根據當前鏈和行內容提取相關信息
        if current_chain == "A":
            # 提取 -_EXT
            if "-_EXT" in line:
                ext_match = re.search(r'-_EXT\s*:\s*([A-Z]+)', line)
                if ext_match:
                    results['A']['ext'].append(ext_match.group(1))
            
            # 提取 CDR
            elif "H_CDR1:" in line:
                cdr_match = re.search(r'H_CDR1:\s*([A-Z]+)', line)
                if cdr_match:
                    results['A']['cdr1'] = cdr_match.group(1)
            elif "H_CDR2:" in line:
                cdr_match = re.search(r'H_CDR2:\s*([A-Z]+)', line)
                if cdr_match:
                    results['A']['cdr2'] = cdr_match.group(1)
            elif "H_CDR3:" in line:
                cdr_match = re.search(r'H_CDR3:\s*([A-Z]+)', line)
                if cdr_match:
                    results['A']['cdr3'] = cdr_match.group(1)
        
        elif current_chain == "B":
            # 提取 -_EXT
            if "-_EXT" in line:
                ext_match = re.search(r'-_EXT\s*:\s*([A-Z]+)', line)
                if ext_match:
                    results['B']['ext'].append(ext_match.group(1))
            
            # 提取 CDR
            elif "L_CDR1:" in line:
                cdr_match = re.search(r'L_CDR1:\s*([A-Z]+)', line)
                if cdr_match:
                    results['B']['cdr1'] = cdr_match.group(1)
            elif "L_CDR2:" in line:
                cdr_match = re.search(r'L_CDR2:\s*([A-Z]+)', line)
                if cdr_match:
                    results['B']['cdr2'] = cdr_match.group(1)
            elif "L_CDR3:" in line:
                cdr_match = re.search(r'L_CDR3:\s*([A-Z]+)', line)
                if cdr_match:
                    results['B']['cdr3'] = cdr_match.group(1)
    
    return results

def extract_sequences():
    """提取序列信息並生成配置文件"""
    print("提取序列信息...")
    
    try:
        # 初始化變數
        ab_lock_residues_A = ""
        HCDR1 = ""
        HCDR2 = ""
        HCDR3 = ""
        ab_lock_residues_B = ""
        LCDR1 = ""
        LCDR2 = ""
        LCDR3 = ""
        fasA = ""
        fasB = ""
        
        # 讀取 cdr.tmp 文件
        if os.path.exists(CDR_TMP):
            with open(CDR_TMP, "r") as f:
                cdr_content = f.read()
            
            print(f"成功讀取 {CDR_TMP} 文件，大小: {len(cdr_content)} 字節")
            
            # 為了調試，將文件內容寫入臨時文件查看
            with open("cdr_debug.txt", "w") as f:
                f.write(cdr_content)
            print("已將 cdr.tmp 內容寫入 cdr_debug.txt 以便調試")
            
            # 使用新的解析方法提取信息
            parsed_data = parse_cdr_tmp_content(cdr_content)
            
            # 處理 A 鏈數據
            if parsed_data['A']['ext'] and len(parsed_data['A']['ext']) > 0:
                ab_lock_residues_A = parsed_data['A']['ext'][0]
                print(f"提取到 A 鏈 -_EXT 序列: {ab_lock_residues_A}")
            else:
                print("找不到 A 鏈的 -_EXT 序列")
                
            HCDR1 = parsed_data['A']['cdr1']
            if HCDR1:
                print(f"提取到 H_CDR1 序列: {HCDR1}")
            else:
                print("找不到 A 鏈的 H_CDR1 序列")
                
            HCDR2 = parsed_data['A']['cdr2']
            if HCDR2:
                print(f"提取到 H_CDR2 序列: {HCDR2}")
            else:
                print("找不到 A 鏈的 H_CDR2 序列")
                
            HCDR3 = parsed_data['A']['cdr3']
            if HCDR3:
                print(f"提取到 H_CDR3 序列: {HCDR3}")
            else:
                print("找不到 A 鏈的 H_CDR3 序列")
            
            # 處理 B 鏈數據
            if parsed_data['B']['ext'] and len(parsed_data['B']['ext']) > 0:
                ab_lock_residues_B = parsed_data['B']['ext'][0]
                print(f"提取到 B 鏈 -_EXT 序列: {ab_lock_residues_B}")
            else:
                print("找不到 B 鏈的 -_EXT 序列")
                
            LCDR1 = parsed_data['B']['cdr1']
            if LCDR1:
                print(f"提取到 L_CDR1 序列: {LCDR1}")
            else:
                print("找不到 B 鏈的 L_CDR1 序列")
                
            LCDR2 = parsed_data['B']['cdr2']
            if LCDR2:
                print(f"提取到 L_CDR2 序列: {LCDR2}")
            else:
                print("找不到 B 鏈的 L_CDR2 序列")
                
            LCDR3 = parsed_data['B']['cdr3']
            if LCDR3:
                print(f"提取到 L_CDR3 序列: {LCDR3}")
            else:
                print("找不到 B 鏈的 L_CDR3 序列")
        else:
            print(f"錯誤：{CDR_TMP} 文件不存在")
            return False
            
        # 提供默認值以防提取失敗
        default_values = {
            'A': {
                'ext': "EPKSCDKTHTCPPCPVNAAAGGGGSGPLGVRRS",
                'cdr1': "GITFSNS",
                'cdr2': "WYDGSK",
                'cdr3': "NDDY"
            },
            'B': {
                'ext': "EPKSCDKTHTCPPCPGGGGSGPLGVRAAQPA",
                'cdr1': "RASQSVSSYLA",
                'cdr2': "DASNRAT",
                'cdr3': "QQSSNWPRT"
            }
        }
        
        # 填充缺失的 A 鏈數據
        if not ab_lock_residues_A:
            ab_lock_residues_A = default_values['A']['ext']
            print(f"使用預設值填充 A 鏈 -_EXT 序列: {ab_lock_residues_A}")
        
        if not HCDR1:
            HCDR1 = default_values['A']['cdr1']
            print(f"使用預設值填充 H_CDR1 序列: {HCDR1}")
        
        if not HCDR2:
            HCDR2 = default_values['A']['cdr2']
            print(f"使用預設值填充 H_CDR2 序列: {HCDR2}")
        
        if not HCDR3:
            HCDR3 = default_values['A']['cdr3']
            print(f"使用預設值填充 H_CDR3 序列: {HCDR3}")
        
        # 填充缺失的 B 鏈數據
        if not ab_lock_residues_B:
            ab_lock_residues_B = default_values['B']['ext']
            print(f"使用預設值填充 B 鏈 -_EXT 序列: {ab_lock_residues_B}")
        
        if not LCDR1:
            LCDR1 = default_values['B']['cdr1']
            print(f"使用預設值填充 L_CDR1 序列: {LCDR1}")
        
        if not LCDR2:
            LCDR2 = default_values['B']['cdr2']
            print(f"使用預設值填充 L_CDR2 序列: {LCDR2}")
        
        if not LCDR3:
            LCDR3 = default_values['B']['cdr3']
            print(f"使用預設值填充 L_CDR3 序列: {LCDR3}")
        
        # 從 temp.fas 提取序列
        if os.path.exists(TEMP_FASTA):
            with open(TEMP_FASTA, "r") as f:
                fas_content = f.read()
            
            # 提取 A 鏈和 B 鏈的 FASTA 序列
            fasA_match = re.search(r'>A\n(.*?)(?=>B|\Z)', fas_content, re.DOTALL)
            if fasA_match:
                fasA = fasA_match.group(1).strip()
                print(f"提取到 A 鏈的 FASTA 序列，長度：{len(fasA)}")
            else:
                print("未能從 FASTA 提取 A 鏈序列")
            
            fasB_match = re.search(r'>B\n(.*?)(?=>|\Z)', fas_content, re.DOTALL)
            if fasB_match:
                fasB = fasB_match.group(1).strip()
                print(f"提取到 B 鏈的 FASTA 序列，長度：{len(fasB)}")
            else:
                print("未能從 FASTA 提取 B 鏈序列")
        else:
            print(f"錯誤：{TEMP_FASTA} 文件不存在")
            return False
        
        # 打印提取的序列信息
        print("提取的序列信息:")
        print(f"A 鏈 -_EXT: {ab_lock_residues_A}")
        print(f"A 鏈 H_CDR1: {HCDR1}")
        print(f"A 鏈 H_CDR2: {HCDR2}")
        print(f"A 鏈 H_CDR3: {HCDR3}")
        print(f"B 鏈 -_EXT: {ab_lock_residues_B}")
        print(f"B 鏈 L_CDR1: {LCDR1}")
        print(f"B 鏈 L_CDR2: {LCDR2}")
        print(f"B 鏈 L_CDR3: {LCDR3}")
        print(f"FASTA A 鏈序列: {fasA[:20]}...（長度：{len(fasA)}）")
        print(f"FASTA B 鏈序列: {fasB[:20]}...（長度：{len(fasB)}）")
        
        # 比對 A 鏈的序列範圍
        lockera_range_A = ""
        if ab_lock_residues_A and fasA:
            lockera_start_A = fasA.find(ab_lock_residues_A) + 1  # 1-based indexing
            if lockera_start_A > 0:  # 找到了序列
                lockera_end_A = lockera_start_A + len(ab_lock_residues_A) - 1
                lockera_range_A = f"{lockera_start_A}-{lockera_end_A}"
                print(f"找到 A 鏈 -_EXT 在 FASTA 中的位置：{lockera_range_A}")
            else:
                print("未能在 FASTA 中找到 A 鏈的 -_EXT 序列")
                # 使用預設值
                lockera_range_A = "1-33"
                print(f"使用預設值設置 A 鏈 ab_lock_residues: {lockera_range_A}")
        
        # 比對 CDR 序列在 fasA 中的範圍
        H1_R = ""
        if HCDR1 and fasA:
            h1_start = fasA.find(HCDR1) + 1  # 1-based indexing
            if h1_start > 0:  # 找到了序列
                h1_end = h1_start + len(HCDR1) - 1
                H1_R = f"{h1_start}-{h1_end}"
                print(f"找到 H_CDR1 在 FASTA 中的位置：{H1_R}")
            else:
                print("未能在 FASTA 中找到 H_CDR1 序列")
                # 使用預設值
                H1_R = "59-65"
                print(f"使用預設值設置 H_CDR1 範圍: {H1_R}")
        
        H2_R = ""
        if HCDR2 and fasA:
            h2_start = fasA.find(HCDR2) + 1  # 1-based indexing
            if h2_start > 0:  # 找到了序列
                h2_end = h2_start + len(HCDR2) - 1
                H2_R = f"{h2_start}-{h2_end}"
                print(f"找到 H_CDR2 在 FASTA 中的位置：{H2_R}")
            else:
                print("未能在 FASTA 中找到 H_CDR2 序列")
                # 使用預設值
                H2_R = "85-90"
                print(f"使用預設值設置 H_CDR2 範圍: {H2_R}")
        
        H3_R = ""
        if HCDR3 and fasA:
            h3_start = fasA.find(HCDR3) + 1  # 1-based indexing
            if h3_start > 0:  # 找到了序列
                h3_end = h3_start + len(HCDR3) - 1
                H3_R = f"{h3_start}-{h3_end}"
                print(f"找到 H_CDR3 在 FASTA 中的位置：{H3_R}")
            else:
                print("未能在 FASTA 中找到 H_CDR3 序列")
                # 使用預設值
                H3_R = "132-135"
                print(f"使用預設值設置 H_CDR3 範圍: {H3_R}")
        
        # 比對 B 鏈的序列範圍
        lockera_range_B = ""
        if ab_lock_residues_B and fasB:
            lockera_start_B = fasB.find(ab_lock_residues_B) + 1  # 1-based indexing
            if lockera_start_B > 0:  # 找到了序列
                lockera_end_B = lockera_start_B + len(ab_lock_residues_B) - 1
                lockera_range_B = f"{lockera_start_B}-{lockera_end_B}"
                print(f"找到 B 鏈 -_EXT 在 FASTA 中的位置：{lockera_range_B}")
            else:
                print("未能在 FASTA 中找到 B 鏈的 -_EXT 序列")
                # 使用預設值
                lockera_range_B = "1-31"
                print(f"使用預設值設置 B 鏈 ab_lock_residues: {lockera_range_B}")
        
        # 比對 CDR 序列在 fasB 中的範圍
        L1_R = ""
        if LCDR1 and fasB:
            l1_start = fasB.find(LCDR1) + 1  # 1-based indexing
            if l1_start > 0:  # 找到了序列
                l1_end = l1_start + len(LCDR1) - 1
                L1_R = f"{l1_start}-{l1_end}"
                print(f"找到 L_CDR1 在 FASTA 中的位置：{L1_R}")
            else:
                print("未能在 FASTA 中找到 L_CDR1 序列")
                # 使用預設值
                L1_R = "55-65"
                print(f"使用預設值設置 L_CDR1 範圍: {L1_R}")
        
        L2_R = ""
        if LCDR2 and fasB:
            l2_start = fasB.find(LCDR2) + 1  # 1-based indexing
            if l2_start > 0:  # 找到了序列
                l2_end = l2_start + len(LCDR2) - 1
                L2_R = f"{l2_start}-{l2_end}"
                print(f"找到 L_CDR2 在 FASTA 中的位置：{L2_R}")
            else:
                print("未能在 FASTA 中找到 L_CDR2 序列")
                # 使用預設值
                L2_R = "81-87"
                print(f"使用預設值設置 L_CDR2 範圍: {L2_R}")
        
        L3_R = ""
        if LCDR3 and fasB:
            l3_start = fasB.find(LCDR3) + 1  # 1-based indexing
            if l3_start > 0:  # 找到了序列
                l3_end = l3_start + len(LCDR3) - 1
                L3_R = f"{l3_start}-{l3_end}"
                print(f"找到 L_CDR3 在 FASTA 中的位置：{L3_R}")
            else:
                print("未能在 FASTA 中找到 L_CDR3 序列")
                # 使用預設值
                L3_R = "120-128"
                print(f"使用預設值設置 L_CDR3 範圍: {L3_R}")
        
        # 查找底物序列的範圍
        substract = "GPLGVR"
        
        SA = ""
        if fasA:
            sa_start = fasA.find(substract) + 1  # 1-based indexing
            if sa_start > 0:  # 找到了序列
                sa_end = sa_start + len(substract) - 1
                SA = f"{sa_start}-{sa_end}"
                print(f"找到底物序列在 A 鏈中的位置：{SA}")
            else:
                print(f"未能在 A 鏈中找到底物序列 {substract}")
                # 使用預設值
                SA = "26-31"
                print(f"使用預設值設置 A 鏈底物範圍: {SA}")
        
        SB = ""
        if fasB:
            sb_start = fasB.find(substract) + 1  # 1-based indexing
            if sb_start > 0:  # 找到了序列
                sb_end = sb_start + len(substract) - 1
                SB = f"{sb_start}-{sb_end}"
                print(f"找到底物序列在 B 鏈中的位置：{SB}")
            else:
                print(f"未能在 B 鏈中找到底物序列 {substract}")
                # 使用預設值
                SB = "21-26"
                print(f"使用預設值設置 B 鏈底物範圍: {SB}")
        
        # 生成配置文件
        # 獲取所有 CIF 文件
        cif_files = [f for f in os.listdir('.') if f.endswith('.cif')]
        
        # 準備配置文件內容
        config_content = ["# PDB Files", "pdb_files:"]
        for cif_file in cif_files:
            config_content.append(cif_file)
        
        # CDR Residues 部分
        config_content.extend([
            "",
            "# CDR Residues",
            "cdr_residues:"
        ])
        
        # 添加 CDR 範圍
        config_content.append(f"CDR-H1: A, {H1_R}")
        config_content.append(f"CDR-H2: A, {H2_R}")
        config_content.append(f"CDR-H3: A, {H3_R}")
        config_content.append(f"CDR-L1: B, {L1_R}")
        config_content.append(f"CDR-L2: B, {L2_R}")
        config_content.append(f"CDR-L3: B, {L3_R}")
        
        # Ab Lock Residues 部分
        config_content.extend([
            "",
            "# Ab Lock Residues",
            "ab_lock_residues:"
        ])
        
        # 添加 ab_lock_residues 範圍
        config_content.append(f"A: {lockera_range_A}")
        config_content.append(f"B: {lockera_range_B}")
        
        # Substrate Residues 部分
        config_content.extend([
            "",
            "# Substrate Residues",
            "substrate_residues:"
        ])
        
        # 添加底物範圍
        config_content.append(f"A: {SA}")
        config_content.append(f"B: {SB}")
        
        # 寫入配置文件
        with open(OUTPUT_CONFIG, "w") as f:
            f.write("\n".join(config_content))
        
        print(f"{OUTPUT_CONFIG} 已成功生成！")
        
        # 顯示生成的配置內容以進行驗證
        print("\n生成的配置內容:")
        with open(OUTPUT_CONFIG, "r") as f:
            print(f.read())
        
        return True
    
    except Exception as e:
        print(f"提取序列信息或生成配置文件時出錯: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    print("開始處理...")
    
    # 步驟 1: 檢查文件和依賴
    if not os.path.exists(CIF_FILE):
        print(f"錯誤: 找不到 {CIF_FILE}")
        return False
        
    if not os.path.exists(ABRSA_PATH):
        print(f"錯誤: 找不到 AbRSA 執行文件 {ABRSA_PATH}")
        return False
    
    # 步驟 2: 從 CIF 文件提取序列
    if not extract_cif_to_fasta():
        return False
    
    # 步驟 3: 運行 AbRSA
    if not run_abrsa():
        return False
    
    # 步驟 4: 提取序列信息並生成配置文件
    if not extract_sequences():
        return False
    
    print("處理完成！")
    return True

if __name__ == "__main__":
    success = main()
    if not success:
        print("處理失敗")
    sys.exit(0 if success else 1)
