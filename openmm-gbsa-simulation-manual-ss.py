#!/usr/bin/env python

# 示例：如何手動定義二硫鍵位置並運行模擬

# 導入修改後的主程式
from openmm_gbsa_simulation_final import run_gbsa_simulation, find_disulfide_bonds, add_disulfide_bonds

# 定義已知的二硫鍵位置（如果您事先知道二硫鍵的位置）
KNOWN_SS_BONDS = [
    ('A', '22', 'A', '58'),  # 鏈A的第22個CYS與鏈A的第58個CYS之間的二硫鍵
    ('A', '76', 'A', '94'),  # 鏈A的第76個CYS與鏈A的第94個CYS之間的二硫鍵
    # 根據需要添加更多
]

def run_simulation_with_manual_ss_bonds(input_file, output_prefix, ss_bonds=None):
    """
    運行模擬並使用手動指定的二硫鍵
    
    Parameters:
    -----------
    input_file : str
        輸入結構檔案（CIF 或 PDB）
    output_prefix : str
        輸出檔案前綴
    ss_bonds : list of tuples
        手動指定的二硫鍵列表，每個元素為 (chain1, resid1, chain2, resid2)
    """
    # 使用主程式進行模擬，但覆蓋自動檢測到的二硫鍵
    # 將 ss_bonds 參數傳遞給修改後的 run_gbsa_simulation 函數
    
    print("運行模擬並使用手動指定的二硫鍵:")
    if ss_bonds:
        for bond in ss_bonds:
            print(f"  - 鏈 {bond[0]} CYS {bond[1]} - 鏈 {bond[2]} CYS {bond[3]}")
    
    # 呼叫主程式函數
    run_gbsa_simulation(
        input_file=input_file,
        output_prefix=output_prefix,
        simulation_steps=500000,
        report_interval=1000,
        gb_model='OBC2',
        ss_cutoff=0.3,  # 仍然使用自動檢測作為備用
        manual_ss_bonds=ss_bonds  # 傳遞手動指定的二硫鍵
    )

if __name__ == "__main__":
    # 使用示例
    structure_file = "pred.rank_0.cif"
    output_prefix = "results/manual_ss_simulation"
    
    # 使用已知的二硫鍵運行模擬
    run_simulation_with_manual_ss_bonds(
        input_file=structure_file,
        output_prefix=output_prefix,
        ss_bonds=KNOWN_SS_BONDS
    )
    
    # 或者使用自動檢測的二硫鍵
    # run_simulation_with_manual_ss_bonds(
    #     input_file=structure_file,
    #     output_prefix=output_prefix
    # )
