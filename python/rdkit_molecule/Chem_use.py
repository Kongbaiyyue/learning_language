from calendar import c
import rdkit
from rdkit import Chem
from rdkit.Chem import rdChemReactions
# from pysmilesutils.augment import MolAugmenter

# 随机生成smiles
def random_smiles(smi_):
    mol_ = Chem.MolFromSmiles(smi_.replace(" ", "").replace("\n", ""))
    
    flag = 0
    while flag < 10:
        smi__ = Chem.MolToSmiles(mol_, doRandom=True)
        flag += 1
        if smi__ != smi_:
            flag = 10
    
    return smi__

# 示例用法
# smi = "C(CO)OC1CCN(C(OC(C)(C)C)=O)CC1"
# smi2 = "C(O)COC1CCN(C(=O)OC(C)(C)C)CC1"
# smi3 = Chem.MolToSmiles(Chem.MolFromSmiles(smi2), canonical=True)
# smi4 = Chem.MolToSmiles(Chem.MolFromSmiles(smi), canonical=True)
# aug = MolAugmenter()
# smi_list = aug([smi3])
# print(len(smi_list))
# mol_str = Chem.MolToSmiles(smi_list[0], canonical=False)
# print(mol_str)
# # print(smi3)
# # print(smi4)
# # print(random_smiles(smi))

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from rdkit import Chem

def generate_all_smiles(mol):
    """
    生成一个分子的所有可能的SMILES表示。
    
    参数：
    mol (Chem.Mol): 要生成SMILES表示的分子。
    
    返回：
    smiles_list (list): 包含所有可能的SMILES表示的列表。
    """
    smiles_list = []
    
    # 生成具有立体化学信息的SMILES表示
    smiles_isomeric = Chem.MolToSmiles(mol, isomericSmiles=True)
    smiles_list.append(smiles_isomeric)
    
    # 生成Kekulé形式的SMILES表示
    kekule_mol = Chem.MolFromSmiles(smiles_isomeric)
    if kekule_mol is not None:
        smiles_kekule = Chem.MolToSmiles(kekule_mol, kekuleSmiles=True)
        smiles_list.append(smiles_kekule)
    
    return smiles_list


from rdkit import Chem
from random import choice

# 随机断键
def delete_bond(mol):
    bonds = mol.GetBonds()
    random_bond = choice(bonds)
    mol.RemoveBond(random_bond.GetBeginAtomIdx(), random_bond.GetEndAtomIdx())
    smi = Chem.MolToSmiles(mol)
    print(smi)

# smiles = "C O C ( = O ) C C C ( = O ) c 1 c c c ( O C 2 C C C C O 2 ) c c 1 O"
# smiles = smiles.replace(" ", "").replace("\n", "")
# mol = Chem.MolFromSmiles(smiles)
# delete_bond(mol)
import pandas as pd
# df = pd.read_pickle("data/retro_uspto_50_template.pickle")
df = pd.read_csv("data/retro_uspto_50_template.csv")
