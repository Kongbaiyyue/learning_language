import rdkit
from rdkit import Chem
from rdkit.Chem import rdChemReactions

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
# smi = "C O c 1 c c c ( C N 2 C C c 3 n c ( N C ( = O ) N [C@H] ( C ) c 4 c c c c c 4 ) c c 4 c 3 c 2 n n 4 C ( c 2 c c c c c 2 ) ( c 2 c c c c c 2 ) c 2 c c c c c 2 ) c c 1"
# print(random_smiles(smi))

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

# # 示例用法
# smiles = "C O C ( = O ) C C C ( = O ) c 1 c c c ( O C 2 C C C C O 2 ) c c 1 O"
# smiles = smiles.replace(" ", "").replace("\n", "")
# mol = Chem.MolFromSmiles(smiles)

# if mol is not None:
#     all_smiles = generate_all_smiles(mol)
#     print(all_smiles)
# else:
#     print("无效的SMILES表示")

from rdkit import Chem
from random import choice

# 随机断键
def delete_bond(mol):
    bonds = mol.GetBonds()
    random_bond = choice(bonds)
    mol.RemoveBond(random_bond.GetBeginAtomIdx(), random_bond.GetEndAtomIdx())
    smi = Chem.MolToSmiles(mol)
    print(smi)

smiles = "C O C ( = O ) C C C ( = O ) c 1 c c c ( O C 2 C C C C O 2 ) c c 1 O"
smiles = smiles.replace(" ", "").replace("\n", "")
mol = Chem.MolFromSmiles(smiles)
delete_bond(mol)