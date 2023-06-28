from math import isnan, nan
from rdkit import Chem
from rdkit.Chem import Draw

import re

from rdchiral.template_extractor import extract_from_reaction, get_changed_atoms, mols_from_smiles_list, \
    replace_deuterated, get_tetrahedral_atoms
from sklearn.metrics import recall_score

from molecule import get_graph_features_from_smi


# 去掉原子编号
def function_1(smi):
    # 将带有原子编号的SMILES转换为分子对象
    mol = Chem.MolFromSmiles(smi)
    # img = Draw.MolToImage(mol)
    # img.show()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    # 将分子对象转换为没有原子编号的SMILES
    smiles = Chem.MolToSmiles(mol, canonical=True)
    
    return smiles
    # print(smiles)
    
    # img.save("mol.jpg")

# function_1()

# 将一个 mol 中的多个分子分割开
def function_2(smi):

    # 创建包含两个分子的Mol对象
    mol = Chem.MolFromSmiles(smi)

    # 将Mol对象转换为Mol块字符串
    molblock = Chem.MolToMolBlock(mol)

    # 使用Chem.SDMolSupplier()函数读取Mol块字符串，并将其转换为多个Mol对象
    suppl = Chem.SDMolSupplier()
    suppl.SetData(molblock)
    mols = [mol for mol in suppl]

    return mols

    # # 输出每个分子的SMILES
    # for mol in mols:
    #     print(Chem.MolToSmiles(mol))



# function_2()

def smi_tokenizer(smi):
    """Tokenize a SMILES sequence or reaction"""
    pattern = "(\[[^\]]+]|Bi|Br?|Ge|Te|Mo|K|Ti|Zr|Y|Na|125I|Al|Ce|Cr|Cl?|Ni?|O|S|Pd?|Fe?|I|b|c|Mn|n|o|s|<unk>|>>|Li|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    if smi != ''.join(tokens):
        print('ERROR:', smi, ''.join(tokens))
    assert smi == ''.join(tokens)
    return tokens


def get_nonreactive_mask(cano_prod_am, raw_prod, raw_reacts, radius=0):
    """Retrieve the ground truth reaction center by RDChiral"""
    reactants = mols_from_smiles_list(replace_deuterated(raw_reacts).split('.'))
    # print(reactants)
    products = mols_from_smiles_list(replace_deuterated(raw_prod).split('.'))
    changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
    # print(changed_atom_tags)

    for _ in range(radius):
        mol = Chem.MolFromSmiles(cano_prod_am)
        changed_atom_tags_neighbor = []
        for atom in mol.GetAtoms():
            if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags:
                for n_atom in atom.GetNeighbors():
                    changed_atom_tags_neighbor.append(n_atom.GetSmarts().split(':')[1][:-1])
        changed_atom_tags = list(set(changed_atom_tags + changed_atom_tags_neighbor))
    # print(changed_atom_tags)
    tokens_list = []
    nonreactive_mask = []
    for i, token in enumerate(smi_tokenizer(cano_prod_am)):
        if token[0] == '[' and token[-1] == ']' and re.match('.*:([0-9]+)]', token):
            am = re.match('.*:([0-9]+)]', token).group(1)
            if am in changed_atom_tags:
                # print(token)
                # tokens_list.append(token)
                # nonreactive_mask.append(False)
                nonreactive_mask.append(True)
                continue
        # nonreactive_mask.append(True)
        nonreactive_mask.append(False)

    if sum(nonreactive_mask) == len(nonreactive_mask):  # if the reaction center is not detected
        # nonreactive_mask = [False] * len(nonreactive_mask)
        nonreactive_mask = [True] * len(nonreactive_mask)
    # print(nonreactive_mask)
    # print(tokens_list)
    return nonreactive_mask, changed_atom_tags


## 1. 从反应中提取模板
# smi = "[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])=[O:7].[NH2:8][c:9]1[cH:10][cH:11][cH:12][cH:13][c:14]1[F:15]"
# # mol = Chem.MolFromSmiles(smi)
# smi2 = "[CH3:1][C:2](=[O:3])[NH:4][c:5]1[cH:6][cH:7][cH:8][cH:9][c:10]1[F:11]"
# # mol2 = Chem.MolFromSmiles(smi2)
# # # smi2 = "CS(=O)(=O)OC[C@H]1CCC(=O)O1.Fc1ccc(Nc2ncnc3cc(OCCN4CCNCC4)c(OC4CCCC4)cc23)cc1Cl"
# # print(len(smi_tokenizer(smi2)))
# # # smi2 = "CS(=O)(=O)O[CH2:1][C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1.[F:8][c:9]1[cH:10][cH:11][c:12]([NH:13][c:14]2[n:15][cH:16][n:17][c:18]3[cH:19][c:20]([O:21][CH2:22][CH2:23][N:24]4[CH2:25][CH2:26][NH:27][CH2:28][CH2:29]4)[c:30]([O:31][CH:32]4[CH2:33][CH2:34][CH2:35][CH2:36]4)[cH:37][c:38]23)[cH:39][c:40]1[Cl:41]"
# # # print(len(smi_tokenizer(smi2)))
# print(len(smi_tokenizer(smi2)))
# print(len(get_nonreactive_mask(smi2, smi2, smi)))


# import pandas as pd
# import numpy as np
# path = 'python/rdkit_molecule/data/retro_uspto_50_template_5.csv'
# df = pd.read_csv(path)
# df_save = {
#     "products": [],
#     "reactants": [],
#     "products_mol": [],
#     "reactants_mol": [],
#     "reactive_smi" : [],
#     "set": [],
#     "reaction_type": []
# }

# for i in range(df.shape[0]):
#     p_smi = df["products"][i]
#     r_smi = df["reactants"][i]
#     p_smi_c = df["products_mol"][i]
#     r_smi_c = df["reactants_mol"][i]

#     df_save["products"].append(p_smi)
#     df_save["reactants"].append(r_smi)
#     df_save["products_mol"].append(p_smi_c)
#     df_save["reactants_mol"].append(r_smi_c)
#     df_save["set"].append(df["set"][i])
#     df_save["reaction_type"].append(df["reaction_type"][i])

#     if type(p_smi) != type(0.) and p_smi != "":
#         nonreactive, changed_atom_tags = get_nonreactive_mask(p_smi, p_smi, r_smi)
#         reactive = np.array(smi_tokenizer(p_smi_c))
#         df_save["reactive_smi"].append("".join(reactive[nonreactive].tolist()))
#     else:
#         df_save["reactive_smi"].append(p_smi_c)

# df_save_df = pd.DataFrame(df_save)

# df_save_df.to_csv("python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.csv")
# df_save_df.to_pickle("python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.pickle")
    

# 获取反应中心的原子
def get_reactive_atom():
    import pandas as pd
    import numpy as np
    path = 'python/rdkit_molecule/data/retro_uspto_50_template_4.csv'
    df = pd.read_csv(path)
    for i in range(df.shape[0]):
        p_smi = df["products"][i]
        r_smi = df["reactants"][i]
        p_smi_c = df["products_mol"][i]
        r_smi_c = df["reactants_mol"][i]

        nonreactive, changed_atom_tags = get_nonreactive_mask(p_smi, p_smi, r_smi)
        print(nonreactive)
        changed_atom_tag = [int(x) for x in changed_atom_tags]
        # print(changed_atom_tag)
        # print("p_smi_c", p_smi_c)
        # print("smi_tokenizer(p_smi_c)", len(smi_tokenizer(p_smi_c)))
        tokens = np.array(smi_tokenizer(p_smi_c))
        # print(tokens)
        print(tokens[nonreactive])
        # print(tokens[changed_atom_tag])
        break

# 将带原子映射的 SMILES 转换为标准 SMILES
def get_canonical_atom_nums(smiles_a_nums, smiles_b_canonical):
    # 转换为 RDKit 的 Molecule 对象
    mol_a = Chem.MolFromSmiles(smiles_a_nums)
    mol_b = Chem.MolFromSmiles(smiles_b_canonical)

    # 在 mol_a 中查找与 sub_mol_b 中的原子相对应的原子
    matches = mol_b.GetSubstructMatches(mol_a)

    # 输出匹配结果
    for match in matches:
        for i, atom_idx in enumerate(match):
            atom = mol_b.GetAtomWithIdx(atom_idx)
            atom.SetAtomMapNum(i + 1)
            # print(f'Atom {i}: {atom.GetSymbol()}, idx: {atom_idx}')

    return Chem.MolToSmiles(mol_b, canonical=False)



# mol = function_2(smi)
# mol2 = function_2(smi2)
# print(get_tetrahedral_atoms(mol, mol2))
# img = Draw.MolToImage(mol)
# img.show()

# img = Draw.MolToImage(mol2)
# img.show()

# from rdkit import Chem

# smiles = 'C1CC1C(Cl)(Br)F'
# mol = Chem.MolFromSmiles(smiles)

# # 获取SMILES字符串中原子的下标
# atom_indexes = [int(x) for x in re.findall(r'\d+', smiles)]
# print(atom_indexes)

# # 将SMILES字符串中的原子下标与分子对象中的索引对应
# atom_index_map = {atom_indexes[i]: i for i in range(len(atom_indexes))}

# # 打印原子下标和对应的索引号
# for i, atom in enumerate(mol.GetAtoms()):
#     print('Atom index in SMILES: %d, Atom index in Mol object: %d' % (atom_indexes[i], atom.GetIdx()))

def smi_tokens(smi):
    """
    Tokenize a SMILES molecule or reaction
    """
    import re

    pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return tokens

from rdkit import Chem

# 将 SMILES 中的原子符号下标与 Mol 中的原子序号对应
def smi_to_mol_index(smi):
    # 将 SMILES 转换为 Mol 对象
    mol = Chem.MolFromSmiles(smi)
    # 存储原子符号和序号的对应关系
    index_map = {}
    for atom in mol.GetAtoms():
        index_map[atom.GetSymbol()] = atom.GetIdx() + 1
    # 遍历 SMILES 中的原子符号，查找对应的原子序号，并替换原子符号下标
    for i in range(len(smi)):
        if smi[i].isdigit() and smi[i-1].isalpha():
            symbol = smi[i-1]
            index = index_map.get(symbol)
            if index:
                smi = smi[:i] + str(index) + smi[i+1:]
    return smi

import torch

def get_s2n():
    # 示例
    smiles = 'CC(C)C1CC[N:1]1C'
    mol = Chem.MolFromSmiles(smiles)
    nodes = len(mol.GetAtoms())
    # print()
    print(len(mol.GetAtoms()))
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    tokens = smi_tokens(smiles)

    smile2node = torch.zeros(nodes, len(tokens))
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    smi = Chem.MolToSmiles(mol)
    smi_ts = smi_tokens(smi)
    for i, tok in enumerate(smi_ts):
        if len(tok.strip().split(":")) == 2:
            smile2node[int(tok.strip().split(":")[1].split("]")[0]) - 1][i] = 1
    # img = Draw.MolToImage(mol)
    # img.show()
    print(smile2node)
    # smi = Chem.MolToSmiles(mol)
    # print(smi)

# smi = '[CH2:1]([C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1)[N:8]1[CH2:9][CH2:10][N:11]([CH2:12][CH2:13][O:14][c:15]2[cH:16][c:17]3[n:18][cH:19][n:20][c:21]([NH:22][c:23]4[cH:24][cH:25][c:26]([F:27])[c:28]([Cl:29])[cH:30]4)[c:31]3[cH:32][c:33]2[O:34][CH:35]2[CH2:36][CH2:37][CH2:38][CH2:39]2)[CH2:40][CH2:41]1'

from rdkit import Chem

def match_num_smi():
    # SMILES 字符串 a 和 b
    smiles_a = '[CH2:1]([C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1)[N:8]1[CH2:9][CH2:10][N:11]([CH2:12][CH2:13][O:14][c:15]2[cH:16][c:17]3[n:18][cH:19][n:20][c:21]([NH:22][c:23]4[cH:24][cH:25][c:26]([F:27])[c:28]([Cl:29])[cH:30]4)[c:31]3[cH:32][c:33]2[O:34][CH:35]2[CH2:36][CH2:37][CH2:38][CH2:39]2)[CH2:40][CH2:41]1'
    # smiles_b = 'O=C1CC[C@H](CN2CCN(CCOc3cc4ncnc(Nc5ccc(F)c(Cl)c5)c4cc3OC3CCCC3)CC2)O1'

    # 转换为 RDKit 的 Molecule 对象
    mol_a = Chem.MolFromSmiles(smiles_a)
    # mol_b = Chem.MolFromSmiles(smiles_b)

    # 获取 SMILES b 中的子图，即 c1ccccc1CCN 中的苯环和 N 原子
    sub_mol_b = Chem.MolFromSmarts('O=C1CC[C@H](CN2CCN(CCOc3cc4ncnc(Nc5ccc(F)c(Cl)c5)c4cc3OC3CCCC3)CC2)O1')

    # 在 mol_a 中查找与 sub_mol_b 中的原子相对应的原子
    matches = mol_a.GetSubstructMatches(sub_mol_b)

    # 输出匹配结果
    for match in matches:
        for i, atom_idx in enumerate(match):
            atom = mol_a.GetAtomWithIdx(atom_idx)
            print(f'Atom {i}: {atom.GetSymbol()}, idx: {atom_idx}')





def get_mol(smiles, kekulize=False):
    """SMILES string to Mol.

    Parameters
    ----------
    smiles: str,
        SMILES string for molecule
    kekulize: bool,
        Whether to kekulize the molecule
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None and kekulize:
        Chem.Kekulize(mol)
    return mol

def get_bond_info(mol):
    """Get information on bonds in the molecule.

    Parameters
    ----------
    mol: Chem.Mol
        Molecule
    """
    if mol is None:
        return {}

    bond_info = {}
    for bond in mol.GetBonds():
        a_start = bond.GetBeginAtom().GetAtomMapNum()
        a_end = bond.GetEndAtom().GetAtomMapNum()

        key_pair = sorted([a_start, a_end])
        bond_info[tuple(key_pair)] = [bond.GetBondTypeAsDouble(), bond.GetIdx()]

    return bond_info

def get_reaction_core(r, p, kekulize=False, use_h_labels=False):
    """Get the reaction core and edits for given reaction

    Parameters
    ----------
    r: str,
        SMILES string representing the reactants
    p: str,
        SMILES string representing the product
    kekulize: bool,
        Whether to kekulize molecules to fetch minimal set of edits
    use_h_labels: bool,
        Whether to use change in hydrogen counts in edits
    """
    reac_mol = get_mol(r)
    prod_mol = get_mol(p)

    if reac_mol is None or prod_mol is None:
        return set(), []

    prod_bonds = get_bond_info(prod_mol)
    p_amap_idx = {atom.GetAtomMapNum(): atom.GetIdx() for atom in prod_mol.GetAtoms()}

    max_amap = max([atom.GetAtomMapNum() for atom in reac_mol.GetAtoms()])
    for atom in reac_mol.GetAtoms():
        if atom.GetAtomMapNum() == 0:
            atom.SetAtomMapNum(max_amap + 1)
            max_amap += 1

    reac_bonds = get_bond_info(reac_mol)
    reac_amap = {atom.GetAtomMapNum(): atom.GetIdx() for atom in reac_mol.GetAtoms()}

    rxn_core = set()
    core_edits = []

    for bond in prod_bonds:
        if bond in reac_bonds and prod_bonds[bond][0] != reac_bonds[bond][0]:
            a_start, a_end = bond
            prod_bo, reac_bo = prod_bonds[bond][0], reac_bonds[bond][0]

            a_start, a_end = sorted([a_start, a_end])
            edit = f"{a_start}:{a_end}:{prod_bo}:{reac_bo}"
            core_edits.append(edit)
            rxn_core.update([a_start, a_end])

        if bond not in reac_bonds:
            a_start, a_end = bond
            reac_bo = 0.0
            prod_bo = prod_bonds[bond][0]

            start, end = sorted([a_start, a_end])
            edit = f"{a_start}:{a_end}:{prod_bo}:{reac_bo}"
            core_edits.append(edit)
            rxn_core.update([a_start, a_end])

    for bond in reac_bonds:
        if bond not in prod_bonds:
            amap1, amap2 = bond

            if (amap1 in p_amap_idx) and (amap2 in p_amap_idx):
                a_start, a_end = sorted([amap1, amap2])
                reac_bo = reac_bonds[bond][0]
                edit = f"{a_start}:{a_end}:{0.0}:{reac_bo}"
                core_edits.append(edit)
                rxn_core.update([a_start, a_end])

    if use_h_labels:
        if len(rxn_core) == 0:
            for atom in prod_mol.GetAtoms():
                amap_num = atom.GetAtomMapNum()

                numHs_prod = atom.GetTotalNumHs()
                numHs_reac = reac_mol.GetAtomWithIdx(reac_amap[amap_num]).GetTotalNumHs()

                if numHs_prod != numHs_reac:
                    edit = f"{amap_num}:{0}:{1.0}:{0.0}"
                    core_edits.append(edit)
                    rxn_core.add(amap_num)

    return rxn_core, core_edits

# reaction_smi = "[NH2:3][c:4]1[cH:5][cH:6][cH:7][c:8]2[cH:9][n:10][cH:11][cH:12][c:13]12.[O:1]=[C:2]([c:14]1[cH:15][c:16]([N+:17](=[O:18])[O-:19])[c:20]([S:21][c:22]2[c:23]([Cl:24])[cH:25][n:26][cH:27][c:28]2[Cl:29])[s:30]1)[OH:31]>>[O:1]=[C:2]([NH:3][c:4]1[cH:5][cH:6][cH:7][c:8]2[cH:9][n:10][cH:11][cH:12][c:13]12)[c:14]1[cH:15][c:16]([N+:17](=[O:18])[O-:19])[c:20]([S:21][c:22]2[c:23]([Cl:24])[cH:25][n:26][cH:27][c:28]2[Cl:29])[s:30]1"
# r, p = reaction_smi.split(">>")
# r = function_1(r)
# p = function_1(p)
# # rxn_core, core_edits = get_reaction_core(r=r, p=p)

# # print(rxn_core)
# # print(core_edits)

# reactants = mols_from_smiles_list(replace_deuterated(r).split('.'))
# products = mols_from_smiles_list(replace_deuterated(p).split('.'))
# print("products[0]", Chem.MolToSmiles(products[0]))
# print(len(products[0].GetAtoms()))
# print(len(smi_tokenizer(Chem.MolToSmiles(products[0]))))
# changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
# print(changed_atom_tags)
# print(changed_atoms)

def get_simility():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem import Draw
    from rdkit.Chem import rdFMCS
    mol1 = Chem.MolFromSmiles("c1ccc(C[N:8]2[CH2:9][CH2:10][N:11]([S:14](=[O:15])(=[O:16])[CH3:17])[CH2:12][CH2:13]2)cc1")
    mol2 = Chem.MolFromSmiles("c1ccc(C[N:6]2[CH2:5][CH:4]([CH2:3][N:2]([CH3:1])[CH3:16])[CH2:8][CH2:7]2)cc1")
    mols = [mol1, mol2]

    res = rdFMCS.FindMCS(mols).queryMol
    mols.append(res)

    mol1_ba_num = len(mol1.GetAtoms()) + len(mol1.GetBonds())
    mol2_ba_num = len(mol2.GetAtoms()) + len(mol2.GetBonds())
    res_ba_num = len(res.GetAtoms()) + len(res.GetBonds())

    print(mol1_ba_num)
    print(mol2_ba_num)
    print(res_ba_num)

    # img = Draw.MolsToGridImage(mols)
    # img

# get_simility()
# mol = Chem.MolFromSmiles(r)
# img = Draw.MolToImage(mol)
# img.show()

# mol = Chem.MolFromSmiles(p)
# img = Draw.MolToImage(mol)
# img.show()

# count = 0
# result_num = dict()
# with open("result_reactant.txt", "r") as f:
#     lines = f.readlines()

#     for i in range(len(lines)):
#         if i % 3 == 0:
#             true_smi = lines[i].strip()
#         elif i % 3 == 1:
#             pred_smi = lines[i].strip()
#         elif i % 3 == 2:
#             label = lines[i].strip().split("(")[1].split(",")[0]
#             if label == "False":
#                 for j in range(len(true_smi)):
#                     if j >= len(pred_smi):
#                         num = result_num.get(j-1, 0)
#                         result_num[j-1] = num + 1
#                         break
#                     if true_smi[j] != pred_smi[j]:
#                         num = result_num.get(j, 0)
#                         result_num[j] = num + 1
#                         break
                    
#                 count += 1

# key_num = [k for k in result_num.keys()]
# key_num.sort()

# print(count)
# for k in key_num:
#     print(k, result_num[k])
# print("result_num", result_num)


def match_retro_uspto():
    import pandas as pd
    import numpy as np
    path = 'python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.csv'
    df = pd.read_csv(path)

    df_train = pd.read_csv("python/rdkit_molecule/data/raw_train.csv")
    df_test = pd.read_csv("python/rdkit_molecule/data/raw_test.csv")
    df_val = pd.read_csv("python/rdkit_molecule/data/raw_val.csv")

    # flag = False
    for i in range(df.shape[0]):
        a = df["products"][i]
        
        # if type(a) == type(0.) and df["set"][i] == "train":
        if type(a) == type(0.):
            p_smi_c = df["products_mol"][i]
            r_smi_c = df["reactants_mol"][i]

            seek_flag = False

            for j in range(df_train.shape[0]):
                reaction = df_train["reactants>reagents>production"][j]
                r, _, p = reaction.split(">")
                p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
                r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
                if seek_flag:
                    break
                if p_smi_c == p_smi and r_smi_c == r_smi:
                    df["products"][i] = p
                    df["reactants"][i] = r
                    # df["products"][i] = get_canonical_atom_nums(p, p_smi)
                    # df["reactants"][i] = get_canonical_atom_nums(r, r_smi)
                    seek_flag = True
                    print(i)
                    break

            for j in range(df_test.shape[0]):
                if seek_flag:
                    break
                reaction = df_test["reactants>reagents>production"][j]
                r, _, p = reaction.split(">")
                p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
                r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
                if p_smi_c == p_smi and r_smi_c == r_smi:
                    df["products"][i] = p
                    df["reactants"][i] = r
                    # df["products"][i] = get_canonical_atom_nums(p, p_smi)
                    # df["reactants"][i] = get_canonical_atom_nums(r, r_smi)
                    seek_flags = True
                    print(i)
                    break

            for j in range(df_val.shape[0]):
                if seek_flag:
                    break
                reaction = df_val["reactants>reagents>production"][j]
                r, _, p = reaction.split(">")
                p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
                r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
                if p_smi_c == p_smi and r_smi_c == r_smi:
                    df["products"][i] = p
                    df["reactants"][i] = r
                    # df["products"][i] = get_canonical_atom_nums(p, p_smi)
                    # df["reactants"][i] = get_canonical_atom_nums(r, r_smi)
                    # flag = True
                    print(i)
                    break
            # if flag:
            #     break
        # if i % 1000 == 0:
        #     print(i)
    
    df.to_csv("python/rdkit_molecule/data/retro_uspto_50_template_reaction_retro_2.csv")
    df.to_pickle("python/rdkit_molecule/data/retro_uspto_50_template_reaction_retro_2.pickle")

# match_retro_uspto()

# import pandas as pd
# import numpy as np
# path = 'python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.csv'
# df = pd.read_csv(path)

# train_count = 0
# valid_count = 0
# test_count = 0
# for i in range(df.shape[0]):
#     if df["products"][i] == "nan" or type(df["products"][i]) == type(0.) or df["products"][i] == "":
#         if df["set"][i] == "train":
#             train_count += 1
#         elif df["set"][i] == "valid":
#             valid_count += 1
#         elif df["set"][i] == "test":
#             test_count += 1
# print(train_count)
# print(valid_count)
# print(test_count)

# df_train = pd.read_csv("python/rdkit_molecule/data/raw_train.csv")

# # p_smi_c = df["products_mol"][263]
# # r_smi_c = df["reactants_mol"][263]
# count = 0
# for i in range(df.shape[0]):
#     if type(df["products"][i]) == type(0.) and df["set"][i] == "train":
#         # if i == 13548:
#         #     print(count)
#         count += 1
# # print(count)
#         reaction = df_train["reactants>reagents>production"][j]
#         r, _, p = reaction.split(">")
#         p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
#         r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
    #     if p_smi_c == p_smi and r_smi_c == r_smi:
    #         df["products"][263] = get_canonical_atom_nums(p, p_smi)
    #         df["reactants"][263] = get_canonical_atom_nums(r, r_smi)
    #         # if p_smi != p_smi_c or r_smi != r_smi_c:
    #         print(p)
    #         print(r_smi_c)
    #         print(df["products"][263])
    #         print(df["reactants"][263])

def detele_50k():
    import pandas as pd
    path1 = "data/uspto_full_deal.pickle"
    path2 = "F:/code/learning_language/python/rdkit_molecule/data/uspto_50.pickle"

    df1 = pd.read_pickle(path1)
    df2 = pd.read_pickle(path2)

    set_50k = set()

    dict_full = {
        "reactants_mol": [],
        "products_mol": [],
        "reaction_type": [],
        "set": []
    }

    for i in range(df2.shape[0]):
        if df2["set"][i] == "train":
            continue
        set_50k.add(Chem.MolToSmiles(df2["reactants_mol"][i]) + "!" + Chem.MolToSmiles(df2["products_mol"][i]))
    count = 0
    for i in range(df1.shape[0]):
        if df1["set"][i] != "train":
            continue
        smi = df1["reactants_mol"][i] + "!" +  df1["products_mol"][i]
        if smi not in set_50k:
            dict_full["reactants_mol"].append(df1["reactants_mol"][i])
            dict_full["products_mol"].append(df1["products_mol"][i])
            dict_full["reaction_type"].append(df1["reaction_type"][i])
            dict_full["set"].append(df1["set"][i])
        else:
            count += 1
    print(count)
    for i in range(df2.shape[0]):
        if df2["set"][i] != "train":
            dict_full["reactants_mol"].append(df2["reactants_mol"][i])
            dict_full["products_mol"].append(df2["products_mol"][i])
            dict_full["reaction_type"].append(df2["reaction_type"][i])
            dict_full["set"].append(df2["set"][i])
    
    df_full = pd.DataFrame(dict_full)
    df_full.to_pickle("data/uspto_full_no50_1m.pickle")

# detele_50k()

def deal_pretrain_cross_graph_text():
    import pandas as pd
    df_full = pd.read_pickle("data/uspto_full_no50_1m.pickle")
    df_50k = pd.read_pickle("data/uspto_50.pickle")

    dict_pretrain = {
        "smiles": [],
        "atom_features": [],
        "edges": [],
        "adjacencys": [],
        "node2smiles": [],
        "set": [],
    }
    smiles_set = set()

    for i in range(df_50k.shape[0]):
        reactants_smi = Chem.MolToSmiles(df_50k["reactants_mol"][i])
        products_smi = Chem.MolToSmiles(df_50k["products_mol"][i])
        if reactants_smi not in smiles_set:
            atom_feature, edge, adj, smile2node = get_graph_features_from_smi(reactants_smi)

            dict_pretrain["smiles"].append(reactants_smi)
            dict_pretrain["atom_features"].append(atom_feature)
            dict_pretrain["edges"].append(edge)
            dict_pretrain["adjacencys"].append(adj)
            dict_pretrain["node2smiles"].append(smile2node)
            dict_pretrain["set"].append(df_50k["set"][i])
            smiles_set.add(reactants_smi)
        
        if products_smi not in smiles_set:
            atom_feature, edge, adj, smile2node = get_graph_features_from_smi(products_smi)

            dict_pretrain["smiles"].append(products_smi)
            dict_pretrain["atom_features"].append(atom_feature)
            dict_pretrain["edges"].append(edge)
            dict_pretrain["adjacencys"].append(adj)
            dict_pretrain["node2smiles"].append(smile2node)
            dict_pretrain["set"].append(df_50k["set"][i])
            smiles_set.add(products_smi)
    
    df = pd.DataFrame(dict_pretrain)
    print(df.shape[0])
    df.to_pickle("data/uspto_pretrain_50k.pickle")


# deal_pretrain_cross_graph_text()
# import pandas as pd
# df = pd.read_pickle("data/uspto_full_no50_1m.pickle")
# # for i in range(df.shape[0]):
# #     if df["set"][i] == "train":
# #         df["reactants_mol"][i] = Chem.MolFromSmiles(df["reactants_mol"][i])
# #         df["products_mol"][i] = Chem.MolFromSmiles(df["products_mol"][i])
# # df.to_pickle("data/uspto_full_no50_1m_mol.pickle")
# print(df)

from calendar import c
from math import isnan, nan
from sre_constants import CH_LOCALE
from rdkit import Chem
from rdkit.Chem import Draw

import torch
import re

# from rdchiral.template_extractor import extract_from_reaction, get_changed_atoms, mols_from_smiles_list, \
#     replace_deuterated, get_tetrahedral_atoms
# from sklearn.metrics import recall_score

# from molecule import get_graph_features_from_smi


# 去掉原子编号
def function_1(smi):
    # 将带有原子编号的SMILES转换为分子对象
    mol = Chem.MolFromSmiles(smi)
    # img = Draw.MolToImage(mol)
    # img.show()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    # 将分子对象转换为没有原子编号的SMILES
    # smiles = Chem.MolToSmiles(mol, canonical=True)
    
    # return smiles
    return mol
#     # print(smiles)
    
#     # img.save("mol.jpg")

# # function_1()

# # 将一个 mol 中的多个分子分割开
# def function_2(smi):

#     # 创建包含两个分子的Mol对象
#     mol = Chem.MolFromSmiles(smi)

#     # 将Mol对象转换为Mol块字符串
#     molblock = Chem.MolToMolBlock(mol)

#     # 使用Chem.SDMolSupplier()函数读取Mol块字符串，并将其转换为多个Mol对象
#     suppl = Chem.SDMolSupplier()
#     suppl.SetData(molblock)
#     mols = [mol for mol in suppl]

#     return mols

#     # # 输出每个分子的SMILES
#     # for mol in mols:
#     #     print(Chem.MolToSmiles(mol))



# # function_2()

# def smi_tokenizer(smi):
#     """Tokenize a SMILES sequence or reaction"""
#     pattern = "(\[[^\]]+]|Bi|Br?|Ge|Te|Mo|K|Ti|Zr|Y|Na|125I|Al|Ce|Cr|Cl?|Ni?|O|S|Pd?|Fe?|I|b|c|Mn|n|o|s|<unk>|>>|Li|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
#     regex = re.compile(pattern)
#     tokens = [token for token in regex.findall(smi)]
#     if smi != ''.join(tokens):
#         print('ERROR:', smi, ''.join(tokens))
#     assert smi == ''.join(tokens)
#     return tokens


# def get_nonreactive_mask(cano_prod_am, raw_prod, raw_reacts, radius=0):
#     """Retrieve the ground truth reaction center by RDChiral"""
#     reactants = mols_from_smiles_list(replace_deuterated(raw_reacts).split('.'))
#     # print(reactants)
#     products = mols_from_smiles_list(replace_deuterated(raw_prod).split('.'))
#     changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
#     # print(changed_atom_tags)

#     for _ in range(radius):
#         mol = Chem.MolFromSmiles(cano_prod_am)
#         changed_atom_tags_neighbor = []
#         for atom in mol.GetAtoms():
#             if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags:
#                 for n_atom in atom.GetNeighbors():
#                     changed_atom_tags_neighbor.append(n_atom.GetSmarts().split(':')[1][:-1])
#         changed_atom_tags = list(set(changed_atom_tags + changed_atom_tags_neighbor))
#     # print(changed_atom_tags)
#     tokens_list = []
#     nonreactive_mask = []
#     for i, token in enumerate(smi_tokenizer(cano_prod_am)):
#         if token[0] == '[' and token[-1] == ']' and re.match('.*:([0-9]+)]', token):
#             am = re.match('.*:([0-9]+)]', token).group(1)
#             if am in changed_atom_tags:
#                 # print(token)
#                 # tokens_list.append(token)
#                 # nonreactive_mask.append(False)
#                 nonreactive_mask.append(True)
#                 continue
#         # nonreactive_mask.append(True)
#         nonreactive_mask.append(False)

#     if sum(nonreactive_mask) == len(nonreactive_mask):  # if the reaction center is not detected
#         # nonreactive_mask = [False] * len(nonreactive_mask)
#         nonreactive_mask = [True] * len(nonreactive_mask)
#     # print(nonreactive_mask)
#     # print(tokens_list)
#     return nonreactive_mask, changed_atom_tags


# ## 1. 从反应中提取模板
# # smi = "[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])=[O:7].[NH2:8][c:9]1[cH:10][cH:11][cH:12][cH:13][c:14]1[F:15]"
# # # mol = Chem.MolFromSmiles(smi)
# # smi2 = "[CH3:1][C:2](=[O:3])[NH:4][c:5]1[cH:6][cH:7][cH:8][cH:9][c:10]1[F:11]"
# # # mol2 = Chem.MolFromSmiles(smi2)
# # # # smi2 = "CS(=O)(=O)OC[C@H]1CCC(=O)O1.Fc1ccc(Nc2ncnc3cc(OCCN4CCNCC4)c(OC4CCCC4)cc23)cc1Cl"
# # # print(len(smi_tokenizer(smi2)))
# # # # smi2 = "CS(=O)(=O)O[CH2:1][C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1.[F:8][c:9]1[cH:10][cH:11][c:12]([NH:13][c:14]2[n:15][cH:16][n:17][c:18]3[cH:19][c:20]([O:21][CH2:22][CH2:23][N:24]4[CH2:25][CH2:26][NH:27][CH2:28][CH2:29]4)[c:30]([O:31][CH:32]4[CH2:33][CH2:34][CH2:35][CH2:36]4)[cH:37][c:38]23)[cH:39][c:40]1[Cl:41]"
# # # # print(len(smi_tokenizer(smi2)))
# # print(len(smi_tokenizer(smi2)))
# # print(len(get_nonreactive_mask(smi2, smi2, smi)))


# # import pandas as pd
# # import numpy as np
# # path = 'python/rdkit_molecule/data/retro_uspto_50_template_5.csv'
# # df = pd.read_csv(path)
# # df_save = {
# #     "products": [],
# #     "reactants": [],
# #     "products_mol": [],
# #     "reactants_mol": [],
# #     "reactive_smi" : [],
# #     "set": [],
# #     "reaction_type": []
# # }

# # for i in range(df.shape[0]):
# #     p_smi = df["products"][i]
# #     r_smi = df["reactants"][i]
# #     p_smi_c = df["products_mol"][i]
# #     r_smi_c = df["reactants_mol"][i]

# #     df_save["products"].append(p_smi)
# #     df_save["reactants"].append(r_smi)
# #     df_save["products_mol"].append(p_smi_c)
# #     df_save["reactants_mol"].append(r_smi_c)
# #     df_save["set"].append(df["set"][i])
# #     df_save["reaction_type"].append(df["reaction_type"][i])

# #     if type(p_smi) != type(0.) and p_smi != "":
# #         nonreactive, changed_atom_tags = get_nonreactive_mask(p_smi, p_smi, r_smi)
# #         reactive = np.array(smi_tokenizer(p_smi_c))
# #         df_save["reactive_smi"].append("".join(reactive[nonreactive].tolist()))
# #     else:
# #         df_save["reactive_smi"].append(p_smi_c)

# # df_save_df = pd.DataFrame(df_save)

# # df_save_df.to_csv("python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.csv")
# # df_save_df.to_pickle("python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.pickle")
    

# # 获取反应中心的原子
# def get_reactive_atom():
#     import pandas as pd
#     import numpy as np
#     path = 'python/rdkit_molecule/data/retro_uspto_50_template_4.csv'
#     df = pd.read_csv(path)
#     for i in range(df.shape[0]):
#         p_smi = df["products"][i]
#         r_smi = df["reactants"][i]
#         p_smi_c = df["products_mol"][i]
#         r_smi_c = df["reactants_mol"][i]

#         nonreactive, changed_atom_tags = get_nonreactive_mask(p_smi, p_smi, r_smi)
#         print(nonreactive)
#         changed_atom_tag = [int(x) for x in changed_atom_tags]
#         # print(changed_atom_tag)
#         # print("p_smi_c", p_smi_c)
#         # print("smi_tokenizer(p_smi_c)", len(smi_tokenizer(p_smi_c)))
#         tokens = np.array(smi_tokenizer(p_smi_c))
#         # print(tokens)
#         print(tokens[nonreactive])
#         # print(tokens[changed_atom_tag])
#         break

# # 将带原子映射的 SMILES 转换为标准 SMILES
# def get_canonical_atom_nums(smiles_a_nums, smiles_b_canonical):
#     # 转换为 RDKit 的 Molecule 对象
#     mol_a = Chem.MolFromSmiles(smiles_a_nums)
#     mol_b = Chem.MolFromSmiles(smiles_b_canonical)

#     # 在 mol_a 中查找与 sub_mol_b 中的原子相对应的原子
#     matches = mol_b.GetSubstructMatches(mol_a)

#     # 输出匹配结果
#     for match in matches:
#         for i, atom_idx in enumerate(match):
#             atom = mol_b.GetAtomWithIdx(atom_idx)
#             atom.SetAtomMapNum(i + 1)
#             # print(f'Atom {i}: {atom.GetSymbol()}, idx: {atom_idx}')

#     return Chem.MolToSmiles(mol_b, canonical=False)



# # mol = function_2(smi)
# # mol2 = function_2(smi2)
# # print(get_tetrahedral_atoms(mol, mol2))
# # img = Draw.MolToImage(mol)
# # img.show()

# # img = Draw.MolToImage(mol2)
# # img.show()

# # from rdkit import Chem

# # smiles = 'C1CC1C(Cl)(Br)F'
# # mol = Chem.MolFromSmiles(smiles)

# # # 获取SMILES字符串中原子的下标
# # atom_indexes = [int(x) for x in re.findall(r'\d+', smiles)]
# # print(atom_indexes)

# # # 将SMILES字符串中的原子下标与分子对象中的索引对应
# # atom_index_map = {atom_indexes[i]: i for i in range(len(atom_indexes))}

# # # 打印原子下标和对应的索引号
# # for i, atom in enumerate(mol.GetAtoms()):
# #     print('Atom index in SMILES: %d, Atom index in Mol object: %d' % (atom_indexes[i], atom.GetIdx()))

# def smi_tokens(smi):
#     """
#     Tokenize a SMILES molecule or reaction
#     """
#     import re

#     pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
#     regex = re.compile(pattern)
#     tokens = [token for token in regex.findall(smi)]
#     assert smi == ''.join(tokens)
#     return tokens

# from rdkit import Chem

# # 将 SMILES 中的原子符号下标与 Mol 中的原子序号对应
# def smi_to_mol_index(smi):
#     # 将 SMILES 转换为 Mol 对象
#     mol = Chem.MolFromSmiles(smi)
#     # 存储原子符号和序号的对应关系
#     index_map = {}
#     for atom in mol.GetAtoms():
#         index_map[atom.GetSymbol()] = atom.GetIdx() + 1
#     # 遍历 SMILES 中的原子符号，查找对应的原子序号，并替换原子符号下标
#     for i in range(len(smi)):
#         if smi[i].isdigit() and smi[i-1].isalpha():
#             symbol = smi[i-1]
#             index = index_map.get(symbol)
#             if index:
#                 smi = smi[:i] + str(index) + smi[i+1:]
#     return smi

# import torch

# def get_s2n():
#     # 示例
#     smiles = 'CC(C)C1CC[N:1]1C'
#     mol = Chem.MolFromSmiles(smiles)
#     nodes = len(mol.GetAtoms())
#     # print()
#     print(len(mol.GetAtoms()))
#     for atom in mol.GetAtoms():
#         atom.SetAtomMapNum(atom.GetIdx() + 1)

#     tokens = smi_tokens(smiles)

#     smile2node = torch.zeros(nodes, len(tokens))
#     for atom in mol.GetAtoms():
#         atom.SetAtomMapNum(atom.GetIdx() + 1)
#     smi = Chem.MolToSmiles(mol)
#     smi_ts = smi_tokens(smi)
#     for i, tok in enumerate(smi_ts):
#         if len(tok.strip().split(":")) == 2:
#             smile2node[int(tok.strip().split(":")[1].split("]")[0]) - 1][i] = 1
#     # img = Draw.MolToImage(mol)
#     # img.show()
#     print(smile2node)
#     # smi = Chem.MolToSmiles(mol)
#     # print(smi)

# # smi = '[CH2:1]([C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1)[N:8]1[CH2:9][CH2:10][N:11]([CH2:12][CH2:13][O:14][c:15]2[cH:16][c:17]3[n:18][cH:19][n:20][c:21]([NH:22][c:23]4[cH:24][cH:25][c:26]([F:27])[c:28]([Cl:29])[cH:30]4)[c:31]3[cH:32][c:33]2[O:34][CH:35]2[CH2:36][CH2:37][CH2:38][CH2:39]2)[CH2:40][CH2:41]1'

# from rdkit import Chem

# def match_num_smi():
#     # SMILES 字符串 a 和 b
#     smiles_a = '[CH2:1]([C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1)[N:8]1[CH2:9][CH2:10][N:11]([CH2:12][CH2:13][O:14][c:15]2[cH:16][c:17]3[n:18][cH:19][n:20][c:21]([NH:22][c:23]4[cH:24][cH:25][c:26]([F:27])[c:28]([Cl:29])[cH:30]4)[c:31]3[cH:32][c:33]2[O:34][CH:35]2[CH2:36][CH2:37][CH2:38][CH2:39]2)[CH2:40][CH2:41]1'
#     # smiles_b = 'O=C1CC[C@H](CN2CCN(CCOc3cc4ncnc(Nc5ccc(F)c(Cl)c5)c4cc3OC3CCCC3)CC2)O1'

#     # 转换为 RDKit 的 Molecule 对象
#     mol_a = Chem.MolFromSmiles(smiles_a)
#     # mol_b = Chem.MolFromSmiles(smiles_b)

#     # 获取 SMILES b 中的子图，即 c1ccccc1CCN 中的苯环和 N 原子
#     sub_mol_b = Chem.MolFromSmarts('O=C1CC[C@H](CN2CCN(CCOc3cc4ncnc(Nc5ccc(F)c(Cl)c5)c4cc3OC3CCCC3)CC2)O1')

#     # 在 mol_a 中查找与 sub_mol_b 中的原子相对应的原子
#     matches = mol_a.GetSubstructMatches(sub_mol_b)

#     # 输出匹配结果
#     for match in matches:
#         for i, atom_idx in enumerate(match):
#             atom = mol_a.GetAtomWithIdx(atom_idx)
#             print(f'Atom {i}: {atom.GetSymbol()}, idx: {atom_idx}')





# def get_mol(smiles, kekulize=False):
#     """SMILES string to Mol.

#     Parameters
#     ----------
#     smiles: str,
#         SMILES string for molecule
#     kekulize: bool,
#         Whether to kekulize the molecule
#     """
#     mol = Chem.MolFromSmiles(smiles)
#     if mol is not None and kekulize:
#         Chem.Kekulize(mol)
#     return mol

# def get_bond_info(mol):
#     """Get information on bonds in the molecule.

#     Parameters
#     ----------
#     mol: Chem.Mol
#         Molecule
#     """
#     if mol is None:
#         return {}

#     bond_info = {}
#     for bond in mol.GetBonds():
#         a_start = bond.GetBeginAtom().GetAtomMapNum()
#         a_end = bond.GetEndAtom().GetAtomMapNum()

#         key_pair = sorted([a_start, a_end])
#         bond_info[tuple(key_pair)] = [bond.GetBondTypeAsDouble(), bond.GetIdx()]

#     return bond_info

# def get_reaction_core(r, p, kekulize=False, use_h_labels=False):
#     """Get the reaction core and edits for given reaction

#     Parameters
#     ----------
#     r: str,
#         SMILES string representing the reactants
#     p: str,
#         SMILES string representing the product
#     kekulize: bool,
#         Whether to kekulize molecules to fetch minimal set of edits
#     use_h_labels: bool,
#         Whether to use change in hydrogen counts in edits
#     """
#     reac_mol = get_mol(r)
#     prod_mol = get_mol(p)

#     if reac_mol is None or prod_mol is None:
#         return set(), []

#     prod_bonds = get_bond_info(prod_mol)
#     p_amap_idx = {atom.GetAtomMapNum(): atom.GetIdx() for atom in prod_mol.GetAtoms()}

#     max_amap = max([atom.GetAtomMapNum() for atom in reac_mol.GetAtoms()])
#     for atom in reac_mol.GetAtoms():
#         if atom.GetAtomMapNum() == 0:
#             atom.SetAtomMapNum(max_amap + 1)
#             max_amap += 1

#     reac_bonds = get_bond_info(reac_mol)
#     reac_amap = {atom.GetAtomMapNum(): atom.GetIdx() for atom in reac_mol.GetAtoms()}

#     rxn_core = set()
#     core_edits = []

#     for bond in prod_bonds:
#         if bond in reac_bonds and prod_bonds[bond][0] != reac_bonds[bond][0]:
#             a_start, a_end = bond
#             prod_bo, reac_bo = prod_bonds[bond][0], reac_bonds[bond][0]

#             a_start, a_end = sorted([a_start, a_end])
#             edit = f"{a_start}:{a_end}:{prod_bo}:{reac_bo}"
#             core_edits.append(edit)
#             rxn_core.update([a_start, a_end])

#         if bond not in reac_bonds:
#             a_start, a_end = bond
#             reac_bo = 0.0
#             prod_bo = prod_bonds[bond][0]

#             start, end = sorted([a_start, a_end])
#             edit = f"{a_start}:{a_end}:{prod_bo}:{reac_bo}"
#             core_edits.append(edit)
#             rxn_core.update([a_start, a_end])

#     for bond in reac_bonds:
#         if bond not in prod_bonds:
#             amap1, amap2 = bond

#             if (amap1 in p_amap_idx) and (amap2 in p_amap_idx):
#                 a_start, a_end = sorted([amap1, amap2])
#                 reac_bo = reac_bonds[bond][0]
#                 edit = f"{a_start}:{a_end}:{0.0}:{reac_bo}"
#                 core_edits.append(edit)
#                 rxn_core.update([a_start, a_end])

#     if use_h_labels:
#         if len(rxn_core) == 0:
#             for atom in prod_mol.GetAtoms():
#                 amap_num = atom.GetAtomMapNum()

#                 numHs_prod = atom.GetTotalNumHs()
#                 numHs_reac = reac_mol.GetAtomWithIdx(reac_amap[amap_num]).GetTotalNumHs()

#                 if numHs_prod != numHs_reac:
#                     edit = f"{amap_num}:{0}:{1.0}:{0.0}"
#                     core_edits.append(edit)
#                     rxn_core.add(amap_num)

#     return rxn_core, core_edits

# # reaction_smi = "[NH2:3][c:4]1[cH:5][cH:6][cH:7][c:8]2[cH:9][n:10][cH:11][cH:12][c:13]12.[O:1]=[C:2]([c:14]1[cH:15][c:16]([N+:17](=[O:18])[O-:19])[c:20]([S:21][c:22]2[c:23]([Cl:24])[cH:25][n:26][cH:27][c:28]2[Cl:29])[s:30]1)[OH:31]>>[O:1]=[C:2]([NH:3][c:4]1[cH:5][cH:6][cH:7][c:8]2[cH:9][n:10][cH:11][cH:12][c:13]12)[c:14]1[cH:15][c:16]([N+:17](=[O:18])[O-:19])[c:20]([S:21][c:22]2[c:23]([Cl:24])[cH:25][n:26][cH:27][c:28]2[Cl:29])[s:30]1"
# # r, p = reaction_smi.split(">>")
# # r = function_1(r)
# # p = function_1(p)
# # # rxn_core, core_edits = get_reaction_core(r=r, p=p)

# # # print(rxn_core)
# # # print(core_edits)

# # reactants = mols_from_smiles_list(replace_deuterated(r).split('.'))
# # products = mols_from_smiles_list(replace_deuterated(p).split('.'))
# # print("products[0]", Chem.MolToSmiles(products[0]))
# # print(len(products[0].GetAtoms()))
# # print(len(smi_tokenizer(Chem.MolToSmiles(products[0]))))
# # changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
# # print(changed_atom_tags)
# # print(changed_atoms)

# def get_simility():
#     from rdkit import Chem
#     from rdkit.Chem import AllChem
#     from rdkit.Chem import Draw
#     from rdkit.Chem import rdFMCS
#     mol1 = Chem.MolFromSmiles("c1ccc(C[N:8]2[CH2:9][CH2:10][N:11]([S:14](=[O:15])(=[O:16])[CH3:17])[CH2:12][CH2:13]2)cc1")
#     mol2 = Chem.MolFromSmiles("c1ccc(C[N:6]2[CH2:5][CH:4]([CH2:3][N:2]([CH3:1])[CH3:16])[CH2:8][CH2:7]2)cc1")
#     mols = [mol1, mol2]

#     res = rdFMCS.FindMCS(mols).queryMol
#     mols.append(res)

#     mol1_ba_num = len(mol1.GetAtoms()) + len(mol1.GetBonds())
#     mol2_ba_num = len(mol2.GetAtoms()) + len(mol2.GetBonds())
#     res_ba_num = len(res.GetAtoms()) + len(res.GetBonds())

#     print(mol1_ba_num)
#     print(mol2_ba_num)
#     print(res_ba_num)

#     # img = Draw.MolsToGridImage(mols)
#     # img

# # get_simility()
# # mol = Chem.MolFromSmiles(r)
# # img = Draw.MolToImage(mol)
# # img.show()

# # mol = Chem.MolFromSmiles(p)
# # img = Draw.MolToImage(mol)
# # img.show()

# # count = 0
# # result_num = dict()
# # with open("result_reactant.txt", "r") as f:
# #     lines = f.readlines()

# #     for i in range(len(lines)):
# #         if i % 3 == 0:
# #             true_smi = lines[i].strip()
# #         elif i % 3 == 1:
# #             pred_smi = lines[i].strip()
# #         elif i % 3 == 2:
# #             label = lines[i].strip().split("(")[1].split(",")[0]
# #             if label == "False":
# #                 for j in range(len(true_smi)):
# #                     if j >= len(pred_smi):
# #                         num = result_num.get(j-1, 0)
# #                         result_num[j-1] = num + 1
# #                         break
# #                     if true_smi[j] != pred_smi[j]:
# #                         num = result_num.get(j, 0)
# #                         result_num[j] = num + 1
# #                         break
                    
# #                 count += 1

# # key_num = [k for k in result_num.keys()]
# # key_num.sort()

# # print(count)
# # for k in key_num:
# #     print(k, result_num[k])
# # print("result_num", result_num)


# def match_retro_uspto():
#     import pandas as pd
#     import numpy as np
#     path = 'python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.csv'
#     df = pd.read_csv(path)

#     df_train = pd.read_csv("python/rdkit_molecule/data/raw_train.csv")
#     df_test = pd.read_csv("python/rdkit_molecule/data/raw_test.csv")
#     df_val = pd.read_csv("python/rdkit_molecule/data/raw_val.csv")

#     # flag = False
#     for i in range(df.shape[0]):
#         a = df["products"][i]
        
#         # if type(a) == type(0.) and df["set"][i] == "train":
#         if type(a) == type(0.):
#             p_smi_c = df["products_mol"][i]
#             r_smi_c = df["reactants_mol"][i]

#             seek_flag = False

#             for j in range(df_train.shape[0]):
#                 reaction = df_train["reactants>reagents>production"][j]
#                 r, _, p = reaction.split(">")
#                 p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
#                 r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
#                 if seek_flag:
#                     break
#                 if p_smi_c == p_smi and r_smi_c == r_smi:
#                     df["products"][i] = p
#                     df["reactants"][i] = r
#                     # df["products"][i] = get_canonical_atom_nums(p, p_smi)
#                     # df["reactants"][i] = get_canonical_atom_nums(r, r_smi)
#                     seek_flag = True
#                     print(i)
#                     break

#             for j in range(df_test.shape[0]):
#                 if seek_flag:
#                     break
#                 reaction = df_test["reactants>reagents>production"][j]
#                 r, _, p = reaction.split(">")
#                 p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
#                 r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
#                 if p_smi_c == p_smi and r_smi_c == r_smi:
#                     df["products"][i] = p
#                     df["reactants"][i] = r
#                     # df["products"][i] = get_canonical_atom_nums(p, p_smi)
#                     # df["reactants"][i] = get_canonical_atom_nums(r, r_smi)
#                     seek_flags = True
#                     print(i)
#                     break

#             for j in range(df_val.shape[0]):
#                 if seek_flag:
#                     break
#                 reaction = df_val["reactants>reagents>production"][j]
#                 r, _, p = reaction.split(">")
#                 p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
#                 r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
#                 if p_smi_c == p_smi and r_smi_c == r_smi:
#                     df["products"][i] = p
#                     df["reactants"][i] = r
#                     # df["products"][i] = get_canonical_atom_nums(p, p_smi)
#                     # df["reactants"][i] = get_canonical_atom_nums(r, r_smi)
#                     # flag = True
#                     print(i)
#                     break
#             # if flag:
#             #     break
#         # if i % 1000 == 0:
#         #     print(i)
    
#     df.to_csv("python/rdkit_molecule/data/retro_uspto_50_template_reaction_retro_2.csv")
#     df.to_pickle("python/rdkit_molecule/data/retro_uspto_50_template_reaction_retro_2.pickle")

# # match_retro_uspto()

# # import pandas as pd
# # import numpy as np
# # path = 'python/rdkit_molecule/data/retro_uspto_50_template_reactive_5.csv'
# # df = pd.read_csv(path)

# # train_count = 0
# # valid_count = 0
# # test_count = 0
# # for i in range(df.shape[0]):
# #     if df["products"][i] == "nan" or type(df["products"][i]) == type(0.) or df["products"][i] == "":
# #         if df["set"][i] == "train":
# #             train_count += 1
# #         elif df["set"][i] == "valid":
# #             valid_count += 1
# #         elif df["set"][i] == "test":
# #             test_count += 1
# # print(train_count)
# # print(valid_count)
# # print(test_count)

# # df_train = pd.read_csv("python/rdkit_molecule/data/raw_train.csv")

# # # p_smi_c = df["products_mol"][263]
# # # r_smi_c = df["reactants_mol"][263]
# # count = 0
# # for i in range(df.shape[0]):
# #     if type(df["products"][i]) == type(0.) and df["set"][i] == "train":
# #         # if i == 13548:
# #         #     print(count)
# #         count += 1
# # # print(count)
# #         reaction = df_train["reactants>reagents>production"][j]
# #         r, _, p = reaction.split(">")
# #         p_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(p)))
# #         r_smi = Chem.MolToSmiles(Chem.MolFromSmiles(function_1(r)))
#     #     if p_smi_c == p_smi and r_smi_c == r_smi:
#     #         df["products"][263] = get_canonical_atom_nums(p, p_smi)
#     #         df["reactants"][263] = get_canonical_atom_nums(r, r_smi)
#     #         # if p_smi != p_smi_c or r_smi != r_smi_c:
#     #         print(p)
#     #         print(r_smi_c)
#     #         print(df["products"][263])
#     #         print(df["reactants"][263])

# def detele_50k():
#     import pandas as pd
#     path1 = "data/uspto_full_deal.pickle"
#     path2 = "F:/code/learning_language/python/rdkit_molecule/data/uspto_50.pickle"

#     df1 = pd.read_pickle(path1)
#     df2 = pd.read_pickle(path2)

#     set_50k = set()

#     dict_full = {
#         "reactants_mol": [],
#         "products_mol": [],
#         "reaction_type": [],
#         "set": []
#     }

#     for i in range(df2.shape[0]):
#         if df2["set"][i] == "train":
#             continue
#         set_50k.add(Chem.MolToSmiles(df2["reactants_mol"][i]) + "!" + Chem.MolToSmiles(df2["products_mol"][i]))
#     count = 0
#     for i in range(df1.shape[0]):
#         if df1["set"][i] != "train":
#             continue
#         smi = df1["reactants_mol"][i] + "!" +  df1["products_mol"][i]
#         if smi not in set_50k:
#             dict_full["reactants_mol"].append(df1["reactants_mol"][i])
#             dict_full["products_mol"].append(df1["products_mol"][i])
#             dict_full["reaction_type"].append(df1["reaction_type"][i])
#             dict_full["set"].append(df1["set"][i])
#         else:
#             count += 1
#     print(count)
#     for i in range(df2.shape[0]):
#         if df2["set"][i] != "train":
#             dict_full["reactants_mol"].append(df2["reactants_mol"][i])
#             dict_full["products_mol"].append(df2["products_mol"][i])
#             dict_full["reaction_type"].append(df2["reaction_type"][i])
#             dict_full["set"].append(df2["set"][i])
    
#     df_full = pd.DataFrame(dict_full)
#     df_full.to_pickle("data/uspto_full_no50_1m.pickle")

# # detele_50k()

# def deal_pretrain_cross_graph_text():
#     import pandas as pd
#     df_full = pd.read_pickle("data/uspto_full_no50_1m.pickle")
#     df_50k = pd.read_pickle("data/uspto_50.pickle")

#     dict_pretrain = {
#         "smiles": [],
#         "atom_features": [],
#         "edges": [],
#         "adjacencys": [],
#         "node2smiles": [],
#         "set": [],
#     }
#     smiles_set = set()

#     for i in range(df_50k.shape[0]):
#         reactants_smi = Chem.MolToSmiles(df_50k["reactants_mol"][i])
#         products_smi = Chem.MolToSmiles(df_50k["products_mol"][i])
#         if reactants_smi not in smiles_set:
#             atom_feature, edge, adj, smile2node = get_graph_features_from_smi(reactants_smi)

#             dict_pretrain["smiles"].append(reactants_smi)
#             dict_pretrain["atom_features"].append(atom_feature)
#             dict_pretrain["edges"].append(edge)
#             dict_pretrain["adjacencys"].append(adj)
#             dict_pretrain["node2smiles"].append(smile2node)
#             dict_pretrain["set"].append(df_50k["set"][i])
#             smiles_set.add(reactants_smi)
        
#         if products_smi not in smiles_set:
#             atom_feature, edge, adj, smile2node = get_graph_features_from_smi(products_smi)

#             dict_pretrain["smiles"].append(products_smi)
#             dict_pretrain["atom_features"].append(atom_feature)
#             dict_pretrain["edges"].append(edge)
#             dict_pretrain["adjacencys"].append(adj)
#             dict_pretrain["node2smiles"].append(smile2node)
#             dict_pretrain["set"].append(df_50k["set"][i])
#             smiles_set.add(products_smi)
    
#     df = pd.DataFrame(dict_pretrain)
#     print(df.shape[0])
#     df.to_pickle("data/uspto_pretrain_50k.pickle")


# # deal_pretrain_cross_graph_text()
# # import pandas as pd
# # df = pd.read_pickle("data/uspto_full_no50_1m.pickle")
# # # for i in range(df.shape[0]):
# # #     if df["set"][i] == "train":
# # #         df["reactants_mol"][i] = Chem.MolFromSmiles(df["reactants_mol"][i])
# # #         df["products_mol"][i] = Chem.MolFromSmiles(df["products_mol"][i])
# # # df.to_pickle("data/uspto_full_no50_1m_mol.pickle")
# # print(df)

import pandas as pd

# df_50k = pd.read_pickle("data/uspto_50.pickle")
# df_templates = pd.read_csv("data/templates.csv")

# smi_set = dict()

# result_dict = {
#     "reactants_mol": [],
#     "products_mol": [],
#     "reaction_type": [],
#     "set": [],
#     "reactants": [],
#     "products": [],
# }

# print(df_templates.shape[0])
# for i in range(df_templates.shape[0]):
#     prod = function_1(df_templates["products"][i])
#     reac = function_1(df_templates["reactants"][i])
#     prod_smi = Chem.MolToSmiles(prod, canonical=False)
#     reac_smi = Chem.MolToSmiles(reac, canonical=False)

#     # smi = Chem.MolToSmiles(prod, canonical=True) + "!" + Chem.MolToSmiles(reac, canonical=True)
#     smi = Chem.MolToSmiles(reac, canonical=True)
#     smi_set[smi] = {
#         "products": prod_smi,
#         "reactants": reac_smi,
#     }
    
# # for smi in smi_set.keys():
# #     print(smi)

# count = 0
# for i in range(df_50k.shape[0]):
#     reac = Chem.MolToSmiles(df_50k["reactants_mol"][i], canonical=True)
#     prod = Chem.MolToSmiles(df_50k["products_mol"][i], canonical=True)
    
#     # smi = prod + "!" + reac
#     smi = reac
#     if smi in smi_set.keys():
#         result_dict["reactants"].append(smi_set[smi]["reactants"])
#         result_dict["products"].append(smi_set[smi]["products"])
#         count += 1
    
#     else:
#         result_dict["reactants"].append("")
#         result_dict["products"].append("")
        
#     result_dict["reactants_mol"].append(reac)
#     result_dict["products_mol"].append(prod)
#     result_dict["reaction_type"].append(df_50k["reaction_type"][i])
#     result_dict["set"].append(df_50k["set"][i])

# print(count)
# df = pd.DataFrame(result_dict)
# df.to_csv("data/uspto_50_templates.csv")
# # df.to_pickle("data/uspto_50_templates.pickle")


# smi = "[CH3:17][S:14](=[O:15])(=[O:16])[N:11]1[CH2:10][CH2:9][NH:8][CH2:13][CH2:12]1"
# mol = function_1(smi)

# print(Chem.MolToSmiles(mol, canonical=False))
# print(Chem.MolToSmiles(mol, canonical=True))

# df_50k = pd.read_csv("data/retro_uspto_50_template.csv")
# df_dict = {
#     "reactants_mol":[],
#     "products_mol": [],
#     "reaction_type":[],
#     "set":[],
#     "reactants":[],
#     "products":[],
# }

# for i in range(df_50k.shape[0]):
#     if df_50k["reactants"][i] == "" or df_50k["reactants"][i] is None or type(df_50k["reactants"][i]) == type(0.0):
#         df_dict["reactants"].append(df_50k["reactants_mol"][i])
#         df_dict["products"].append(df_50k["products_mol"][i])
#     else:
#         mol_reac = function_1(df_50k["reactants"][i])
#         mol_prod = function_1(df_50k["products"][i])
#         df_dict["reactants"].append(Chem.MolToSmiles(mol_reac, canonical=False))
#         df_dict["products"].append(Chem.MolToSmiles(mol_prod, canonical=False))
    
#     df_dict["products_mol"].append(df_50k["products_mol"][i])
#     df_dict["reactants_mol"].append(df_50k["reactants_mol"][i])
#     df_dict["reaction_type"].append(df_50k["reaction_type"][i])
#     df_dict["set"].append(df_50k["set"][i])
        
# df = pd.DataFrame(df_dict)
# df.to_pickle("data/retro_uspto_50_template_md.pickle")
# print(df)

# df = pd.read_csv("data/retro_uspto_50_template_md.csv")

# total = 0
# acc_token = 0
# acc_str = 0
# for i in range(df.shape[0]):
#     mol_str = df["products"][i]
#     target_str = df["products_mol"][i]
    
#     len_mol = len(mol_str) if len(mol_str) < len(target_str) else len(target_str)
#     for j in range(len_mol):
        
#         total += 1
        
#         if mol_str[j] == target_str[j]:
#             acc_token += 1
    
#     if mol_str == target_str:
#         acc_str += 1

# print("token_acc", acc_token / total)
# print("str_acc", acc_str / df.shape[0])


def smi_tokens(smi):
    """
    Tokenize a SMILES molecule or reaction
    """
    import re

    # pattern = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    # regex = re.compile(pattern)
    # tokens = [token for token in regex.findall(smi)]
    # assert smi == ''.join(tokens)
    # return tokens
    pattern = "(\[[^\]]+]|Bi|Br?|Ge|Te|Mo|K|Ti|Zr|Y|Na|125I|Al|Ce|Cr|Cl?|Ni?|O|S|Pd?|Fe?|I|b|c|Mn|n|o|s|<unk>|>>|Li|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    if smi != ''.join(tokens):
        print('ERROR:', smi, ''.join(tokens))
    assert smi == ''.join(tokens)
    return tokens


def mol_map_diff_smiles(smi1, smi2):
    from rdkit.Chem.rdFMCS import FindMCS

    src_chars = smi_tokens(smi1)
    tgt_chars = smi_tokens(smi2)

    src_mol = Chem.MolFromSmiles(smi1)
    tgt_mol = Chem.MolFromSmiles(smi2)

    atom_map = torch.zeros(src_mol.GetNumAtoms(), tgt_mol.GetNumAtoms())

    tgt_mol = Chem.MolFromSmiles(smi2)
    mols = [src_mol, tgt_mol]
    result = FindMCS(mols, timeout=10)
    result_mol = Chem.MolFromSmarts(result.smartsString)
    src_mat = src_mol.GetSubstructMatches(result_mol)
    print(src_mat[0])
    tgt_mat = tgt_mol.GetSubstructMatches(result_mol)
    if len(src_mat) > 0 and len(tgt_mat) > 0:
        for i, j in zip(src_mat[0], tgt_mat[0]):
            atom_map[i, j] = 1

    return atom_map

torch.set_printoptions(profile="full")

smi1 = "C([C@H]1CCC(=O)O1)N1CCN(CCOc2cc3ncnc(Nc4ccc(F)c(Cl)c4)c3cc2OC2CCCC2)CC1"
smi2 = "O=C1CC[C@H](CN2CCN(CCOc3cc4ncnc(Nc5ccc(F)c(Cl)c5)c4cc3OC3CCCC3)CC2)O1"

atom_map = mol_map_diff_smiles(smi1, smi2)
