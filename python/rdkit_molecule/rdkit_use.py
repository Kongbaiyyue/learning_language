from rdkit import Chem
from rdkit.Chem import Draw

import re

from rdchiral.template_extractor import extract_from_reaction, get_changed_atoms, mols_from_smiles_list, \
    replace_deuterated, get_tetrahedral_atoms


# 去掉原子编号
def function_1():
    # 将带有原子编号的SMILES转换为分子对象
    mol = Chem.MolFromSmiles('c1ccc(C[N:8]2[CH2:9][CH2:10][N:11]([S:14](=[O:15])(=[O:16])[CH3:17])[CH2:12][CH2:13]2)cc1')
    img = Draw.MolToImage(mol)
    img.show()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)

    # 将分子对象转换为没有原子编号的SMILES
    smiles = Chem.MolToSmiles(mol, canonical=True)
    
    print(smiles)
    
    img.save("mol.jpg")

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
    print(changed_atom_tags)

    for _ in range(radius):
        mol = Chem.MolFromSmiles(cano_prod_am)
        changed_atom_tags_neighbor = []
        for atom in mol.GetAtoms():
            if atom.GetSmarts().split(':')[1][:-1] in changed_atom_tags:
                for n_atom in atom.GetNeighbors():
                    changed_atom_tags_neighbor.append(n_atom.GetSmarts().split(':')[1][:-1])
        changed_atom_tags = list(set(changed_atom_tags + changed_atom_tags_neighbor))
    # print(changed_atom_tags)

    nonreactive_mask = []
    for i, token in enumerate(smi_tokenizer(cano_prod_am)):
        if token[0] == '[' and token[-1] == ']' and re.match('.*:([0-9]+)]', token):
            am = re.match('.*:([0-9]+)]', token).group(1)
            if am in changed_atom_tags:
                nonreactive_mask.append(False)
                continue
        nonreactive_mask.append(True)

    if sum(nonreactive_mask) == len(nonreactive_mask):  # if the reaction center is not detected
        nonreactive_mask = [False] * len(nonreactive_mask)
    # print(nonreactive_mask)
    return nonreactive_mask


# smi = "[CH3:1][C:2](=[O:3])[O:4][C:5]([CH3:6])=[O:7].[NH2:8][c:9]1[cH:10][cH:11][cH:12][cH:13][c:14]1[F:15]"
# mol = Chem.MolFromSmiles(smi)
# smi2 = "[CH3:1][C:2](=[O:3])[NH:4][c:5]1[cH:6][cH:7][cH:8][cH:9][c:10]1[F:11]"
# mol2 = Chem.MolFromSmiles(smi2)
# # smi2 = "CS(=O)(=O)OC[C@H]1CCC(=O)O1.Fc1ccc(Nc2ncnc3cc(OCCN4CCNCC4)c(OC4CCCC4)cc23)cc1Cl"
# print(len(smi_tokenizer(smi2)))
# # smi2 = "CS(=O)(=O)O[CH2:1][C@H:2]1[CH2:3][CH2:4][C:5](=[O:6])[O:7]1.[F:8][c:9]1[cH:10][cH:11][c:12]([NH:13][c:14]2[n:15][cH:16][n:17][c:18]3[cH:19][c:20]([O:21][CH2:22][CH2:23][N:24]4[CH2:25][CH2:26][NH:27][CH2:28][CH2:29]4)[c:30]([O:31][CH:32]4[CH2:33][CH2:34][CH2:35][CH2:36]4)[cH:37][c:38]23)[cH:39][c:40]1[Cl:41]"
# # print(len(smi_tokenizer(smi2)))
# print(get_nonreactive_mask(smi2, smi2, smi))
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

reaction_smi = "[NH2:3][c:4]1[cH:5][cH:6][cH:7][c:8]2[cH:9][n:10][cH:11][cH:12][c:13]12.[O:1]=[C:2]([c:14]1[cH:15][c:16]([N+:17](=[O:18])[O-:19])[c:20]([S:21][c:22]2[c:23]([Cl:24])[cH:25][n:26][cH:27][c:28]2[Cl:29])[s:30]1)[OH:31]>>[O:1]=[C:2]([NH:3][c:4]1[cH:5][cH:6][cH:7][c:8]2[cH:9][n:10][cH:11][cH:12][c:13]12)[c:14]1[cH:15][c:16]([N+:17](=[O:18])[O-:19])[c:20]([S:21][c:22]2[c:23]([Cl:24])[cH:25][n:26][cH:27][c:28]2[Cl:29])[s:30]1"
r, p = reaction_smi.split(">>")
rxn_core, core_edits = get_reaction_core(r=r, p=p)

print(rxn_core)
print(core_edits)

reactants = mols_from_smiles_list(replace_deuterated(r).split('.'))
products = mols_from_smiles_list(replace_deuterated(p).split('.'))
changed_atoms, changed_atom_tags, err = get_changed_atoms(reactants, products)
print(changed_atom_tags)
# mol = Chem.MolFromSmiles(r)
# img = Draw.MolToImage(mol)
# img.show()

# mol = Chem.MolFromSmiles(p)
# img = Draw.MolToImage(mol)
# img.show()