from rdkit import Chem
from functools import partial


def _fix_charge(mol, match, value=1):
    """
    Sets charge of the first atom in the match to specified value
    """
    mol.GetAtomWithIdx(match[0]).SetFormalCharge(value)
    return mol


def fix_charge(value):
    """
    Generates the function to set charge of the first atom in the 
    match to specified value
    """
    return partial(_fix_charge, value=value)


def set_dative_bonds(mol, match, all=False, bond_type=Chem.BondType.DATIVE):
    """
    Sets a bond between first to matched atoms to specified bond type

    The bond is created to preserve the order of the atoms in the matche
    Ie. The dative bond will have first atom as donor and second as acceptor
    """

    emol = Chem.RWMol(mol)
    atom = emol.GetAtomWithIdx(match[0])
    if all:
        bonds = atom.GetBonds()
    else:
        bonds = [emol.GetBondBetweenAtoms(match[0], match[1])]

    for bond in bonds:
        idx1, idx2 = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        emol.RemoveBond(idx1, idx2)

        # Sort bond indexes so our match is second (acceptor) (direction important)
        if match[0] != idx2:
            idx2, idx1 = idx1, idx2

        emol.AddBond(idx1, idx2, bond_type)
        atom = emol.GetAtomWithIdx(idx2)

    return emol


def set_unspecified_bonds(mol, match, all=True):
    """
    Reset the information about bond type
    """
    atom = mol.GetAtomWithIdx(match[0])
    if all:
        bonds = atom.GetBonds()
    else:
        bonds = [mol.GetBondBetweenAtoms(match[0], match[1])]

    for bond in bonds:
        bond.SetBondType(Chem.BondType.UNSPECIFIED)
    return mol


def skip(mol, match):
    return False


fixes = [
    # (fix_name, RdKit SMARTS pattern, fix_function -> fixed molecule )
    # Fixes are applied in the spcified order
    # If the function returns None the fix was unsuccessful
    # If reurns False then it is deliberetly skipped
    (
        "[skp] borane clusters",
        Chem.MolFromSmarts("B(B)(*)(*)"),
        skip,
    ),  # Creates problems with clustering, very few compounds
    ("[fix] tertiary amine", Chem.MolFromSmarts("N(C)(C)C"), fix_charge(1)),
    ("[fix] nitrogen dioxide", Chem.MolFromSmarts("N(=O)(-O)*"), fix_charge(1)),
    ("[fix] group 1 charge", Chem.MolFromSmarts("[Li,Na,K,Rb,Cs]"), fix_charge(1)),
    ("[fix] group 2 charge", Chem.MolFromSmarts("[Be,Mg,Ca,Sr,Ba]"), fix_charge(2)),
    (
        "[wka] metal coord. bonds",
        Chem.MolFromSmarts(
            "[Na,Mg,Al,K,Ca,Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,"
            "Rb,Sr,Y,Zr,Nb,Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,"
            "W,Re,Os,Ir,Pt,Au,Hg]~[O,N,C]"
        ),
        partial(set_unspecified_bonds, all=False),
    ),
    ("[fix] azide charge", Chem.MolFromSmarts("N(=N)=N"), fix_charge(1)),
    (
        "[wka] chlorine dioxide",
        Chem.MolFromSmarts("Cl(=O)(=O)*"),
        set_unspecified_bonds,
    ),
    ("[wka] CuClCu", Chem.MolFromSmarts("[Cl]([Cu])[Cu]"), set_unspecified_bonds),
    ("[wka] BeF3", Chem.MolFromSmarts("[Be](F)(F)F"), set_unspecified_bonds),
    (
        "[wka] ME-C#O",
        Chem.MolFromSmarts(
            "O#C~[Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn,Ga,Rb,Sr,Y,Zr,Nb,"
            "Mo,Tc,Ru,Rh,Pd,Ag,Cd,In,Sn,Hf,Ta,W,Re,Os,Ir,Pt,Au,Hg]"
        ),
        partial(set_unspecified_bonds, all=False),
    ),
]
