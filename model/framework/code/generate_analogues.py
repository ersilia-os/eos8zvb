import os
import sys
import random
import copy

from rdkit import Chem
from rdkit.Chem.inchi import MolToInchi

_ENGINE = None


def _load_engine(checkpoints_dir):
    sys.path.insert(0, checkpoints_dir)

    from pymolgen.fragment_molecule import convert_fragment_molecule_to_mol
    from pymolgen.fragment_molecule_builder import (
        read_fragment_database_graph,
        read_bond_frequencies_dict,
    )
    from pymolgen.fragment_builder import get_fragment_database

    chembl = os.path.join(checkpoints_dir, "chembl")
    pfx = "fragments_30_50k_co_10_l5_5_sorted_filter_copy"

    fragment_db = get_fragment_database(os.path.join(chembl, f"{pfx}.sdf"))
    fragment_db_graph = read_fragment_database_graph(
        os.path.join(chembl, "fragment_database_30_50k_co_10_l5_5_sorted_filter_copy.txt")
    )
    bond_freq_dict = read_bond_frequencies_dict(
        os.path.join(chembl, "bond_frequencies_30_50k_co_10_l5_5_sorted_filter_copy.txt")
    )

    inchi_lookup = {}
    with open(os.path.join(chembl, f"{pfx}.inchi")) as fh:
        for frag_id, line in enumerate(fh):
            inchi = line.strip()
            if inchi and inchi not in inchi_lookup:
                inchi_lookup[inchi] = frag_id

    return fragment_db, fragment_db_graph, bond_freq_dict, inchi_lookup


def _get_engine(checkpoints_dir):
    global _ENGINE
    if _ENGINE is None:
        _ENGINE = _load_engine(checkpoints_dir)
    return _ENGINE


def _find_parent_frag_id(smiles, checkpoints_dir, fragment_db_graph, inchi_lookup):
    """Return (frag_id, attach_atom) for the best matching parent fragment."""
    sys.path.insert(0, checkpoints_dir)
    from pymolgen.molecule_formats import molecule_from_smiles, molecule_to_inchi
    from pymolgen.fragment_mol import get_fragments_dataset
    from pymolgen.molecule import Molecule

    mol = molecule_from_smiles(smiles)
    best_frag_id, best_n_heavy = None, 0

    if mol is not None:
        try:
            frags, _, _ = get_fragments_dataset(mol)
            for fg in (frags or []):
                try:
                    frag_mol = Molecule()
                    frag_mol.graph = fg
                    frag_mol.free_valence_list = []
                    inchi = molecule_to_inchi(frag_mol)
                    frag_id = inchi_lookup.get(inchi)
                    if frag_id is None:
                        continue
                    n_heavy = sum(
                        1 for n in fg.nodes if fg.nodes[n]["element"] != "H"
                    )
                    if n_heavy > best_n_heavy:
                        best_n_heavy = n_heavy
                        best_frag_id = frag_id
                except Exception:
                    continue
        except Exception:
            pass

    if best_frag_id is None:
        # fallback: benzene
        benzene_inchi = "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
        best_frag_id = inchi_lookup.get(benzene_inchi, 0)

    fragment = fragment_db_graph.fragments[best_frag_id]
    if not fragment.attachment_points:
        best_frag_id = 0
        fragment = fragment_db_graph.fragments[0]

    ap = fragment.attachment_points[0]
    cm = fragment.get_canonical_mapping()[ap]
    return best_frag_id, ap, cm


def generate_analogues(smiles, checkpoints_dir, n=100):
    """Generate n drug-like analogues of input SMILES. Returns list of length n (None for failed slots)."""
    from pymolgen.fragment_molecule import FragmentMolecule, convert_fragment_molecule_to_mol
    from pymolgen.fragment_molecule_builder import extend_molecule_random
    from pymolgen.molecule_formats import molecule_to_smiles

    fragment_db, fragment_db_graph, bond_freq_dict, inchi_lookup = _get_engine(
        checkpoints_dir
    )

    frag_id, ap, cm = _find_parent_frag_id(
        smiles, checkpoints_dir, fragment_db_graph, inchi_lookup
    )

    parent = FragmentMolecule()
    parent.add_fragment(frag_id, [ap], {ap: cm})
    parent._graph._build_probability2 = 1.0

    generated = []
    seen = set()
    max_attempts = 5000

    for _ in range(max_attempts):
        if len(generated) >= n:
            break
        try:
            for mol in extend_molecule_random(
                FragmentMolecule=parent,
                bond_frequencies=bond_freq_dict,
                fragment_database_graph=fragment_db_graph,
                depth=10,
                depth_min=3,
            ):
                try:
                    mol_obj = convert_fragment_molecule_to_mol(mol, fragment_db)
                    smi = molecule_to_smiles(mol_obj)
                    if smi and smi not in seen:
                        rdmol = Chem.MolFromSmiles(smi)
                        if rdmol is not None:
                            seen.add(smi)
                            generated.append(smi)
                except Exception:
                    pass
                if len(generated) >= n:
                    break
        except Exception:
            pass

    while len(generated) < n:
        generated.append(None)

    return generated[:n]
