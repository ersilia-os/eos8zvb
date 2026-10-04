import os
import sys
import random
import copy

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.inchi import MolToInchi

_ENGINE = None

_BENZENE_INCHI = "InChI=1S/C6H6/c1-2-4-6-5-3-1/h1-6H"
_MAX_ATTEMPTS = 5000


def _similarity(ref_fp, mol):
    """Morgan (radius 2, 2048 bits) Tanimoto between `ref_fp` and `mol`; 0.0 if there is no reference."""
    if ref_fp is None:
        return 0.0
    return DataStructs.TanimotoSimilarity(
        ref_fp, AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    )


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


def _fragment_inchi(frag_mol):
    """InChI of a fragment, or None. A fragment cut out of a larger molecule can have an aromatic
    ring nitrogen that lost its hydrogen with the substituent (purine, indole, quinolone, ...), and
    the plain conversion then fails to kekulize: a hydrogen is put back on one aromatic nitrogen at
    a time until it works (InChI's mobile-hydrogen normalisation makes the choice irrelevant)."""
    from pymolgen.molecule_formats import molecule_to_inchi, molecule_to_rdkit

    try:
        return molecule_to_inchi(frag_mol)
    except Exception:
        pass
    try:
        rdmol = molecule_to_rdkit(frag_mol)
    except Exception:
        return None
    for atom in rdmol.GetAtoms():
        if atom.GetSymbol() == "N" and atom.GetIsAromatic():
            repaired = Chem.RWMol(rdmol)
            repaired.GetAtomWithIdx(atom.GetIdx()).SetNumExplicitHs(1)
            try:
                return Chem.MolToInchi(repaired, options="-SNon")
            except Exception:
                continue
    return None


def _parent_candidates(smiles, checkpoints_dir, inchi_lookup):
    """Fragment ids to grow from, best first: the fragments of the input that are in the fragment
    database, largest first. Benzene is added as the last resort (it is the only candidate for an
    unparseable input or one with no fragment in the database)."""
    sys.path.insert(0, checkpoints_dir)
    from pymolgen.molecule_formats import molecule_from_smiles
    from pymolgen.fragment_mol import get_fragments_dataset
    from pymolgen.molecule import Molecule

    try:
        mol = molecule_from_smiles(smiles)
    except Exception:
        # molecule_from_smiles calls Chem.AddHs() unconditionally before
        # checking for a parse failure, so an invalid SMILES raises instead
        # of returning None. Treat it the same as an unparseable molecule.
        mol = None

    heavy_atoms = {}  # fragment id -> number of heavy atoms
    if mol is not None:
        try:
            frags, _, _ = get_fragments_dataset(mol)
        except Exception:
            frags = []
        for fg in frags or []:
            try:
                frag_mol = Molecule()
                frag_mol.graph = fg
                frag_id = inchi_lookup.get(_fragment_inchi(frag_mol))
            except Exception:
                continue
            if frag_id is not None:
                heavy_atoms[frag_id] = sum(
                    1 for n in fg.nodes if fg.nodes[n]["element"] != "H"
                )

    candidates = sorted(heavy_atoms, key=lambda frag_id: -heavy_atoms[frag_id])
    benzene_id = inchi_lookup.get(_BENZENE_INCHI, 0)
    if benzene_id not in candidates:
        candidates.append(benzene_id)
    return candidates


def _parent_molecule(frag_id, fragment_db_graph):
    """A FragmentMolecule made of fragment `frag_id` alone, or None if it has no attachment point."""
    from pymolgen.fragment_molecule import FragmentMolecule

    fragment = fragment_db_graph.fragments[frag_id]
    if not fragment.attachment_points:
        return None
    ap = fragment.attachment_points[0]
    cm = fragment.get_canonical_mapping()[ap]

    parent = FragmentMolecule()
    parent.add_fragment(frag_id, [ap], {ap: cm})
    parent._graph._build_probability2 = 1.0
    return parent


def _grow(parent, n, generated, seen, engine, ref_fp, input_flat):
    """Randomly extend `parent`, appending (tanimoto, canonical SMILES) to `generated` until it holds `n`
    molecules or _MAX_ATTEMPTS is used up. Every attempt is one independent random growth, and it adds
    one molecule: extend_molecule_random yields the growing molecule after each added fragment, and
    keeping all of those stages would fill the output with nested copies of the same growth, so one stage
    is picked at random. Duplicates (canonical isomeric SMILES, shared `seen` set) and echoes of the
    input (equal to `input_flat` ignoring stereochemistry) are skipped."""
    from pymolgen.fragment_molecule import convert_fragment_molecule_to_mol
    from pymolgen.fragment_molecule_builder import extend_molecule_random
    from pymolgen.molecule_formats import molecule_to_smiles

    fragment_db, fragment_db_graph, bond_freq_dict, _ = engine

    for _ in range(_MAX_ATTEMPTS):
        if len(generated) >= n:
            return
        stages = []
        try:
            for mol in extend_molecule_random(
                FragmentMolecule=parent,
                bond_frequencies=bond_freq_dict,
                fragment_database_graph=fragment_db_graph,
                depth=10,
                depth_min=3,
            ):
                stages.append(mol)
        except Exception:
            pass
        if not stages:
            continue
        try:
            mol_obj = convert_fragment_molecule_to_mol(random.choice(stages), fragment_db)
            smi = molecule_to_smiles(mol_obj)
            rdmol = Chem.MolFromSmiles(smi) if smi else None
            if rdmol is not None:
                key = Chem.MolToSmiles(rdmol)
                if (
                    key not in seen
                    and Chem.MolToSmiles(rdmol, isomericSmiles=False) != input_flat
                ):
                    seen.add(key)
                    generated.append((_similarity(ref_fp, rdmol), key))
        except Exception:
            pass


def generate_analogues(smiles, checkpoints_dir, n=100):
    """Generate n drug-like analogues of input SMILES, ordered from most to least similar
    (Morgan Tanimoto) to the input. Returns list of length n (None for failed slots, at the end).

    The largest fragment of the input found in the fragment database is grown by random fragment
    additions, one independent growth per output; if it cannot give n molecules, the next largest
    fragment is used, and finally benzene. Outputs are canonical SMILES, unique, and never equal to
    the input."""
    engine = _get_engine(checkpoints_dir)
    fragment_db_graph, inchi_lookup = engine[1], engine[3]

    # an unparseable input has no reference: no similarity order and no echo check
    input_mol = Chem.MolFromSmiles(smiles)
    if input_mol is not None:
        ref_fp = AllChem.GetMorganFingerprintAsBitVect(input_mol, 2, nBits=2048)
        input_flat = Chem.MolToSmiles(input_mol, isomericSmiles=False)
    else:
        ref_fp, input_flat = None, None

    generated = []
    seen = set()
    for frag_id in _parent_candidates(smiles, checkpoints_dir, inchi_lookup):
        if len(generated) >= n:
            break
        parent = _parent_molecule(frag_id, fragment_db_graph)
        if parent is not None:
            _grow(parent, n, generated, seen, engine, ref_fp, input_flat)

    # most similar first; the sort is stable, so ties keep their generation order
    generated.sort(key=lambda item: -item[0])
    generated = [smi for _, smi in generated]

    while len(generated) < n:
        generated.append(None)

    return generated[:n]
