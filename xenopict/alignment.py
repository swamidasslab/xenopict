"""Functions for aligning 2D molecular depictions using RDKit.

This module provides several methods for aligning molecules to template structures:
- Basic alignment using automatic substructure matching
- Manual alignment using explicit atom pairs
- Alignment using atom indices
- Alignment using atom map IDs
"""

from rdkit import Chem  # type: ignore
from rdkit.Chem import rdDepictor  # type: ignore
from typing import List, Dict, Tuple, Union
from rdkit.Chem.rdFMCS import FindMCS
from typing import NamedTuple, Optional


class Alignment(NamedTuple):
    aligned_atoms: list[tuple[int, int]]
    score: float
    from_mol: Chem.Mol
    to_mol: Chem.Mol

    def reverse(self) -> "Alignment":
        return Alignment(
            aligned_atoms=[(b, a) for a, b in self.aligned_atoms],
            score=self.score,
            from_mol=self.to_mol,
            to_mol=self.from_mol,
        )


def auto_align_molecules(
    mols: list[Chem.Mol], hints: Optional[list[Alignment]] = None
) -> list[Chem.Mol]:
    """
    Align a list of molecules using automatic alignment and using networkx to find a maximum
    spanning tree of the alignment graph.

    Args:
        mols: list of molecules to align
        hints: list of Alignment objects to use as hints for the alignment
    """
    from networkx import maximum_spanning_tree, Graph
    import  networkx as nx

    hints = hints or []

    g = Graph()
    mol2node = {}
    for i, mol in enumerate(mols):
        g.add_node(i, mol=mol)
        mol2node[id(mol)] = i

    hinted_edges = set()
    for alignment in hints:
        i = mol2node[id(alignment.from_mol)]
        j = mol2node[id(alignment.to_mol)]
        hinted_edges.add(frozenset({i, j}))

        # add hinted alignments graph, with  very high weight so it is preferred
        g.add_edge(i, j, weight=alignment.score + 1000, alignment=alignment, template=j)


    print(hinted_edges)
    # iterate over all pairs
    for i, mol in enumerate(mols):
        for j, other_mol in enumerate(mols):
            # if we already considered this pair
            if i >= j:
                continue
              
            # don't auto align pairs that are aligned in hints
            if frozenset({i, j}) in hinted_edges:
              continue

            alignment = auto_alignment(mol, other_mol)

            g.add_edge(i, j, weight=alignment.score, alignment=alignment, template=j)

    # pick the alignments based on the maximum spanning tree
    t = maximum_spanning_tree(g)

    # unrooted tree so any node can be the root
    # we will pick the largest molecule
    max_mol_size = -1
    max_mol_idx = -1
    for node in t.nodes:
        m = t.nodes[node]["mol"]
        s = m.GetNumAtoms()
        if s > max_mol_size:
            max_mol_size = s
            max_mol_idx = node

    # now iterate over the tree in BFS order
    for edge in nx.traversal.bfs_edges(t, max_mol_idx):
        alignment = g.edges[edge]["alignment"]

        # if the first idx of the edge isn't the template swap the alignment
        if edge[0] != g.edges[edge]["template"]:
            alignment = alignment.reverse()

        template = g.nodes[edge[0]]["mol"]
        mol = g.nodes[edge[1]]["mol"]

        align_to_template_manual(mol, template, alignment.aligned_atoms)

    return [g.nodes[i]["mol"] for i in sorted(g.nodes)]


def align_to_template(mol: Chem.Mol, template: Chem.Mol) -> Chem.Mol:
    """Align a molecule to a template molecule using automatic substructure matching.

    This function attempts to align the input molecule to match the 2D depiction
    of the template molecule. The alignment is done automatically using RDKit's
    substructure matching.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to

    Returns:
        The aligned molecule (same object as input mol)
    """

    alignment = auto_alignment(mol, template)

    if not alignment.aligned_atoms:
        raise ValueError("No alignment found")
    return align_to_template_manual(mol, template, alignment.aligned_atoms)


def auto_alignment(mol: Chem.Mol, template: Chem.Mol) -> Alignment:
    """Align a molecule to a template molecule using maximum common substructure.

    Args:
        mol: The molecule to align
        template: The template molecule to align to

    Returns:
        An Alignment namedtuple containing the aligned atoms, score, and the two molecules
    """

    exact_mcs = FindMCS([mol, template])
    exact_query = exact_mcs.queryMol

    inexact_mcs = FindMCS(
        [mol, template], bondCompare=Chem.rdFMCS.BondCompare.CompareAny
    )
    inexact_query = inexact_mcs.queryMol

    score = (
        exact_mcs.numAtoms
        + exact_mcs.numBonds
        + inexact_mcs.numAtoms
        + inexact_mcs.numBonds
    ) / 4

    if exact_mcs.numAtoms == inexact_mcs.numAtoms:
        query = exact_query
    else:
        query = inexact_query

    if score == 0:
        return Alignment([], 0, mol, template)

    aligned_atoms = list(
        zip(
            Chem.Mol(mol).GetSubstructMatch(query),
            Chem.Mol(template).GetSubstructMatch(query),
        )
    )

    return Alignment(aligned_atoms, score, mol, template)


def align_to_template_manual(
    mol: Chem.Mol, template: Chem.Mol, aligned_atoms: List[Tuple[int, int]]
) -> Chem.Mol:
    """Align a molecule to a template molecule using explicit atom pairs.

    This function aligns the input molecule to match the 2D depiction of the template
    molecule using explicitly specified atom pairs.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        aligned_atoms: List of (mol_atom_idx, template_atom_idx) pairs specifying
                      which atoms should be aligned to each other

    Returns:
        The aligned molecule (same object as input mol)

    Raises:
        AssertionError: If any atom indices are out of range
    """
    _ensure_coords(template)

    # Validate atom indices are within bounds
    for mol_idx, template_idx in aligned_atoms:
        assert mol_idx < mol.GetNumAtoms(), (
            f"Molecule atom index {mol_idx} out of range (max {mol.GetNumAtoms() - 1})"
        )
        assert template_idx < template.GetNumAtoms(), (
            f"Template atom index {template_idx} out of range (max {template.GetNumAtoms() - 1})"
        )
        assert mol_idx >= 0, f"Negative molecule atom index {mol_idx}"
        assert template_idx >= 0, f"Negative template atom index {template_idx}"

    # bizarrely, convention of GenerateDepictionMatching2DStructure is opposite order
    aligned_atoms = [(t, m) for m, t in aligned_atoms]

    if len(aligned_atoms):  # empty list causes segfault
        rdDepictor.GenerateDepictionMatching2DStructure(mol, template, aligned_atoms)
    else:
        _ensure_coords(mol)
    return mol


def align_to_template_with_indices(
    mol: Chem.Mol, template: Chem.Mol, mol_to_template: Union[List[int], Dict[int, int]]
) -> Chem.Mol:
    """Align a molecule to a template molecule using atom index mappings.

    This function aligns the input molecule to match the 2D depiction of the template
    molecule using a mapping between atom indices.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        mol_to_template: Either:
            - A list where index i contains the template atom index that mol atom i maps to
              (-1 for unmapped atoms). Length must match number of atoms in template.
            - A dict mapping mol atom indices to template atom indices
              (omit or use -1 for unmapped atoms)

    Returns:
        The aligned molecule (same object as input mol)

    Raises:
        ValueError: If mol_to_template is not a list or dict, or if indices are out of range
        AssertionError: If mol_to_template is a list with wrong length
    """
    _ensure_coords(template)

    if isinstance(mol_to_template, list):
        assert len(mol_to_template) == mol.GetNumAtoms(), (
            f"Length of mol_to_template must equal number of mol atoms: "
            f"{len(mol_to_template)} != {mol.GetNumAtoms()}"
        )

        # Convert list format to explicit atom pairs, skipping unmapped atoms (-1)
        aligned_atoms = []
        for mol_idx, template_idx in enumerate(mol_to_template):
            if template_idx >= 0:
                if template_idx >= template.GetNumAtoms():
                    raise ValueError(
                        "Reference atom index in mol_to_template out of range"
                    )
                aligned_atoms.append((mol_idx, template_idx))

    elif isinstance(mol_to_template, dict):
        # Convert dict format to explicit atom pairs, skipping unmapped atoms (-1)
        aligned_atoms = []
        for mol_idx, template_idx in mol_to_template.items():
            if template_idx < -1:  # Allow -1 for unmapped atoms
                raise ValueError("Reference atom index in mol_to_template out of range")
            if template_idx >= template.GetNumAtoms():
                raise ValueError("Reference atom index in mol_to_template out of range")
            if template_idx >= 0:
                aligned_atoms.append((mol_idx, template_idx))
    else:
        raise ValueError(
            f"mol_to_template must be a list or dictionary, got {type(mol_to_template)}"
        )

    return align_to_template_manual(mol, template, aligned_atoms)


def align_to_template_by_mapids(
    mol: Chem.Mol, template: Chem.Mol, mol_map: Chem.Mol, template_map: Chem.Mol
) -> Chem.Mol:
    """Align a molecule to a template using atom map IDs.

    This function aligns the input molecule to match the 2D depiction of the template
    molecule using atom map IDs to identify corresponding atoms. The mapping is done
    through two auxiliary molecules (mol_map and template_map) that contain the map IDs.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        mol_map: Molecule with map IDs matching atoms in mol
        template_map: Molecule with map IDs matching atoms in template

    Returns:
        The aligned molecule (same object as input mol)
    """
    _ensure_coords(template)

    # Get map ID to atom index mappings for both molecules
    mapid2idx_mol = _get_matched_mapids(mol, mol_map)
    mapid2idx_template = _get_matched_mapids(template, template_map)

    # Create list of aligned atom pairs from matching map IDs
    aligned_atoms = [
        (mapid2idx_mol[map_id], mapid2idx_template[map_id])
        for map_id in set(mapid2idx_mol) & set(mapid2idx_template)
    ]

    return align_to_template_manual(mol, template, aligned_atoms)


def GetCoords(mol: Chem.Mol, i: int) -> tuple[float, float, float]:
    c = mol.GetConformer(0).GetAtomPosition(i)
    return c.x, c.y, c.z


def _ensure_coords(mol: Chem.Mol) -> Chem.Mol:
    """Ensure that the molecule has 2D coordinates.

    Args:
        mol: The molecule to check/compute coordinates for

    Returns:
        The input molecule (modified in place with 2D coords if needed)
    """
    if mol.GetNumConformers() == 0:
        rdDepictor.Compute2DCoords(mol)
    return mol


def _extract_map_ids(mol: Chem.Mol) -> Dict[int, int]:
    """Extract atom map IDs from a molecule.

    Args:
        mol: The molecule to extract map IDs from

    Returns:
        Dictionary mapping atom map IDs to atom indices

    Raises:
        ValueError: If mol is None (invalid SMILES)
    """
    if mol is None:
        raise ValueError("Invalid SMILES")

    return {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() != 0
    }


def _get_matched_mapids(mol: Chem.Mol, mol_map: Chem.Mol) -> Dict[int, int]:
    """Get mapping between atom map IDs and atom indices for matched atoms.

    Args:
        mol: The molecule to match against
        mol_map: The molecule containing the atom mapping IDs to match

    Returns:
        Dictionary mapping atom map IDs to matched atom indices in mol

    Raises:
        ValueError: If mol_map is None or no substructure match is found
    """
    # Get map ID to index mapping from the mapping molecule
    mapid2idx = _extract_map_ids(mol_map)

    # Find substructure match between molecules
    match = mol.GetSubstructMatch(mol_map)
    if not match:
        # If no match found, just ensure coords and return empty mapping
        _ensure_coords(mol)
        return {}

    # Map the map IDs to the matched atom indices
    return {
        map_id: match[idx]
        for map_id, idx in mapid2idx.items()
        if idx < len(match) and match[idx] >= 0
    }
