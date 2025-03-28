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
    atom_pairs: list[tuple[int, int]]
    score: float
    source_mol: Chem.Mol
    template_mol: Chem.Mol

    @staticmethod
    def from_aligned_atoms(
        source_mol: Chem.Mol, template_mol: Chem.Mol, atom_pairs: list[tuple[int, int]]
    ) -> "Alignment":
        return Alignment(
            atom_pairs=atom_pairs,
            score=0,
            source_mol=source_mol,
            template_mol=template_mol,
        ).validate()

    @staticmethod
    def from_mcs(source_mol: Chem.Mol, template_mol: Chem.Mol) -> "Alignment":
        """Create alignment using maximum common substructure (MCS) matching.
        
        Uses both exact and inexact MCS matching to find the best alignment between
        molecules. The inexact matching allows for different bond types.
        """
        exact_mcs = FindMCS([source_mol, template_mol])
        exact_query = exact_mcs.queryMol

        inexact_mcs = FindMCS(
            [source_mol, template_mol], bondCompare=Chem.rdFMCS.BondCompare.CompareAny
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
            return Alignment([], 0, source_mol, template_mol)

        atom_pairs = list(
            zip(
                Chem.Mol(source_mol).GetSubstructMatch(query),
                Chem.Mol(template_mol).GetSubstructMatch(query),
            )
        )

        return Alignment(atom_pairs, score, source_mol, template_mol)

    @staticmethod
    def from_mapids(
        source_mol: Chem.Mol, template_mol: Chem.Mol, source_map: Chem.Mol, template_map: Chem.Mol
    ) -> "Alignment":
        """Create alignment using atom map IDs.
        
        Args:
            source_mol: The source molecule to align
            template_mol: The template molecule to align to
            source_map: The source molecule with atom map IDs
            template_map: The template molecule with atom map IDs
            
        Returns:
            An Alignment object
            
        Raises:
            ValueError: If no atom maps are found in either molecule or no matching map IDs exist
        """
        mapid2idx_source = _get_matched_mapids(source_mol, source_map, is_template=False)
        mapid2idx_template = _get_matched_mapids(template_mol, template_map, is_template=True)

        # Find common map IDs between source and template
        common_map_ids = set(mapid2idx_source) & set(mapid2idx_template)
        if not common_map_ids:
            raise ValueError("No matching atom map IDs found")

        atom_pairs = [
            (mapid2idx_source[map_id], mapid2idx_template[map_id])
            for map_id in common_map_ids
        ]

        return Alignment.from_aligned_atoms(source_mol, template_mol, atom_pairs)

    @staticmethod
    def from_indices(
        source_mol: Chem.Mol, template_mol: Chem.Mol, index_mapping: list[int] | dict[int, int]
    ) -> "Alignment":
        """Align a molecule to a template molecule using atom index mappings."""

        if not isinstance(index_mapping, dict) and not isinstance(
            index_mapping, list
        ):
            raise ValueError(
                f"index_mapping must be a list or dictionary, got {type(index_mapping)}"
            )

        # convert list to dict format
        if isinstance(index_mapping, list):
            assert len(index_mapping) == source_mol.GetNumAtoms(), (
                f"Length of index_mapping must equal number of source_mol atoms: "
                f"{len(index_mapping)} != {source_mol.GetNumAtoms()}"
            )

            index_mapping = {k: v for k, v in enumerate(index_mapping) if v >= 0}

        # convert dict format to explicit atom pairs
        atom_pairs = list(index_mapping.items())

        return Alignment.from_aligned_atoms(source_mol, template_mol, atom_pairs)

    def reverse(self) -> "Alignment":
        return Alignment(
            atom_pairs=[(b, a) for a, b in self.atom_pairs],
            score=self.score,
            source_mol=self.template_mol,
            template_mol=self.source_mol,
        )

    def validate(self) -> "Alignment":
        source_mol = self.source_mol
        template_mol = self.template_mol
        atom_pairs = self.atom_pairs

        # Validate atom indices are within bounds
        for source_idx, template_idx in atom_pairs:
            assert source_idx < source_mol.GetNumAtoms(), (
                f"Source molecule atom index {source_idx} out of range (max {source_mol.GetNumAtoms() - 1})"
            )
            assert template_idx < template_mol.GetNumAtoms(), (
                f"Template atom index {template_idx} out of range (max {template_mol.GetNumAtoms() - 1})"
            )
            assert source_idx >= 0, f"Negative source molecule atom index {source_idx}"
            assert template_idx >= 0, f"Negative template atom index {template_idx}"

        return self

    def apply(self) -> "Alignment":
        """Apply the alignment to the molecule so that the source_mol coordinates are changed to
        match the template_mol coordinates."""

        source_mol = self.source_mol
        template_mol = self.template_mol
        atom_pairs = self.atom_pairs

        _ensure_coords(template_mol)

        # bizarrely, convention of GenerateDepictionMatching2DStructure is opposite order
        atom_pairs = [(t, m) for m, t in atom_pairs]

        if len(atom_pairs):  # empty list causes segfault
            rdDepictor.GenerateDepictionMatching2DStructure(
                source_mol, template_mol, atom_pairs
            )
        else:
            _ensure_coords(source_mol)

        return self


def auto_align_molecules(
    molecules: list[Chem.Mol], hints: Optional[list[Alignment]] = None
) -> list[Chem.Mol]:
    """Align a list of molecules using automatic alignment and using networkx to find a maximum
    spanning tree of the alignment graph.

    Args:
        molecules: list of molecules to align
        hints: list of Alignment objects to use as hints for the alignment
    """
    from networkx import maximum_spanning_tree, Graph
    import networkx as nx

    hints = hints or []

    g = Graph()
    mol2node = {}
    for i, mol in enumerate(molecules):
        g.add_node(i, mol=mol)
        mol2node[id(mol)] = i

    hinted_edges = set()
    for alignment in hints:
        i = mol2node[id(alignment.source_mol)]
        j = mol2node[id(alignment.template_mol)]
        hinted_edges.add(frozenset({i, j}))

        # add hinted alignments graph, with very high weight so it is preferred
        g.add_edge(i, j, weight=alignment.score + 1000, alignment=alignment, template=j)

    # iterate over all pairs
    for i, source_mol in enumerate(molecules):
        for j, template_mol in enumerate(molecules):
            # if we already considered this pair
            if i >= j:
                continue

            # don't auto align pairs that are aligned in hints
            if frozenset({i, j}) in hinted_edges:
                continue

            alignment = Alignment.from_mcs(source_mol, template_mol)
            g.add_edge(i, j, weight=alignment.score, alignment=alignment, template=j)

    # pick the alignments based on the maximum spanning tree
    t = maximum_spanning_tree(g)

    # unrooted tree so any node can be the root
    # we will pick the largest molecule
    max_mol_size = -1
    max_mol_idx = -1
    for node in t.nodes:
        m = t.nodes[node]["mol"]  # type: ignore
        s = m.GetNumAtoms()
        if s > max_mol_size:
            max_mol_size = s
            max_mol_idx = node

    # now iterate over the tree in BFS order
    for edge in nx.traversal.bfs_edges(t, max_mol_idx):  # type: ignore
        alignment = g.edges[edge]["alignment"]

        # if the first idx of the edge isn't the template swap the alignment
        if edge[0] != g.edges[edge]["template"]:
            alignment = alignment.reverse()

        alignment.apply()

    return [g.nodes[i]["mol"] for i in sorted(g.nodes)]


def align_from_mcs(source_mol: Chem.Mol, template_mol: Chem.Mol) -> Chem.Mol:
    """Align a molecule to a template molecule using maximum common substructure (MCS).
    
    This is the main high-level function for automatic alignment. It finds the largest
    common substructure between the molecules and uses it to determine the alignment.
    Both exact and inexact matching are used to find the best possible alignment.
    """
    alignment = Alignment.from_mcs(source_mol, template_mol)
    if not alignment.atom_pairs:
        raise ValueError("No alignment found")
    alignment.apply()
    return source_mol


def align_from_atom_pairs(source_mol: Chem.Mol, template_mol: Chem.Mol, atom_pairs: List[Tuple[int, int]]) -> Chem.Mol:
    """Align a molecule to a template molecule using explicit atom pairs."""
    alignment = Alignment.from_aligned_atoms(source_mol, template_mol, atom_pairs)
    alignment.apply()
    return source_mol


def align_from_indices(source_mol: Chem.Mol, template_mol: Chem.Mol, index_mapping: Union[List[int], Dict[int, int]]) -> Chem.Mol:
    """Align a molecule to a template molecule using atom index mappings."""
    alignment = Alignment.from_indices(source_mol, template_mol, index_mapping)
    alignment.apply()
    return source_mol


def align_from_mapids(source_mol: Chem.Mol, template_mol: Chem.Mol, source_map: Chem.Mol, template_map: Chem.Mol) -> Chem.Mol:
    """Align a molecule to a template using atom map IDs."""
    alignment = Alignment.from_mapids(source_mol, template_mol, source_map, template_map)
    alignment.apply()
    return source_mol


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


def _get_matched_mapids(mol: Chem.Mol, mol_map: Chem.Mol, is_template: bool = False) -> Dict[int, int]:
    """Get mapping between atom map IDs and atom indices for matched atoms.

    Args:
        mol: The molecule to match against
        mol_map: The molecule containing the atom mapping IDs to match
        is_template: Whether this is the template molecule (affects error message)

    Returns:
        Dictionary mapping atom map IDs to matched atom indices in mol

    Raises:
        ValueError: If mol_map is None, no atom maps are found, or no substructure match exists
    """
    # Get map ID to index mapping from the mapping molecule
    mapid2idx = _extract_map_ids(mol_map)
    if not mapid2idx:
        msg = "No atom maps found in template molecule" if is_template else "No atom maps found in source molecule"
        raise ValueError(msg)

    # Find substructure match between molecules
    match = mol.GetSubstructMatch(mol_map)
    if not match:
        msg = "No matching atom map IDs found"
        raise ValueError(msg)

    # Map the map IDs to the matched atom indices
    return {
        map_id: match[idx]
        for map_id, idx in mapid2idx.items()
        if idx < len(match) and match[idx] >= 0
    }
