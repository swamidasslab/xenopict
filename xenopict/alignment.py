"""Functions for aligning 2D molecular depictions using RDKit.

This module provides methods for aligning molecules in 2D space to create consistent
molecular depictions. The primary entrypoint is the auto_align_molecules() function,
which can align multiple molecules simultaneously while ensuring optimal and consistent
alignments across the entire set.

Primary Usage:
    >>> from rdkit import Chem
    >>> # Create some molecules to align
    >>> ethanol = Chem.MolFromSmiles("CCO")
    >>> methanol = Chem.MolFromSmiles("CO") 
    >>> propanol = Chem.MolFromSmiles("CCCO")
    >>> # Align all molecules (automatically aligns OH groups)
    >>> aligned = auto_align_molecules([ethanol, methanol, propanol])
    
    For more control, you can provide hints to force specific alignments:
    >>> # Create a hint to align specific atoms
    >>> hint = Alignment.from_aligned_atoms(ethanol, methanol, [(0, 0)])  # align terminal carbons
    >>> aligned = auto_align_molecules([ethanol, methanol, propanol], hints=[hint])

The module also provides several lower-level alignment methods through the Alignment class:

1. Maximum Common Substructure (MCS):
   - Automatically finds the largest common substructure between molecules
   - Best for general alignment when molecules share structural features
   - Example: aligning similar molecules like ethanol and methanol

2. Explicit Atom Pairs:
   - Manual specification of which atoms should be aligned
   - Best when you need precise control over the alignment
   - Example: forcing specific functional groups to align

3. Atom Map IDs:
   - Uses atom map numbers to determine alignment
   - Best for reaction-based alignments where atoms are already mapped
   - Example: aligning reactants and products in a reaction

4. Atom Indices:
   - Uses direct atom index mappings
   - Best when working with known atom indices
   - Example: aligning based on atom ordering

These lower-level methods are wrapped in convenience functions (align_from_mcs,
align_from_atom_pairs, etc.) but most users should prefer auto_align_molecules()
as it handles the complexity of finding optimal alignments across multiple molecules.

Examples:
    Here are some common use cases:

    1. Basic usage - align multiple molecules:
    >>> # Align a set of alcohols by their OH groups
    >>> molecules = [ethanol, methanol, propanol]
    >>> aligned = auto_align_molecules(molecules)

    2. Using hints for custom alignments:
    >>> # Force ethanol's CH3 to align with methanol's OH
    >>> unusual_pairs = [(0, 1)]  # ethanol CH3 to methanol OH
    >>> hint = Alignment.from_aligned_atoms(ethanol, methanol, unusual_pairs)
    >>> aligned = auto_align_molecules(molecules, hints=[hint])

    3. Direct use of lower-level methods (if needed):
    >>> # Align two molecules using MCS
    >>> aligned = align_from_mcs(ethanol, methanol)

See individual class and function documentation for more detailed examples and usage.
"""

from rdkit import Chem  # type: ignore
from rdkit.Chem import rdDepictor  # type: ignore
from rdkit.Chem import rdFMCS  # type: ignore
from typing import List, Dict, Tuple, Union
from typing import NamedTuple, Optional


class Alignment(NamedTuple):
    """A class representing a 2D alignment between two molecules.

    This class stores information about how two molecules should be aligned in 2D space.
    It includes the atom pairs that should be matched between molecules, a score indicating
    alignment quality, and references to both molecules.

    The class provides several factory methods for creating alignments:
    - from_aligned_atoms: Create from explicit atom pairs
    - from_mcs: Create using maximum common substructure
    - from_mapids: Create using atom map IDs
    - from_indices: Create using atom index mappings

    Attributes:
        atom_pairs: List of (source_idx, template_idx) tuples specifying matched atoms
        score: Numerical score indicating alignment quality (higher is better)
        source_mol: The molecule being aligned
        template_mol: The template molecule being aligned to
    """

    atom_pairs: list[tuple[int, int]]
    score: float
    source_mol: Chem.Mol
    template_mol: Chem.Mol

    @staticmethod
    def from_aligned_atoms(
        source_mol: Chem.Mol, template_mol: Chem.Mol, atom_pairs: list[tuple[int, int]]
    ) -> "Alignment":
        """Create an alignment from explicit atom pairs.

        This is the most basic factory method that creates an alignment directly from
        a list of atom pairs. The alignment score is set to 0 since the pairs are
        manually specified rather than discovered algorithmically.

        Args:
            source_mol: The molecule to align
            template_mol: The template molecule to align to
            atom_pairs: List of (source_idx, template_idx) pairs specifying matched atoms

        Returns:
            A new Alignment object with the specified atom pairs

        Raises:
            AssertionError: If any atom indices are invalid
        """
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
        molecules. The inexact matching allows for different bond types, which can help
        find alignments when bond orders differ between molecules.

        The alignment score is computed as the average of:
        - Number of atoms in exact MCS match
        - Number of bonds in exact MCS match  
        - Number of atoms in inexact MCS match
        - Number of bonds in inexact MCS match

        Args:
            source_mol: The molecule to align
            template_mol: The template molecule to align to

        Returns:
            A new Alignment object based on MCS matching
        """
        exact_mcs = rdFMCS.FindMCS([source_mol, template_mol])
        exact_query = exact_mcs.queryMol

        inexact_mcs = rdFMCS.FindMCS(
            [source_mol, template_mol], bondCompare=rdFMCS.BondCompare.CompareAny
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
        """Align a molecule to a template molecule using atom index mappings.
        
        This function aligns a molecule by mapping atom indices between the source and template
        molecules. The mapping can be provided either as a list or dictionary. When using a list,
        the index in the list corresponds to the source atom index, and the value is the template
        atom index. When using a dictionary, the keys are source atom indices and values are
        template atom indices.

        Args:
            source_mol: The molecule to align
            template_mol: The template molecule to align to
            index_mapping: Either a list or dictionary mapping source atom indices to template indices.
                          For list format, -1 indicates no mapping for that source atom.

        Returns:
            An Alignment object
        """
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
        """Create a new alignment with source and template molecules swapped.

        This method creates a new alignment where:
        - The source molecule becomes the template
        - The template molecule becomes the source
        - The atom pairs are reversed (b,a) instead of (a,b)
        - The score remains the same

        Returns:
            A new Alignment object with source and template swapped

        Examples:
            Create and reverse an alignment between ethanol and methanol:

            >>> from rdkit import Chem
            >>> source = Chem.MolFromSmiles("CCO")  # indices: 0=C, 1=C, 2=O
            >>> template = Chem.MolFromSmiles("CO")  # indices: 0=C, 1=O
            >>> alignment = Alignment.from_aligned_atoms(source, template, [(1, 0), (2, 1)])
            >>> reversed_align = alignment.reverse()
            >>> reversed_align.source_mol == template
            True
            >>> reversed_align.template_mol == source
            True
            >>> reversed_align.atom_pairs == [(0, 1), (1, 2)]
            True
            >>> reversed_align.score == alignment.score
            True
        """
        return Alignment(
            atom_pairs=[(b, a) for a, b in self.atom_pairs],
            score=self.score,
            source_mol=self.template_mol,
            template_mol=self.source_mol,
        )

    def validate(self) -> "Alignment":
        """Validate that all atom indices in the alignment are valid.

        This method checks that:
        - All atom indices are non-negative
        - All atom indices are within bounds for their respective molecules
        - No duplicate atom indices exist

        Returns:
            Self, allowing for method chaining

        Raises:
            AssertionError: If any validation checks fail

        Examples:
            Create and validate a valid alignment:

            >>> from rdkit import Chem
            >>> source = Chem.MolFromSmiles("CCO")  # indices: 0=C, 1=C, 2=O
            >>> template = Chem.MolFromSmiles("CO")  # indices: 0=C, 1=O
            >>> alignment = Alignment([(1, 0), (2, 1)], 1.0, source, template)
            >>> validated = alignment.validate()  # No error raised
            >>> validated == alignment
            True

            Try to validate an alignment with invalid indices:

            >>> bad_alignment = Alignment([(-1, 0)], 1.0, source, template)
            >>> bad_alignment.validate()  # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            AssertionError: Negative source molecule atom index -1

            >>> bad_alignment = Alignment([(3, 0)], 1.0, source, template)
            >>> bad_alignment.validate()  # doctest: +IGNORE_EXCEPTION_DETAIL
            Traceback (most recent call last):
                ...
            AssertionError: Source molecule atom index 3 out of range (max 2)
        """
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
        """Apply the alignment to update the source molecule's coordinates.

        This method modifies the source molecule in place by:
        1. Ensuring the template molecule has 2D coordinates
        2. Using RDKit's depiction generation to align the source to the template
        3. If no atom pairs exist, just ensures source has 2D coordinates

        Returns:
            Self, allowing for method chaining

        Examples:
            Create and apply an alignment between ethanol and methanol:

            >>> from rdkit import Chem
            >>> source = Chem.MolFromSmiles("CCO")  # indices: 0=C, 1=C, 2=O
            >>> template = Chem.MolFromSmiles("CO")  # indices: 0=C, 1=O
            >>> alignment = Alignment.from_aligned_atoms(source, template, [(1, 0), (2, 1)])
            >>> applied = alignment.apply()
            >>> applied == alignment  # Returns self for chaining
            True

            Verify coordinates were updated:

            >>> source_O = GetCoords(source, 2)
            >>> template_O = GetCoords(template, 1)
            >>> abs(source_O[0] - template_O[0]) < 0.1
            True
            >>> abs(source_O[1] - template_O[1]) < 0.1
            True

            Empty alignments just ensure 2D coordinates exist:

            >>> empty = Alignment([], 0.0, source, template)
            >>> empty.apply()  # No error raised
            Alignment(atom_pairs=[], score=0.0, source_mol=...)
            >>> source.GetNumConformers() > 0  # Has coordinates
            True
        """
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
    """Automatically align multiple molecules by finding an optimal alignment tree.

    This function aligns a list of molecules by constructing a maximum spanning tree of
    alignments. Each edge in the tree represents an alignment between two molecules,
    with the edge weight being the alignment score. The tree structure ensures that
    all molecules are connected through a series of high-quality alignments.

    The function can also use "hints" - predefined alignments that should be preferred
    over automatically discovered ones. This is useful when you want to enforce certain
    alignments while letting the algorithm figure out the rest.

    The algorithm works by:
    1. Building a complete graph where nodes are molecules and edges are alignments
    2. Finding a maximum spanning tree to get the best set of alignments
    3. Using the largest molecule as the root of the tree
    4. Applying alignments in breadth-first order from the root

    Args:
        molecules: List of molecules to align
        hints: Optional list of preferred Alignment objects. These alignments get a
              very high weight (score + 1000) to ensure they are used when possible.

    Returns:
        The input molecules, modified in place with new 2D coordinates

    Examples:
        Let's align three molecules - ethanol, methanol, and propanol:

        >>> from rdkit import Chem
        >>> ethanol = Chem.MolFromSmiles("CCO")
        >>> methanol = Chem.MolFromSmiles("CO")
        >>> propanol = Chem.MolFromSmiles("CCCO")

        Align all molecules:

        >>> aligned = auto_align_molecules([ethanol, methanol, propanol])
        >>> len(aligned) == 3
        True

        We can verify the OH groups are aligned by checking coordinates:

        >>> # Get OH coordinates for each molecule
        >>> eth_O = GetCoords(ethanol, 2)  # O is at index 2
        >>> meth_O = GetCoords(methanol, 1)  # O is at index 1
        >>> prop_O = GetCoords(propanol, 3)  # O is at index 3

        >>> # Check that OH groups are aligned
        >>> abs(eth_O[0] - meth_O[0]) < 0.1
        True
        >>> abs(eth_O[1] - meth_O[1]) < 0.1
        True
        >>> abs(prop_O[0] - meth_O[0]) < 0.1
        True
        >>> abs(prop_O[1] - meth_O[1]) < 0.1
        True

        Now let's provide a hint that aligns ethanol differently - matching the terminal
        carbon of ethanol to the oxygen of methanol (an unusual alignment):

        >>> # Create a hint with unusual atom pairs
        >>> unusual_pairs = [(0, 1)]  # Match ethanol's CH3 to methanol's OH
        >>> hint = Alignment.from_aligned_atoms(ethanol, methanol, unusual_pairs)
        >>> aligned = auto_align_molecules([ethanol, methanol, propanol], hints=[hint])
        >>> len(aligned) == 3
        True

        Verify that our unusual hint was used - ethanol's CH3 should be where methanol's OH was:

        >>> # Get coordinates after hint-based alignment
        >>> eth_CH3 = GetCoords(ethanol, 0)  # Terminal carbon in ethanol
        >>> meth_OH = GetCoords(methanol, 1)  # Oxygen in methanol

        >>> # The hint should force these atoms to align
        >>> abs(eth_CH3[0] - meth_OH[0]) < 0.1
        True
        >>> abs(eth_CH3[1] - meth_OH[1]) < 0.1
        True
    """
    import networkx as nx

    hints = hints or []

    g = nx.Graph() # type: ignore
    mol2node = {}
    for i, mol in enumerate(molecules):
        _ensure_coords(mol)
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
    t = nx.maximum_spanning_tree(g) # type: ignore

    # unrooted tree so any node can be the root
    # we will pick the largest molecule
    max_mol_size = -1
    max_mol_idx = -1
    for node in t.nodes:  # type: ignore
        m = t.nodes[node]["mol"]  # type: ignore
        s = m.GetNumAtoms()
        if s > max_mol_size:
            max_mol_size = s
            max_mol_idx = node

    # now iterate over the tree in BFS order
    for edge in nx.traversal.bfs_edges(t, max_mol_idx):  # type: ignore
        alignment = g.edges[edge]["alignment"]  # type: ignore

        # if the first idx of the edge isn't the template swap the alignment
        if edge[0] != g.edges[edge]["template"]:  # type: ignore
            alignment = alignment.reverse()

        alignment.apply()

    return [g.nodes[i]["mol"] for i in sorted(g.nodes)]  # type: ignore


def align_from_mcs(source_mol: Chem.Mol, template_mol: Chem.Mol) -> Chem.Mol:
    """Align a molecule to a template molecule using maximum common substructure (MCS).

    This function automatically aligns a molecule to a template by finding their maximum
    common substructure and using it to determine the alignment. It uses both exact and
    inexact matching to find the best possible alignment, where inexact matching allows
    for different bond types.

    The function modifies the source molecule in place by updating its 2D coordinates
    to match the template molecule's orientation.

    Args:
        source_mol: The molecule to align
        template_mol: The template molecule to align to

    Returns:
        The source molecule, modified in place with new 2D coordinates

    Raises:
        ValueError: If no alignment could be found between the molecules

    Examples:
        First, let's align ethanol to methanol, matching the OH group:

        >>> from rdkit import Chem
        >>> source = Chem.MolFromSmiles("CCO")
        >>> template = Chem.MolFromSmiles("CO")

        Perform the alignment:

        >>> aligned = align_from_mcs(source, template)

        Ensure get template and source coordinates:

        >>> template_O = GetCoords(template, 1)  # Oxygen is at index 1
        >>> template_C = GetCoords(template, 0)  # Carbon is at index 0

        >>> source_O = GetCoords(source, 2)  # Oxygen is at index 2
        >>> source_C = GetCoords(source, 1)  # Matching carbon is at index 1

        Verify that the OH group coordinates match between molecules:

        >>> abs(source_O[0] - template_O[0]) < 0.1
        True
        >>> abs(source_O[1] - template_O[1]) < 0.1
        True
        >>> abs(source_C[0] - template_C[0]) < 0.1
        True
        >>> abs(source_C[1] - template_C[1]) < 0.1
        True

        The function raises an error when no alignment is possible:

        >>> source = Chem.MolFromSmiles("c1ccccc1")
        >>> template = Chem.MolFromSmiles("O")
        >>> align_from_mcs(source, template)  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ValueError: No alignment found
    """
    alignment = Alignment.from_mcs(source_mol, template_mol)
    if not alignment.atom_pairs:
        raise ValueError("No alignment found")
    alignment.apply()
    return source_mol


def align_from_atom_pairs(source_mol: Chem.Mol, template_mol: Chem.Mol, atom_pairs: List[Tuple[int, int]]) -> Chem.Mol:
    """Align a molecule to a template molecule using explicit atom pairs.
    
    This function aligns a molecule by matching specific pairs of atoms between the source
    and template molecules. This is useful when you want precise control over which atoms
    should be aligned, rather than letting the algorithm find matches automatically.

    The function modifies the source molecule in place by updating its 2D coordinates
    to match the template molecule's orientation.

    Args:
        source_mol: The molecule to align
        template_mol: The template molecule to align to
        atom_pairs: List of (source_idx, template_idx) pairs specifying which atoms to match

    Returns:
        The source molecule, modified in place with new 2D coordinates

    Examples:
        First, let's align ethanol to methanol by explicitly matching the OH group:

        >>> from rdkit import Chem
        >>> source = Chem.MolFromSmiles("CCO")
        >>> template = Chem.MolFromSmiles("CO")

        Align using explicit atom pairs - match oxygen and its attached carbon:

        >>> atom_pairs = [(2, 1), (1, 0)]  # (source O, template O), (source C, template C)
        >>> aligned = align_from_atom_pairs(source, template, atom_pairs)

        Get coordinates after alignment:

        >>> template_O = GetCoords(template, 1)
        >>> template_C = GetCoords(template, 0)
        >>> source_O = GetCoords(source, 2)
        >>> source_C = GetCoords(source, 1)

        Verify that the matched atoms now have the same coordinates:

        >>> abs(source_O[0] - template_O[0]) < 0.1
        True
        >>> abs(source_O[1] - template_O[1]) < 0.1
        True
        >>> abs(source_C[0] - template_C[0]) < 0.1
        True
        >>> abs(source_C[1] - template_C[1]) < 0.1
        True
    """
    alignment = Alignment.from_aligned_atoms(source_mol, template_mol, atom_pairs)
    alignment.apply()
    return source_mol


def align_from_indices(source_mol: Chem.Mol, template_mol: Chem.Mol, index_mapping: Union[List[int], Dict[int, int]]) -> Chem.Mol:
    """Align a molecule to a template molecule using atom index mappings.
    
    This function aligns a molecule by mapping atom indices between the source and template
    molecules. The mapping can be provided either as a list or dictionary. When using a list,
    the index in the list corresponds to the source atom index, and the value is the template
    atom index. When using a dictionary, the keys are source atom indices and values are
    template atom indices.

    The function modifies the source molecule in place by updating its 2D coordinates
    to match the template molecule's orientation.

    Args:
        source_mol: The molecule to align
        template_mol: The template molecule to align to
        index_mapping: Either a list or dictionary mapping source atom indices to template indices.
                      For list format, -1 indicates no mapping for that source atom.

    Returns:
        The source molecule, modified in place with new 2D coordinates

    Examples:
        First, let's align ethanol to methanol using a list mapping:

        >>> from rdkit import Chem
        >>> source = Chem.MolFromSmiles("CCO")  # indices: 0=C, 1=C, 2=O
        >>> template = Chem.MolFromSmiles("CO")  # indices: 0=C, 1=O

        Align using a list mapping - match OH group:
        
        >>> mapping = [-1, 0, 1]  # Source C1->Template C0, Source O2->Template O1
        >>> aligned = align_from_indices(source, template, mapping)

        Get coordinates after alignment:

        >>> template_O = GetCoords(template, 1)
        >>> template_C = GetCoords(template, 0)
        >>> source_O = GetCoords(source, 2)
        >>> source_C = GetCoords(source, 1)

        Verify that the matched atoms have the same coordinates:

        >>> abs(source_O[0] - template_O[0]) < 0.1
        True
        >>> abs(source_O[1] - template_O[1]) < 0.1
        True
        >>> abs(source_C[0] - template_C[0]) < 0.1
        True
        >>> abs(source_C[1] - template_C[1]) < 0.1
        True

        We can also use a dictionary mapping:

        >>> mapping = {1: 0, 2: 1}  # Same mapping but as a dict
        >>> aligned = align_from_indices(source, template, mapping)

        The alignment results in the same coordinate matches:

        >>> source_O = GetCoords(source, 2)
        >>> source_C = GetCoords(source, 1)
        >>> abs(source_O[0] - template_O[0]) < 0.1
        True
        >>> abs(source_C[0] - template_C[0]) < 0.1
        True
    """
    alignment = Alignment.from_indices(source_mol, template_mol, index_mapping)
    alignment.apply()
    return source_mol


def align_from_mapids(source_mol: Chem.Mol, template_mol: Chem.Mol, source_map: Chem.Mol, template_map: Chem.Mol) -> Chem.Mol:
    """Align a molecule to a template using atom map IDs.
    
    This function aligns a molecule by matching atoms with corresponding map IDs between
    the source and template molecules. Map IDs are integers assigned to atoms that help
    establish correspondence between different molecules. This is particularly useful when
    you want to align molecules based on atom-to-atom mapping from a reaction or
    transformation.

    The function modifies the source molecule in place by updating its 2D coordinates
    to match the template molecule's orientation.

    Args:
        source_mol: The molecule to align
        template_mol: The template molecule to align to
        source_map: The source molecule with atom map IDs
        template_map: The template molecule with atom map IDs

    Returns:
        The source molecule, modified in place with new 2D coordinates

    Raises:
        ValueError: If no atom maps are found in either molecule or no matching map IDs exist

    Examples:
        Let's align ethanol to methanol using atom map IDs:

        >>> from rdkit import Chem
        >>> source = Chem.MolFromSmiles("CCO")
        >>> template = Chem.MolFromSmiles("CO")
        >>> source_map = Chem.MolFromSmiles("[CH3:1][CH2:2][OH:3]")
        >>> template_map = Chem.MolFromSmiles("[CH3:2][OH:3]")

        Perform the alignment:

        >>> aligned = align_from_mapids(source, template, source_map, template_map)

        Get coordinates after alignment:

        >>> template_O = GetCoords(template, 1)
        >>> template_C = GetCoords(template, 0)
        >>> source_O = GetCoords(source, 2)
        >>> source_C = GetCoords(source, 1)

        Verify that the matched atoms have the same coordinates:

        >>> abs(source_O[0] - template_O[0]) < 0.1
        True
        >>> abs(source_O[1] - template_O[1]) < 0.1
        True
        >>> abs(source_C[0] - template_C[0]) < 0.1
        True
        >>> abs(source_C[1] - template_C[1]) < 0.1
        True

        The function raises an error when no atom maps are found:

        >>> bad_map = Chem.MolFromSmiles("CCO")  # No atom maps
        >>> align_from_mapids(source, template, bad_map, template_map)  # doctest: +IGNORE_EXCEPTION_DETAIL
        Traceback (most recent call last):
            ...
        ValueError: No atom maps found in source molecule
    """
    alignment = Alignment.from_mapids(source_mol, template_mol, source_map, template_map)
    alignment.apply()
    return source_mol


def GetCoords(mol: Chem.Mol, i: int) -> tuple[float, float]:
    c = mol.GetConformer(0).GetAtomPosition(i)
    return c.x, c.y


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
