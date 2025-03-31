---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: 'Python 3.9.13 (''python-3.9.13'': venv)'
  language: python
  name: python3
---

# Quickstart

Xenopict is a python library for depicting molecules.

## Installation

Xenopict installs from pypi.

```console
pip install xenopict
```

The package is only tested on Python 3.10 and 3.9. It may run on earlier versions of python, but we do not intend to maintain compatibility.

## Reading in a Molecule

The main entrypoint to xenopict is the Xenopict class, which is initialized with a small molecule. The class can take a SMILES String

```{code-cell} ipython3
from xenopict import Xenopict     
mol = Xenopict('O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl')
```

It can also be initialized with a RDKit mol.

```{code-cell} ipython3
from rdkit import Chem
mol = Chem.MolFromSmiles('O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl')
diclofenac = mol = Xenopict(mol)
```

Now that a molecule is loaded, we can display it easily in a notebook.

```{code-cell} ipython3
from xenopict import magic # the magic package enables several integrations with notebooks
mol
```

## Atom and Bond Shading

Xenopict can shade molecules by with an intensity, which can range from 0 to 1 for sequential colormaps, and -1 to 1 for diverging colormaps. This shading is  implemented in concise vector graphics. Xenosite produces an SVG, which can be accessed with the 'to_svg' method.

To demonstrate, we will shade a molecule by its absolute partial charge [0 to 1], which we compute from the RDKit library. Each atom is assigned a number in a vector.

```{code-cell} ipython3
import rdkit.Chem.rdPartialCharges
import numpy as np
rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol.mol)
shading = np.array([a.GetDoubleProp("_GasteigerCharge")  for a in mol.GetAtoms()])
shading = abs(shading / abs(shading).max() ) # partial charge (scaled to [0, 1])
shading
```

The shading is applied to molecule using the "shade" method. The "halo" method adds a halo effect that aids visibility. 

```{code-cell} ipython3
Xenopict(mol).shade(shading).halo()
```

We can shade molecules by bond too. To demonstrate, we will create a shading vector for bonds that is the average absolute charge of the atom at each end. 

```{code-cell} ipython3
a1 = [b.GetBeginAtomIdx() for b in mol.GetBonds()]
a2 = [b.GetEndAtomIdx() for b in mol.GetBonds()]
bshading = (shading[a1] + shading[a2])  / 2
bshading
```

```{code-cell} ipython3
Xenopict(mol).shade(bond_shading=(a1, a2, bshading)).halo()
```

Atoms and bonds can be shaded together at the same time.

```{code-cell} ipython3
x = Xenopict(mol).shade(shading, bond_shading=(a1, a2, bshading)).halo()
x
```

## Marking, Shading and Focusing Substructures

At times we want to mark a substructure. We can do so by supplying a list of atom idxs.

```{code-cell} ipython3
x.mark_substructure([0,1,2,3])
```

We can also shade substructures, by specifying them as lists of IDs, wiht a list of shading values. Substructures with intensities closer to zero are painted behind those with higher intensities. 

```{code-cell} ipython3
substr1 = [0,1,2,3]
substr2 = [7,8,9,10]
substr3 = [3,4,5,9]

shading = [1.0, 0.5, 0.6]

Xenopict(mol).shade_substructure([substr1, substr2, substr3], shading)
```

By default, xenopict adds in all the bonds joining a group of atoms. But sometimes this includes
bonds that were not wanted. Substructures can be more precisely specified by including a list of bonds. 

```{code-cell} ipython3
mol = Chem.MolFromSmiles("c1ccccc1CCC")

atoms = [0,1,2,3,4,5]
bonds = [[0,1], [2,1], [2,3], [3,4],[4, 5], [5,6]]

Xenopict(mol).shade_substructure([atoms], [1.0]), \
Xenopict(mol).shade_substructure([atoms], [1.0], [bonds])
```

Alongside shadings or alternately to them, we can also depict a substructure in isolation.

```{code-cell} ipython3
Xenopict(mol).substructure_focus(atoms), \
Xenopict(mol).substructure_focus(atoms, bonds)
```

## Adding Atom Indices

Molecules can be drawn with indices too.

```{code-cell} ipython3
x = Xenopict(mol,scale=30,  add_atom_indices=True)
x
```

Atom IDs are carried over to substructure depictions, but they are not filtered out yet.

```{code-cell} ipython3
x.substructure_focus([4,5,6])
```
