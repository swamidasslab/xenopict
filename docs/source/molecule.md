---
jupyter:
  jupytext:
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.1
  kernelspec:
    language: python
    name: python3
---
# Example

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]

from xenopict import Xenopict, magic
from rdkit import Chem
import rdkit.Chem.rdPartialCharges

import numpy as np


diclofenac = mol = Xenopict('O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl')
rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol.mol)
shading = np.array([a.GetDoubleProp("_GasteigerCharge")  for a in mol.GetAtoms()])
shading = shading / abs(shading).max()  # partial charge (scaled to [-1, 1])

# increase the default scale
Xenopict.scale = 30
```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# Atom shading 
drawer = Xenopict(mol, cmap="xenosite_bwr")
drawer.shade(shading).halo()
```

```python

```

```python
Xenopict(mol).shade(shading).halo()


```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
a1 = [b.GetBeginAtomIdx() for b in mol.GetBonds()]
a2 = [b.GetEndAtomIdx() for b in mol.GetBonds()]
bshading = (shading[a1] + shading[a2])  / 2

# Bond shading
drawer = Xenopict(mol)
drawer.shade(bond_shading=(a1, a2, bshading)).halo()
```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# Atom and bond shading togetehr
drawer = Xenopict(mol)
drawer.shade(shading, bond_shading=(a1, a2, bshading))
drawer.halo()
```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# Mark substructures (defined as a list of atom IDs)
drawer.mark_substructure([0,1,2,3, 10])
```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# Shade substructures
drawer = Xenopict(mol)
drawer.shade_substructure([[0,1,2,3], [7,8,9,10], [3,4,5,9]], [1, -0.65, 0.6])
drawer.halo()

```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# Another one of the xenosite colormaps
Xenopict(mol, cmap="xenosite_pwo").shade(shading).halo()
```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# any matplotlib colormap works, but most won't look as good as the default
Xenopict(mol, cmap="RdBu", scale=30).shade(shading).halo()
```

```python papermill={"duration": null, "end_time": null, "exception": null, "start_time": null, "status": "completed"} tags=[]
# depeict a substructure
drawer = Xenopict(mol)
drawer.substructure_focus([0,1,2,3])
```

```python
Xenopict('O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl').mol
```

```python
mol.mol.AdjustQueryProperties
```

```python

```
