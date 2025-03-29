from __future__ import annotations

import pickle
import re

import pandas as pd
import pytest
from rdkit import Chem

from xenopict import Xenopict


def test_pandas_df_style():
    df = pd.DataFrame.from_records([["O=C(O)C"], ["CCC"]], columns=["Smiles"])
    df["Xenopict"] = df["Smiles"].apply(Xenopict)  # type: ignore

    pandas_html = df.style.to_html()

    x1_html = df["Xenopict"].iloc[0].to_html()
    assert x1_html in pandas_html

    x2_html = df["Xenopict"].iloc[0].to_html()
    assert x2_html in pandas_html


def test_svg_ids_unique_and_stable():
    x = Xenopict("CCC")

    # all ids and hrefs are uniq
    x_id_href = _get_ids_and_hrefs(x.to_svg())
    assert len(x_id_href) == len(set(x_id_href))

    # same molecule produces same ids/hrefs
    y_id_href = _get_ids_and_hrefs(x.to_svg())
    assert x_id_href == y_id_href

    # when changes made to depiction, no ids/hrefs are shared
    x.shade([0.4, 0.3, 0.5])
    y1_id_href = _get_ids_and_hrefs(x.to_svg())
    assert len(set(x_id_href) & set(y1_id_href)) == 0

    # with different molecule, no ids/hrefs are shared
    z = Xenopict("CCC=O")
    z_id_href = _get_ids_and_hrefs(z.to_svg())
    assert len(set(z_id_href) & set(x_id_href)) == 0


def test_pickle_xenopict():
    x = Xenopict("CCC")

    # pickle and load
    pkx = pickle.dumps(x)
    pkx = pickle.loads(pkx)

    # same svg produced
    assert x.to_svg() == pkx.to_svg()


def test_smarts_sanitization_failure():
    m = Chem.MolFromSmarts("c1cc([NH2])ccc1")
    Xenopict(m)


def test_smarts_aromatic():
    # second bond should be dotted
    m = Chem.MolFromSmarts("c:c:c")
    assert "stroke-dasharray" in Xenopict(m).to_svg()


@pytest.mark.xfail
def test_clipping():
    m = Chem.MolFromSmarts("[CX4][Cl,Br,I]")
    Xenopict(m)
    assert False  # must manually check for now


def test_set_backbone_color():
    # Create a simple molecule
    x = Xenopict("CC(=O)O")  # Acetic acid

    # Set backbone color to red
    x.set_backbone_color("#FF0000")
    svg = x.to_svg()

    # Check that both lines and text have the correct color in the style attributes
    assert 'style="stroke:#FF0000' in svg  # Check bond color
    assert 'style="fill:#FF0000' in svg  # Check atom label color

    # Set backbone color to blue and verify change
    x.set_backbone_color("#0000FF")
    svg = x.to_svg()
    assert 'style="stroke:#0000FF' in svg
    assert 'style="fill:#0000FF' in svg


def test_filter_molecule():
    # Create a molecule with multiple atoms (benzene)
    x = Xenopict("c1ccccc1")

    # Filter to show only first three atoms and their connecting bonds
    atoms = [0, 1, 2]
    bonds = [(0, 1), (1, 2)]
    x.filter(atoms, bonds)
    svg = x.to_svg()

    # Check that filtered bonds are present (bonds contain atom references)
    for bond in bonds:
        assert f'class="bond-{bond[0]} atom-{bond[0]} atom-{bond[1]}"' in svg

    # Check that other bonds are not present
    assert "bond-3" not in svg
    assert "bond-4" not in svg
    assert "bond-5" not in svg


def test_from_smarts():
    # Test valid SMARTS pattern
    smarts = "C[OH]"  # Methanol pattern
    x = Xenopict.from_smarts(smarts)
    assert x is not None

    # Test that the SVG is generated correctly
    svg = x.to_svg()
    assert 'class="bond-0' in svg

    # Test invalid SMARTS pattern
    with pytest.raises(ValueError, match="Invalid SMARTS pattern"):
        Xenopict.from_smarts("invalid[pattern")


def test_substructure_focus():
    # Create a molecule with a substructure (phenol)
    x = Xenopict("Oc1ccccc1")

    # Focus on the hydroxyl group (O and connected C)
    substructure = [0, 1]  # O and C indices
    x.substructure_focus(substructure)
    svg = x.to_svg()

    # Check that focused atoms have normal opacity
    assert "opacity:1" in svg

    # Check that non-focused atoms have reduced opacity
    assert "opacity:0.2" in svg

    # Test with empty substructure (should not modify the molecule)
    x2 = Xenopict("Oc1ccccc1")
    x2.substructure_focus([])
    svg2 = x2.to_svg()
    assert "opacity:0.2" not in svg2


def test_mark_substructure():
    # Create a molecule with a substructure (ethanol)
    x = Xenopict("CCO")

    # Mark the hydroxyl group (O and connected C)
    atoms = [1, 2]  # C-O group
    bonds = [(1, 2)]  # C-O bond
    x.mark_substructure(atoms, bonds)
    svg = x.to_svg()

    # Check that the mark layer was created
    assert 'class="mark"' in svg

    # Check that the bond was marked with correct atom references
    assert 'class="bond-1 atom-1 atom-2"' in svg

    # Test with empty atoms list (should return without changes)
    x2 = Xenopict("CCO")
    x2.mark_substructure([])
    svg2 = x2.to_svg()
    assert 'class="mark"' not in svg2

    # Test with implicit bonds (not specified)
    x3 = Xenopict("CCO")
    x3.mark_substructure([1, 2])  # Same atoms, no explicit bonds
    svg3 = x3.to_svg()
    assert 'class="bond-1 atom-1 atom-2"' in svg3


def test_svg_path_optimization():
    # Create a molecule with repeated patterns (naphthalene)
    x = Xenopict("c1ccc2ccccc2c1")

    # Enable SVG optimization
    x.optimize_svg = True
    svg = x.to_svg()

    # Check that defs section was created for reusable elements
    assert "<defs>" in svg

    # Create a simple molecule (methane) where optimization should not create defs
    x2 = Xenopict("C")
    x2.optimize_svg = True
    svg2 = x2.to_svg()

    # Check that defs section was not created for simple molecule
    assert "<defs>" not in svg2


def _get_ids_and_hrefs(svg: str) -> list[str]:
    return list(re.findall(r'href=".+?"|id=".+?"', svg))
