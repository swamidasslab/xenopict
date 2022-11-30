from xenopict import Xenopict, magic
import pandas as pd
import pickle
import re
from rdkit import Chem
import pytest


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
    # second bond is not dotted
    m = Chem.MolFromSmarts("[CX4][Cl,Br,I]")
    Xenopict(m)
    assert False  # must manually check for now


def _get_ids_and_hrefs(svg: str) -> list[str]:
    return list(re.findall(r'href=".+?"|id=".+?"', svg))
