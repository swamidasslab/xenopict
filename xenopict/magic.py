from rdkit.Chem import rdchem
from xenopict import Xenopict
from xenopict.monkey import BoostModulePatcher

import xml.dom.minidom
from pml import HTML


def install():
    register_rdkit()
    register_minidom()
    register_list_mol()

    patcher = patch_rdkit()
    patcher.install()

    return patcher


#
# Register minidom formatter so SVG doms display as images.
#


def register_minidom():
    from IPython import get_ipython

    if ip := get_ipython():
        formatter = ip.display_formatter.formatters[  # type: ignore
            "image/svg+xml"
        ]  # ignore
        formatter.for_type(xml.dom.minidom.Document, _minidom_repr_svg)


def _minidom_repr_svg(doc):
    if doc.firstChild.tagName == "svg":
        return doc.toxml()


#
# Patch/register RdKit
#


def _rdkit_repr_html(mol):
    """Formatter that uses Xenopict for rdchem.Mols"""
    return Xenopict(mol).to_html() if isinstance(mol, rdchem.Mol) else mol


def _rdkit_repr_svg(mol):
    """Formatter that uses Xenopict for rdchem.Mols"""
    return Xenopict(mol).to_svg() if isinstance(mol, rdchem.Mol) else mol


def register_rdkit():
    from IPython import get_ipython

    if ip := get_ipython():
        formatter = ip.display_formatter.formatters[  # type: ignore
            "image/svg+xml"
        ]  # ignore
        formatter.for_type(rdchem.Mol, _rdkit_repr_svg)


def patch_rdkit():
    patcher = BoostModulePatcher()
    patcher.replace(rdchem.Mol, "__str__", _rdkit_repr_html)
    patcher.replace(rdchem.Mol, "_repr_svg_", _rdkit_repr_svg)
    # rdchem.Mol.__str__ = _rdkit_repr_html
    # rdchem.Mol._repr_svg_ = _rdkit_repr_svg

    return patcher


#
# Register list and tuple
#


def _list_mol_html(input):
    if not input or len(input) > 50:
        return repr(input)

    if hasattr(input[0], "_repr_svg_"):
        h = HTML().div(style="display:flex;flex-wrap:wrap;align-items:flex-start")  # type: ignore

        for item in input:
            r = item._repr_svg_() if hasattr(item, "_repr_svg_") else repr(item)
            h.div(style="margin:0.1em;background:white").div(r, escape=False)  # type: ignore
        return str(h)

    return repr(input)


def register_list_mol():
    from IPython import get_ipython

    if ip := get_ipython():
        formatter = get_ipython().display_formatter.formatters[  # type: ignore
            "text/html"
        ]  # ignore
        formatter.for_type(tuple, _list_mol_html)
        formatter.for_type(list, _list_mol_html)


#
# Patch Pandas
#


def _pandas_mol_repr_html(style):
    style.format(_rdkit_repr_html)
    return style


def patch_pandas():
    import pandas.core.frame

    if not hasattr(pandas.core.frame.DataFrame, "_xenopict"):
        pandas.core.frame.DataFrame._repr_html_orig_ = (  # type: ignore
            pandas.core.frame.DataFrame._repr_html_  # type: ignore
        )
        pandas.core.frame.DataFrame._repr_html_ = _pandas_mol_repr_html  # type: ignore


patcher = install()
uninstall = patcher.uninstall
