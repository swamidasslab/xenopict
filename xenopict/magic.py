"""IPython magic commands for xenopict."""

import sys
import xml.dom.minidom

from rdkit.Chem import rdchem

from xenopict.monkey import BoostModulePatcher

from .drawer import Xenopict


def get_ipython():
    if "IPython" in sys.modules:
        return sys.modules["IPython"].core.getipython.get_ipython()  # type: ignore
    return None


def get_pandas():
    if get_ipython() or "pandas" in sys.modules:  # import pandas if IPython is available
        try:
            import pandas

            return pandas
        except ImportError:
            pass


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
    if ip := get_ipython():
        formatter = ip.display_formatter.formatters[  # type: ignore
            "image/svg+xml"
        ]  # ignore
        formatter.for_type(xml.dom.minidom.Document, _minidom_repr_svg)


def _minidom_repr_svg(doc):
    if doc.firstChild.tagName == "svg":
        return doc.toxml()
    raise NotImplementedError


#
# Patch/register RdKit
#


def _rdkit_repr_html(mol):
    """Formatter that uses Xenopict for rdchem.Mols"""
    return Xenopict(mol).to_html() if isinstance(mol, rdchem.Mol) else mol


def _rdkit_repr_svg(mol):
    """Formatter that uses Xenopict for rdchem.Mols"""
    return Xenopict(mol)._repr_svg_() if isinstance(mol, rdchem.Mol) else mol


def register_rdkit():
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


def _list_mol_html(input):
    if not input or len(input) > 50:
        return repr(input)

    if not isinstance(input[0], Xenopict):
        raise NotImplementedError

    divs = [f"<div style='border:solid 1px black;'>{item.to_html()}</div>" for item in input]

    return f'<div style="display:flex;flex-wrap:wrap;align-items:flex-start">{"".join(divs)}</div>'


def register_list_mol():
    if ip := get_ipython():
        formatter = ip.display_formatter.formatters[  # type: ignore
            "text/html"
        ]
        formatter.for_type(tuple, _list_mol_html)
        formatter.for_type(list, _list_mol_html)


#
# Patch Pandas


def _pandas_mol_repr_html(style):
    style.format(_rdkit_repr_html)
    return style


def patch_pandas():
    if pandas := get_pandas():
        if not hasattr(pandas.core.frame.DataFrame, "_xenopict"):
            pandas.core.frame.DataFrame._repr_html_orig_ = (  # type: ignore
                pandas.core.frame.DataFrame._repr_html_  # type: ignore
            )
            pandas.core.frame.DataFrame._repr_html_ = _pandas_mol_repr_html  # type: ignore


patcher = install()
uninstall = patcher.uninstall
