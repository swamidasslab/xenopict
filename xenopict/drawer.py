from __future__ import annotations
import numpy as np
from xml.dom.minidom import parseString, Element
import contextlib
from six.moves.collections_abc import Sequence, Mapping  # type: ignore
from rdkit.Chem.Draw import rdMolDraw2D, rdDepictor
from rdkit.Chem.rdchem import Mol
from rdkit.Chem import MolFromSmiles, MolFromSmarts  # type: ignore
from .colormap import install_colormaps
from .plotdot import PlotDot
from .alignment import align_from_mcs

from urllib.parse import quote
from collections import defaultdict
from typing import Optional, Union, Tuple
import simplejson as json
import hashlib
import re
import os

from warnings import warn

with contextlib.suppress(ImportError):
    from matplotlib.colors import Colormap

from shapely.geometry import LineString, Point

install_colormaps()

__all__ = ["Xenopict"]

_DEBUG = os.environ.get("XENOPICT_DEBUG", False)

if _DEBUG:
    from icecream import ic
else:
    ic = lambda x: x


AtomIdx = int

import sys
MINOR_VERSION = int(sys.version.split(".")[1])


if MINOR_VERSION >= 9:
  BondShading = tuple[Sequence[AtomIdx], Sequence[AtomIdx], Sequence[float]]
  AtomShading = Sequence[float]
SVG = str


def _style2dict(s: str) -> dict[str, str]:
    if not s:
        return {}
    s = s.strip().strip(";")
    return dict([e.split(":") for e in s.split(";")]) if s else {}


def _dict2style(s: Mapping[str, str]) -> str:
    return ";".join(["%s:%s" % i for i in s.items()])


class Xenopict:
    """
    This class draws an RDK molecule with sensible defaults,
    cleaning up the output SVG for easier modification and
    reduced size.

    >>> from rdkit import Chem
    >>> import rdkit.Chem.rdPartialCharges
    >>> diclofenac = mol = Chem.MolFromSmiles('O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl')
    >>> rdkit.Chem.rdPartialCharges.ComputeGasteigerCharges(mol)

    Partial charge (scaled to be in [-1, 1])

    >>> shading = np.array([a.GetDoubleProp("_GasteigerCharge")  for a in mol.GetAtoms()])
    >>> shading = shading / abs(shading).max()

    SVG of molecule shaded by partial charge,

    >>> drawer = Xenopict(mol)
    >>> drawer.shade(shading)
    <xenopict.drawer.Xenopict ...>
    >>> str(drawer)
    '<...>'

    Atoms can also be annotated with a circle,

    >>> drawer.mark_atoms([1, 2])
    <xenopict.drawer.Xenopict ...>

    The underlying svg dom (xml.dom.minidom) is accessible:

    >>> drawer.svgdom
    <xml.dom.minidom.Document ...>
    """

    down_scale: float = 0.7
    mark_down_scale: float = 1.0
    shapely_resolution: int = 6
    scale: float = 20
    diverging_cmap: bool = False
    add_atom_indices: bool = False
    add_bond_indices: bool = False
    optimize_svg: bool = True
    embed_script: bool = False
    #kekulize : bool = False
    dummies_are_attachments : bool = False
    plot_dot: PlotDot = PlotDot()
    cmap: Union[str, "Colormap"] = "xenosite"

    def set_backbone_color(self, color: str) -> "Xenopict":
        """Set the color of the molecule's backbone and atom labels.
        
        Args:
            color: Hex color code (e.g. '#FF0000' for red)
            
        Returns:
            self for method chaining
        """
        # Set color for lines (bonds)
        lines_group = self.groups["lines"]
        if lines_group:
            style = lines_group.getAttribute("style")
            style_dict = _style2dict(style)
            style_dict["stroke"] = color
            lines_group.setAttribute("style", _dict2style(style_dict))
            
        # Set color for text (atom labels)
        text_group = self.groups["text"]
        if text_group:
            style = text_group.getAttribute("style") or ""
            style_dict = _style2dict(style)
            style_dict["fill"] = color
            style_dict["stroke"] = "none"  # Ensure no stroke on text
            text_group.setAttribute("style", _dict2style(style_dict))
            
        return self

    @classmethod
    def from_smarts(cls, smarts: str, **kwargs) -> "Xenopict":
        """Create a Xenopict object from a SMARTS pattern.
        
        Args:
            smarts: SMARTS pattern string
            **kwargs: Additional arguments passed to Xenopict constructor
            
        Returns:
            Xenopict object initialized with the SMARTS pattern
        """
        mol = MolFromSmarts(smarts)
        if mol is None:
            raise ValueError(f"Invalid SMARTS pattern: {smarts}")
        return cls(mol, **kwargs)

    def __init__(
        self, input_mol: Union[str, Mol, "Xenopict"], **kwargs
    ):  # sourcery skip: use-dict-items
        if isinstance(input_mol, Mol):
            mol: Mol = input_mol
        elif isinstance(input_mol, str):
            mol: Mol = MolFromSmiles(input_mol)
        elif isinstance(input_mol, Xenopict):
            mol: Mol = input_mol.mol
        else:
            raise ValueError("Input must be RDMol or smiles string.")

        self.__dict__.update(kwargs)

        self.mol: Mol = mol

        self.groups: dict[str, Element] = {}
        self.draw_mol()

    def draw_mol(self, mol: Optional[Mol] = None):
        self.mol: Mol = mol or self.mol

        self._filter = None
        d2d = rdMolDraw2D.MolDraw2DSVG(-1, -1)

        rdDepictor.SetPreferCoordGen(False)
        dopt = d2d.drawOptions()
        dopt.fixedBondLength = self.scale
        dopt.scalingFactor = self.scale
        dopt.fixedScale = True
        dopt.addAtomIndices = self.add_atom_indices
        dopt.addBondIndices = self.add_bond_indices
        dopt.dummiesAreAttachments = self.dummies_are_attachments
        dopt.padding = 0.2
        dopt.useBWAtomPalette()

        try:
          dopt.prepareMolsBeforeDrawing = True
          d2d.DrawMolecule(self.mol)
        except RuntimeError: # necessary for SMARTS strings that fail sanitization
          dopt.prepareMolsBeforeDrawing = False 
          d2d.DrawMolecule(self.mol)


        self.coords = np.array(
            [list(d2d.GetDrawCoords(i)) for i in range(self.mol.GetNumAtoms())]
        )
        d2d.FinishDrawing()

        svg = d2d.GetDrawingText()

        self.svgdom = dom = parseString(str(svg))

        # remove RDKIT namespace, because this xml is heavily modified
        self.svgdom.firstChild.removeAttribute("xmlns:rdkit")

        groups = ["shading", "mol_halo", "lines", "text", "overlay"]
        self.groups = {}
        for g in groups:
            self.groups[g] = dom.createElementNS("http://www.w3.org/2000/svg", "g")
            self.groups[g].setAttribute("class", g)

            gid = dom.createElementNS("http://www.w3.org/2000/svg", "g")
            gid.setAttribute("id", g)
            self.groups[g].appendChild(gid)

        for c in list(dom.childNodes[0].childNodes):
            if c.nodeName == "rect":
                dom.childNodes[0].removeChild(c)
                continue

            dom.childNodes[0].removeChild(c)

            if c.nodeName == "path":
                with contextlib.suppress(Exception):
                    c.removeAttribute("fill")

                if c.getAttribute("style"):
                    s = _style2dict(c.getAttribute("style"))
                    c.removeAttribute("style")

                    with contextlib.suppress(Exception):
                        del s["stroke-linecap"]
                        del s["stroke-linejoin"]
                        if s["fill"] == "none":
                            del s["fill"]

                    if "stroke-dasharray" in s:
                        c.setAttribute(
                            "style", f"stroke-dasharray:{s['stroke-dasharray']}"
                        )

                    self.groups["lines"].setAttribute("style", _dict2style(s))
                    self.groups["lines"].firstChild.appendChild(c)  # type: ignore
                else:
                    self.groups["text"].firstChild.appendChild(c)  # type: ignore
            elif not c.TEXT_NODE:
                self.groups["text"].appendChild(c)

        s = self.groups["lines"].getAttribute("style")
        s = _style2dict(s)
        s["stroke-linecap"] = "round"
        s["stroke-linejoin"] = "round"
        s["stroke-width"] = self.stroke_width = str(self.scale * 0.1)
        self.groups["lines"].setAttribute("style", _dict2style(s))

        for value in self.groups.values():
            self.svgdom.firstChild.appendChild(value)

        json.encoder.FLOAT_REPR = lambda o: format(o, ".1f")  # type: ignore
        json.encoder.c_make_encoder = None  # type: ignore

        JSON = {"coords": self.coords.tolist(), "scale": self.scale}
        JSON = json.dumps(JSON, use_decimal=True)
        if self.embed_script:
            script = dom.createElementNS("http://www.w3.org/2000/svg", "script")
            script.setAttribute("type", "application/json")
            script.appendChild(dom.createTextNode(JSON))
            self.svgdom.firstChild.appendChild(script)

        if self.optimize_svg:
            self._optimize_svg(self.svgdom)

        self.reframe()

        return

    def _optimize_svg(self, svgdom):
        _optimize_svg(svgdom)

    def get_cmap(self) -> "Colormap":
        from matplotlib import colormaps  # type: ignore

        cmap = self.cmap
        return colormaps[cmap] if isinstance(cmap, str) else cmap

    def copy(self) -> "Xenopict":
        return Xenopict(self.mol)

    def color_map(self, color):
        if self.diverging_cmap:
            color = (color + 1.0) / 2
        return self.get_cmap()(color)  # type: ignore

    def _shapely_from_atoms(
        self,
        atoms: Sequence[AtomIdx],
        bonds: Optional[Sequence[Sequence[AtomIdx]]] = None,
        twohop=False,
    ):

        atom_set = set(atoms)

        out = LineString()  # empty set

        for a in atoms:
            xy = self.coords[a]
            out = out.union(Point(*xy))  # type: ignore

        # limit to bonds provided in args (default: obtain bonds from atoms)
        _bonds: Sequence[Sequence[AtomIdx]] = bonds or [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in self.mol.GetBonds()
        ]

        # filter out all bonds where both ends are not in atom_set
        _bonds = [b for b in _bonds if b[0] in atom_set and b[1] in atom_set]

        for a1, a2 in _bonds:
            c1 = self.coords[a1]
            c2 = self.coords[a2]
            out = out.union(LineString([c1, c2]))

        if twohop:
            neighborhood = defaultdict(set)
            for a1, a2 in _bonds:
                neighborhood[a1].add(a2)
                neighborhood[a2].add(a1)

            twohop = defaultdict(set)
            for a in neighborhood:
                for n in neighborhood[a]:
                    for m in neighborhood[n]:
                        twohop[a].add(m)

            for a1, a2s in twohop.items():
                c1 = self.coords[a1]
                for a2 in a2s:
                    c2 = self.coords[a2]
                    out = out.union(LineString([c1, c2]))

        return out

    def mark_substructure(
        self,
        atoms: Sequence[AtomIdx],
        substr_bonds: Optional[Sequence[Sequence[AtomIdx]]] = None,
    ) -> "Xenopict":
        if not atoms:
            return self

        # Get the bonds to mark
        bonds_to_mark = substr_bonds if substr_bonds is not None else [
            (b.GetBeginAtomIdx(), b.GetEndAtomIdx())
            for b in self.mol.GetBonds()
            if b.GetBeginAtomIdx() in atoms and b.GetEndAtomIdx() in atoms
        ]

        # Create a path for each bond
        for bond in bonds_to_mark:
            c1 = self.coords[bond[0]]
            c2 = self.coords[bond[1]]
            line = LineString([c1, c2]).buffer(
                self.scale * self.mark_down_scale, resolution=self.shapely_resolution
            )
            d = _poly_to_path(line)

            mark = self.svgdom.createElementNS("http://www.w3.org/2000/svg", "path")
            mark.setAttribute("d", d)
            mark.setAttribute("class", f"bond-{bond[0]} atom-{bond[0]} atom-{bond[1]}")
            self._append_mark(mark)

        return self

    def _append_mark(self, mark):
        self._init_mark_layers()
        self.groups["mark"].firstChild.appendChild(mark.cloneNode(True))  # type: ignore
        # self.groups["halo"].appendChild(mark.cloneNode(True))

    def shade_substructure(
        self,
        substrs_by_atoms: Sequence[Sequence[AtomIdx]],
        shading: Sequence[float],
        substrs_bonds: Optional[Sequence[Optional[Sequence[Sequence[AtomIdx]]]]] = None,
    ) -> "Xenopict":
        """
        shade_substructure shades a list of substruture, each one defined as a list of atom idxs.

        By default, all the bonds connecting these structures are included. Optionally, a list of bonds can be provided instead.

        Args:
            substrs_by_atoms (Sequence[Sequence[AtomIdx]]): A list of substructures, each one defined as a list of atoms.

            shading (Sequence[float]): A list of shading intensities, one for each substructure.

            substrs_bonds (Optional[Sequence[Sequence[AtomIdx]]], optional): Optionally specify the bonds
            to include for each substructure. Defaults to None.

        Returns:
            Xenopict: Modifies object in place, but returns copy of self to enable chaining.
        """

        assert len(substrs_by_atoms) == len(
            shading
        ), "Number of substructures must equal number of shading values."

        dots = self.plot_dot.all_dots(shading)
        shapes = []

        if substrs_bonds:
            assert len(substrs_by_atoms) == len(
                substrs_bonds
            ), "If provided, nubmer of bond lists must equal number of substructures."
            _substrs_bonds: Sequence[
                Optional[Sequence[Sequence[AtomIdx]]]
            ] = substrs_bonds
        else:
            _substrs_bonds = [None] * len(substrs_by_atoms)

        for atoms, dot, bonds in zip(substrs_by_atoms, dots, _substrs_bonds):
            if not atoms:
                continue
            substr = self._shapely_from_atoms(atoms, bonds)
            shapes.extend((radius, color, substr) for radius, color in dot)

        self.plot_dot._sort_dots(shapes)

        for radius, color, substr in shapes:
            color = self.color_map(color)
            fill = self._color_to_style(color)
            d = _poly_to_path(
                substr.buffer(
                    self.scale * radius * 0.9, resolution=self.shapely_resolution
                )
            )

            shade = self.svgdom.createElementNS("http://www.w3.org/2000/svg", "path")
            shade.setAttribute("d", d)
            shade.setAttribute("style", f"fill:{fill}")
            self.groups["shading"].appendChild(shade)

        return self

    def __getstate__(self):
        state = self.__dict__.copy()
        del state["groups"]
        # Do not uniquify id/href when pickling
        state["svgdom"] = self.to_svg(uniquify_internal_refs=False)
        return state

    def __setstate__(self, state):
        state["svgdom"] = dom = parseString(state["svgdom"])
        state["groups"] = {
            g.getAttribute("class"): g for g in dom.getElementsByTagName("g")
        }
        self.__dict__ = state

    def _color_to_style(self, color: Sequence[float]):
        return "rgb(%g,%g,%g)" % tuple(int(x * 255) for x in color[:3])

    def mark_atoms(self, atoms: Sequence[AtomIdx]) -> "Xenopict":
        def marks():
            for a in atoms:
                xy = self.coords[a]
                circle = self._circle(
                    xy, self.scale * self.mark_down_scale, cls=f"atom-{a}"  # type: ignore
                )
                yield circle

        for mark in marks():
            self._append_mark(mark)

        return self

    def shade(
        self,
        atom_shading: Optional[AtomShading] = None,
        bond_shading: Optional[BondShading] = None,
    ) -> "Xenopict":
        scaling = self.scale * 0.9

        if atom_shading is not None and bond_shading is not None:
            scaling = self.scale * 0.8

        circles = []
        if atom_shading is not None:
            circles.extend(self.plot_dot(atom_shading, self.coords))

        if bond_shading is not None:
            atom1 = bond_shading[0]
            atom2 = bond_shading[1]
            bond_coords = (
                np.take(self.coords, atom1, axis=0)  # type: ignore
                + np.take(self.coords, atom2, axis=0)  # type: ignore
            ) / 2
            circles.extend(self.plot_dot(bond_shading[2], bond_coords))

        self.plot_dot._sort_dots(circles)

        for radius, color, xy in circles:
            c = self._circle(xy, radius * scaling, self.color_map(color))
            self.groups["shading"].appendChild(c)

        return self

    def _frame(self, padding=1.5, atoms=None):
        coords = self.coords

        if atoms is not None:
            coords = coords[atoms]
        elif self._filter:
            coords = coords[self._filter]  # type: ignore

        xy1, xy2 = coords.min(axis=0), coords.max(axis=0)  # type: ignore
        wh = xy2 - xy1 + self.scale * padding * 2
        xy1 = xy1 - self.scale * padding
        return xy1[0], xy1[1], wh[0], wh[1]

    def reframe(self, padding=1.5, atoms=None) -> "Xenopict":
        x, y, w, h = self._frame(padding, atoms)
        self.svgdom.firstChild.setAttribute(
            "viewBox",
            "%0.1f %0.1f %0.1f %0.1f" % (x, y, w, h),
        )
        self.svgdom.firstChild.setAttribute("width", "%0.1f" % w)
        self.svgdom.firstChild.setAttribute("height", "%0.1f" % h)
        return self

    def __str__(self) -> SVG:
        return self.to_html()

    def __ge__(self, other: Union["Xenopict", Mol, str]):
        if isinstance(other, Xenopict):
            patt: Mol = other.mol
        elif isinstance(other, str):
            patt: Mol = MolFromSmarts(other)
        else:
            patt = other

        hit_ats = self.mol.GetSubstructMatch(patt)

        if not hit_ats:
            return False

        self.mark_substructure(hit_ats)
        return True

    def _repr_svg_(self):
        return self.to_svg()

    def to_svg(
        self,
        uniquify_internal_refs: bool = True,
        hash_length: int = 10,
        svg_attributes: dict = {},
    ):
        """
        Convert Xenopict into an svg. By default, this function will uniquify all id/hrefs
        int the svg with the same md5 hash. This prevents id
        clashes in any documents into which svgs are embedded.
        """

        if not uniquify_internal_refs:
            return self.svgdom.toxml()

        svg = self.svgdom.toxml()

        md5 = hashlib.md5(svg.encode("utf-8")).hexdigest()[:hash_length]

        # The 're.subs' below is equivalent to this code:
        #
        # for e in dom.getElementsByTagName("g"):
        #     if i := e.getAttribute("id"):
        #         e.setAttribute("id", f"{i}_u_{md5}")
        #
        # for e in dom.getElementsByTagName("use"):
        #     if i := e.getAttribute("href"):
        #         e.setAttribute("href", f"{i}_u_{md5}")

        def addhash(matchobj):
            return f'{matchobj.group(0)[:-1]}_xeno_{md5}"'

        svg = re.sub(r'href=".+?"|id=".+?"', addhash, svg)

        # add in any SVG attribues
        if svg_attributes:
            attr = " ".join(
                [f"{key}='{value}'" for key, value in svg_attributes.items()]
            )
            svg = svg.replace("<svg ", f"<svg {attr} ", 1)

        return svg

    def to_html(self, svg_datauri=False) -> str:
        """Return the HTML string depicting the molecule, embedding the
        SVG element within a white-background styled div. Optionally,
        the SVG can be placed into an img tag's datauri isntead.
        """

        # Data URI images stripped from GitHub, so this works locally but not on github
        if svg_datauri:
            datauri = f"data:image/svg+xml;utf8,{quote(self.to_svg())}"
            return f"<div style='background:white;width:100%'><img style='display:block;max-width:100%;margin:auto' data-xenopict src={datauri} /></div>"

        # Embedding SVG directly works on GitHub, so this is the default.
        else:
            svg = self.to_svg(
                svg_attributes={"style": "display:block;max-width:100%;margin:auto"}
            )
            return f"<div style='background:white;width:100%'>{svg}</div>"

    def __getattr__(self, key):
        """
        Look up missing members by aliasing to self.mol.
        """

        if key in self.__dict__:
            return self.__dict__[key]

        if hasattr(super(), key):
            return getattr(super(), key)

        return getattr(self.mol, key)

    def halo(self) -> "Xenopict":
        # warn(
        #     "The halo method is depreciated and will be automatically applied in future version.",
        #     DeprecationWarning,
        # )

        lines = self.svgdom.createElementNS("http://www.w3.org/2000/svg", "use")
        lines.setAttribute("href", "#lines")

        text = self.svgdom.createElementNS("http://www.w3.org/2000/svg", "use")
        text.setAttribute("href", "#text")

        self.groups["mol_halo"].appendChild(lines)
        self.groups["mol_halo"].appendChild(text)

        self.groups["mol_halo"].setAttribute(
            "style",
            "stroke:white;opacity:0.5;stroke-linecap:round;stroke-linejoin:round",
        )

        text.setAttribute("style", f"stroke-width:{self.scale * .1}")
        lines.setAttribute("style", f"stroke-width:{self.scale * .2}")

        return self

    def _circle(
        self,
        xy: Sequence[float],
        radius: float,
        fill: Optional[Sequence[float]] = None,
        stroke: Optional[Sequence[float]] = None,
        style: Optional[dict[str, str]] = None,
        cls: Optional[str] = None,
    ):
        c = self.svgdom.createElementNS("http://www.w3.org/2000/svg", "circle")
        c.setAttribute("r", "%.1f" % (radius))
        c.setAttribute("cx", "%.1f" % xy[0])
        c.setAttribute("cy", "%.1f" % xy[1])

        style = style or {}
        if fill is not None:
            style["fill"] = self._color_to_style(fill)
        if stroke is not None:
            style["stroke"] = self._color_to_style(stroke)
        if style:
            c.setAttribute("style", _dict2style(style))
        if cls:
            c.setAttribute("class", cls)

        return c

    def filter(
        self, atoms: Sequence[AtomIdx], bonds: Optional[Sequence[Sequence[AtomIdx]]]
    ) -> "Xenopict":
        atom_set = set(atoms)

        elems = list(self.groups["lines"].firstChild.childNodes)  # type: ignore
        elems += list(self.groups["text"].firstChild.childNodes)  # type: ignore

        _bonds = (
            {self.mol.GetBondBetweenAtoms(*b).GetIdx() for b in bonds}
            if bonds
            else set()
        )

        for elem in elems:
            cls = set(elem.getAttribute("class").split())
            a_cls = {int(c.split("-")[1]) for c in cls if c.startswith("atom-")}

            if len(atom_set & a_cls) != len(a_cls):
                elem.parentNode.removeChild(elem)
                continue

            if bonds:

                # If not a bond, continue
                _b = {int(c.split("-")[1]) for c in cls if c.startswith("bond-")}
                if not _b:
                    continue

                # if not in provided bonds, remove element
                if len(set(_b) & _bonds) == 0:
                    elem.parentNode.removeChild(elem)

        self._filter = atoms

        return self

    def substructure_focus(
        self,
        atoms: Sequence[AtomIdx],
        substr_bonds: Optional[Sequence[Sequence[AtomIdx]]] = None,
    ) -> "Xenopict":
        if not atoms:
            return self

        self.filter(atoms, substr_bonds)
        self.reframe()
        return self

    def _init_mark_layers(self):
        if "mark" in self.groups:
            return

        self.groups["halo"] = h = self.svgdom.createElementNS(
            "http://www.w3.org/2000/svg", "use"
        )
        h.setAttribute("href", "#mark")

        # self.svgdom.createElementNS(
        #    "http://www.w3.org/2000/svg", "g"
        # )
        h.setAttribute("class", "halo")
        h.setAttribute("stroke", "#555")
        h.setAttribute("opacity", "0.45")
        h.setAttribute("style", f"fill:none;stroke-width:{self.scale * 0.2}")

        self.groups["mark"] = m = self.svgdom.createElementNS(
            "http://www.w3.org/2000/svg", "g"
        )
        m.setAttribute("class", "mark")
        m.setAttribute(
            "style", f"fill:none;stroke-width:{self.scale * 0.1};opacity:0.7"
        )

        gid = self.svgdom.createElementNS("http://www.w3.org/2000/svg", "g")
        gid.setAttribute("id", "mark")
        m.appendChild(gid)

        self.groups["overlay"].appendChild(h)
        self.groups["overlay"].appendChild(m)

    def align_to(self, template: Union[str, Mol, "Xenopict"]) -> "Xenopict":
        """Align this molecule to a template molecule using MCS.
        
        This method aligns the current molecule to match the orientation of the template
        molecule by finding their maximum common substructure. The alignment is done
        in-place, modifying the current molecule's coordinates.
        
        Args:
            template: Template molecule to align to. Can be:
                     - SMILES string
                     - RDKit Mol object
                     - Another Xenopict object
                     
        Returns:
            self for method chaining
            
        Examples:
            >>> from rdkit import Chem
            >>> # Create ethanol and align propanol to it
            >>> ethanol = Xenopict("CCO")
            >>> propanol = Xenopict("CCCO")
            >>> propanol.align_to(ethanol)  # Aligns by OH group
            <xenopict.drawer.Xenopict ...>
        """
        # Convert input to RDKit Mol if needed
        if isinstance(template, str):
            template_mol = MolFromSmiles(template)
            if template_mol is None:
                raise ValueError(f"Invalid SMILES: {template}")
        elif isinstance(template, Xenopict):
            template_mol = template.mol
        else:
            template_mol = template
            
        # Perform the alignment
        align_from_mcs(self.mol, template_mol)
        
        # Redraw the molecule with new coordinates
        self.draw_mol()
        
        return self


def _poly_to_path(shape):
    d = ""
    if hasattr(shape.boundary, "geoms"):
        bounds = shape.boundary.geoms
    else:
        bounds = [shape.boundary]

    for s in bounds:
        XY = s.coords
        d += "M %0.1f,%0.1f " % XY[0]
        for xy in XY[1:]:
            d += "L %0.1f,%0.1f " % xy
        d += "Z "
    return d


def load_ipython_extension():
    import xenopict.magic


from collections import defaultdict


def _optimize_svg(svgdom):
    for elem in svgdom.getElementsByTagName("path"):
        if elem.hasAttribute("d"):
            d = elem.getAttribute("d")
            elem.setAttribute("d", _relative_path(d))

    symb = defaultdict(list)

    N = 0
    for elem in svgdom.getElementsByTagName("path"):
        if not elem.hasAttribute("d"):
            continue
        d = elem.getAttribute("d")
        if len(d) < 30:
            continue

        x, y, s = _d_symbol(d)  # type: ignore

        symb[s].append((x, y, elem))

    if symb := [i for i in symb.items() if len(i[1]) > 1]:
        g = svgdom.createElementNS("http://www.w3.org/2000/svg", "defs")
        svgdom.firstChild.appendChild(g)

        for n, (s, xyes) in enumerate(symb):
            e = svgdom.createElementNS("http://www.w3.org/2000/svg", "path")
            n = f"s{n}"
            e.setAttribute("id", n)
            e.setAttribute("d", s)
            g.appendChild(e)

            for x, y, elem in xyes:
                e = svgdom.createElementNS("http://www.w3.org/2000/svg", "use")
                e.setAttribute("href", f"#{n}")
                e.setAttribute("x", x)
                e.setAttribute("y", y)

                for c in ["style", "class", "fill", "id"]:
                    if elem.hasAttribute(c):
                        e.setAttribute(c, elem.getAttribute(c))

                elem.parentNode.replaceChild(e, elem)


def _d_symbol(d):
    m = re.match(r"^[mM](-?[0-9\.]+)[ ,]+(-?[0-9\.]+)", d)
    if m is None:
        raise ValueError(d)
    x = m[1]
    y = m[2]
    symb = f"M0,0{d[m.span()[1]:]}"
    return x, y, symb


def _relative_path(D):
    """Converts absolute path to relative path. This is an incomplete implementation
    narrowly scoped to compress rdkit SVG depictions."""
    xy: np.ndarray = np.array([0, 0])
    D = iter(D.replace(",", " ").split())
    out = ""
    try:
        while True:
            d = next(D)
            if d == "M":
                out += "m"
                xy1 = np.array([float(next(D)), float(next(D))])
                delta = xy1 - xy  # type: ignore
                xy = xy1
                out += "%.1f %.1f" % tuple(delta)
                continue

            if d == "Q":
                out += "q"
                xy1 = np.array([float(next(D)), float(next(D).strip(","))])
                delta = xy1 - xy
                out += "%.1f %.1f " % tuple(delta)

                xy1 = np.array([float(next(D)), float(next(D))])
                delta = xy1 - xy
                xy = xy1
                out += "%.1f %.1f" % tuple(delta)
                continue

            if d == "L":
                out += "l"
                xy1 = np.array([float(next(D)), float(next(D))])
                delta = xy1 - xy
                xy = xy1
                out += "%.1f %.1f " % tuple(delta)
                continue
            
            if d == 'Z':
                out += 'z'
                continue

            raise ValueError(d)

    except StopIteration:
        return f"m{out[1:]}".replace(".0", "").strip()
