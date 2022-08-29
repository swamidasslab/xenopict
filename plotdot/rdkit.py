from plotdot import *

from rdkit import Geometry, Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor

def drawMoleculeCoords(d2d, mol):
        coords = [ tuple(d2d.GetDrawCoords(i)) for i in range(mol.GetNumAtoms())]
        return np.array(coords)

def drawCircle(drawer, xy, radius, colour=(0.8, 0.3, 0.3), rawCoords=True):
        drawer.SetColour(colour)
        drawer.DrawArc(Geometry.Point2D(xy[0], xy[1]), radius, 0, 360, rawCoords)
    
def RKDMol2ShadedSVG(mol, shading=None, scale = 20, pcm = ColorMap(), pld = PlotDot()):


    d2d = rdMolDraw2D.MolDraw2DSVG(-1,-1)

    rdDepictor.SetPreferCoordGen(False)
    d2d.drawOptions().fixedBondLength = scale
    d2d.drawOptions().padding = 0.1
    d2d.drawOptions().useBWAtomPalette()

    # mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
    rdDepictor.Compute2DCoords(mol)

    #Hack to get molecule coordinates
    d2d.DrawMolecule(mol)
    coords = drawMoleculeCoords(d2d, mol)
    d2d.ClearDrawing()


    

    if shading is not None:
        for radius, color, xy in  pld(shading, coords):
            drawCircle(d2d, xy, scale * radius * 0.9, tuple(pcm(color)))

    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()

    return d2d.GetDrawingText()
