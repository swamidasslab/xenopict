

from colorio.cs import ColorCoordinates
import numpy as np
from rdkit import Chem

class PlotDot:
    def __init__(self, levels=4,):
        self.levels = levels
        self.stops = (np.arange(levels) + 1) / levels
        
    def dot_radius(self, z, level):
        z = abs(z)
        if level==0: return self.stops[0] ** 0.5
        offset = 1 - self.stops[level]
        R = z - offset 
        return 0 if R < self.stops[0] else R ** 0.5 

    def dot_color(self, z, level):
        s = -1 if z < 0 else 1
        return z if level == 0 else s * self.stops[-level - 1]

    def zcolor(self, z):
        return z

    def single_dot(self, z):
        radius = [self.dot_radius(z, l) for l in range(self.levels)]
        color = [self.dot_color(z, l) for l in range(self.levels)]
        dot = [d for d in zip(radius, color) if d[0]]
        dot = [(r, self.zcolor(c)) for r, c in dot]
        return dot

    def all_dots(self, zs):
        return [self.single_dot(z) for z in zs]

    def __call__(self, zs, coords):
        '''
        
        Input:
            coords - iterable of coordinates
            zs - iterable of z values (in range [-1, 1]).
            output_colorspace - colorspace to output colors (default: srgbhex)
            mode - colorspace conversation mode (default: clip)

        Output:
            drawing sorted list of circles in tuple form (radius, color, coords)
        '''
        out = []
        for dot, coord in zip(self.all_dots(zs), coords):
            for radius, color in dot:           
                out.append((radius, color, coord))

        out.sort(key=lambda x: abs(x[1]), reverse=False) # sort from largest to smallest radius
        return out


        


class ColorMap:
    def __init__(self,  pos_color=[0.2,0.2,1.], neg_color=[1.,0.2,0.2], zero_color=[1.,1.,1.],  colorspace="oklab", output_colorspace="srgb1"):
        self.pos_color = ColorCoordinates(pos_color, "srgb1")
        self.neg_color = ColorCoordinates(neg_color, "srgb1")
        self.zero_color = ColorCoordinates(zero_color, "srgb1")

        self.pos_color.convert(colorspace)
        self.neg_color.convert(colorspace)
        self.zero_color.convert(colorspace)

        self.colorspace = colorspace
        self.output_colorspace = output_colorspace

    def __call__(self, z):
        az = abs(z)
        c = self.neg_color if z < 0 else self.pos_color
        zcolor = c * az + self.zero_color * (1 - az)
        zcolor.convert(self.output_colorspace, mode="clip")
        return zcolor.data



    
def RKDMol2ShadedSVG(mol, shading=None, scale = 20, pcm = ColorMap(), pld = PlotDot()):
    from rdkit import Geometry, Chem
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem import rdDepictor

    d2d = rdMolDraw2D.MolDraw2DSVG(-1,-1)

    def drawMoleculeCoords(d2d, mol):
        coords = [ tuple(d2d.GetDrawCoords(i)) for i in range(mol.GetNumAtoms())]
        return np.array(coords)

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


    def drawCircle(drawer, xy, radius, colour=(0.8, 0.3, 0.3), rawCoords=True):
        drawer.SetColour(colour)
        drawer.DrawArc(Geometry.Point2D(xy[0], xy[1]), radius, 0, 360, rawCoords)

    if shading is not None:
        for radius, color, xy in  pld(shading, coords):
            drawCircle(d2d, xy, scale * radius * 0.9, tuple(pcm(color)))

    d2d.DrawMolecule(mol)
    d2d.FinishDrawing()

    return d2d.GetDrawingText()
