

from colorio.cs import ColorCoordinates
import numpy as np


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

        out.sort(key=lambda x: abs(x[1]), reverse=False) # sort from largest to smallest color
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

