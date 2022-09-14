from __future__ import annotations

import numpy as np
from six.moves.collections_abc import Sequence


class PlotDot:
    def __init__(self, levels=4):
        self.levels = levels
        self.stops = (np.arange(levels) + 1) / levels

    def dot_radius(self, z: float, level: int) -> float:
        z = abs(z)
        if level == 0:
            return self.stops[0] ** 0.5
        offset = 1 - self.stops[level]
        R = z - offset
        return 0 if R < self.stops[0] else R ** 0.5

    def dot_color(self, z: float, level: int) -> float:
        s = -1 if z < 0 else 1
        return z if level == 0 else s * self.stops[-level - 1]

    def single_dot(self, z: float) -> list[tuple[float, float]]:
        radius = [self.dot_radius(z, l) for l in range(self.levels)]
        color = [self.dot_color(z, l) for l in range(self.levels)]
        return [d for d in zip(radius, color) if d[0]]

    def all_dots(self, zs: Sequence[float]) -> list[list[tuple[float, float]]]:
        return [self.single_dot(z) for z in zs]

    def __call__(self, zs: Sequence[float], coords):
        """
        Input:
            coords - iterable of coordinates
            zs - iterable of z values (in range [-1, 1]).
            output_colorspace - colorspace to output colors (default: srgbhex)
            mode - colorspace conversation mode (default: clip)

        Output:
            drawing sorted list of circles in tuple form (radius, color, coords)
        """
        out = []
        for dot, coord in zip(self.all_dots(zs), coords):
            out.extend((radius, color, coord) for radius, color in dot)

        self._sort_dots(out)
        return out

    def _sort_dots(self, dots):
        dots.sort(key=lambda x: (abs(x[1]), x[0]), reverse=False)
