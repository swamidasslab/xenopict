from __future__ import annotations

import contextlib

from colorcet import LinearSegmentedColormap, register_cmap

from xenopict._cm_def import colormaps

for k in list(colormaps):
    colormaps[f"{k}_r"] = list(reversed(colormaps[k]))


def install_colormaps() -> list[str]:
    for k in colormaps:
        cmap = LinearSegmentedColormap.from_list(k, colormaps[k], len(colormaps[k]))
        with contextlib.suppress(NameError):
            register_cmap(k, cmap)
    return list(colormaps)
