from matplotlib import colormaps as cm  # type: ignore
from typing import Union
from matplotlib.colors import LinearSegmentedColormap
import contextlib
from PIL import Image
from six.moves.collections_abc import Sequence  # type: ignore
import io
import base64

from colorio.cs import ColorCoordinates, HSV, HSL
import numpy as np
import logging

logger = logging.getLogger(__name__)


red = ColorCoordinates([1.0, 0, 0], "srgb1")
blue = ColorCoordinates([0.0, 0, 1], "srgb1")
green = ColorCoordinates([0.0, 0.8, 0], "srgb1")
yellow = ColorCoordinates([0.8, 0.8, 0], "srgb1")
purple = ColorCoordinates([0.5, 0.0, 1], "srgb1")
pink = ColorCoordinates([1.0, 0.5, 1], "srgb1")
white = ColorCoordinates([1.0, 1, 1], "srgb1")
black = ColorCoordinates([0.0, 0, 0], "srgb1")
orange = ColorCoordinates([1.0, 0.6, 0], "srgb1")


def convert(color: ColorCoordinates, to_cs):
    from_cs = color.color_space
    if from_cs == to_cs:
        return color

    color = ColorCoordinates(color.data, from_cs)
    kwargs = {}
    if "rgb" in to_cs:
        kwargs["mode"] = "clip"
    color.convert(to_cs, **kwargs)
    return color


class ColorInterpolator(object):
    def __init__(
        self,
        interpolation_space: str = "oklab",
        perceptual_space: str = "srlab2",
        resolution: int = 256,
    ):
        self.interpolation_space = interpolation_space
        self.perceptual_space = perceptual_space
        self.resolution = resolution

    def many_color_swatch(
        self, colors: list[ColorCoordinates], perceptually_scale=True
    ) -> ColorCoordinates:
        swatches = []
        for c in colors:
            c = c.copy()
            if c.color_space != self.interpolation_space:
                c = convert(c, self.interpolation_space)
            swatches.append(c.data[:, None])

        x0 = np.linspace(0, 1, len(swatches))
        x = np.linspace(0, 1, self.resolution)

        swatches = np.hstack(swatches)

        swatch = self._interpolate(x, x0, swatches)
        swatch = ColorCoordinates(swatch, self.interpolation_space)

        if perceptually_scale:
            swatch = self.perceptually_scale_swatch(swatch)

        return swatch

    def _interpolate(self, x, x0, data):
        return np.array([np.interp(x, x0, data[i]) for i in range(3)])
        # from scipy.interpolate import interp1d
        # interpolator =  interp1d(x0, data, 'quadratic')
        # return interpolator(x)

    def set_lightness(self, color: ColorCoordinates, lightness: float = 50):
        color = color.copy()
        color = convert(color, self.perceptual_space)
        color.lightness

    def perceptually_scale_swatch(self, swatch: ColorCoordinates) -> ColorCoordinates:
        if len(swatch.data.shape) < 2:
            return swatch
        if swatch.data.shape[1] < 2:
            return swatch

        swatch = swatch.copy()
        swatch = convert(swatch, "srgb1")
        swatch = convert(swatch, self.interpolation_space)

        # diff = self.color_diff(swatch)
        diff = abs(swatch.data[0] - np.roll(swatch.data[0], 1))
        diff[0] = 0

        x_scaled = np.cumsum(diff)
        x_scaled /= x_scaled.max()

        x = np.linspace(0, 1, self.resolution)

        scaled_swatch = self._interpolate(x, x_scaled, swatch.data)
        return ColorCoordinates(scaled_swatch, self.interpolation_space)

    def color_diff(self, swatch: ColorCoordinates):
        swatch_pcs = swatch.copy()
        swatch_pcs = convert(swatch_pcs, self.perceptual_space)
        swatch_pcs = swatch_pcs.data

        diff = ((swatch_pcs - np.roll(swatch_pcs, 1, axis=1)) ** 2).sum(axis=0) ** 0.5
        diff[0] = 0.0
        return diff * swatch_pcs.shape[1]

    def lightness_cutoff(self, swatch: ColorCoordinates, lightness: float):
        max_lightness = convert(white, self.perceptual_space).data[0]
        swatch = swatch.copy()
        swatch = convert(swatch, self.perceptual_space)
        mask = swatch.data[0] > lightness * max_lightness
        swatch.data = swatch.data[:, mask]
        return self.perceptually_scale_swatch(swatch)

    def _as_color_list(
        self,
        colors: Union[Sequence[ColorCoordinates], ColorCoordinates],
    ) -> list[ColorCoordinates]:
        return [colors] if isinstance(colors, ColorCoordinates) else list(colors)

    def diverging_swatch(
        self,
        left_colors: Union[Sequence[ColorCoordinates], ColorCoordinates],
        right_colors: Union[Sequence[ColorCoordinates], ColorCoordinates],
        neutral_color: ColorCoordinates = white.copy(),
        extreme_color: ColorCoordinates = black.copy(),
        lightness: float = 0.30,
    ) -> ColorCoordinates:

        left = (
            [extreme_color]
            + list(reversed(self._as_color_list(left_colors)))
            + [neutral_color]
        )
        right = [neutral_color] + self._as_color_list(right_colors) + [extreme_color]

        left = self.many_color_swatch(left, perceptually_scale=False)
        right = self.many_color_swatch(right, perceptually_scale=False)

        left = self.lightness_cutoff(left, lightness)
        right = self.lightness_cutoff(right, lightness)

        return self.concatenate_swatches([left, right])

    def concatenate_swatches(
        self, swatches: list[ColorCoordinates]
    ) -> ColorCoordinates:
        cs = swatches[0].color_space

        cat_list = [swatches[0].data]
        for s in swatches[1:]:
            s = convert(s, cs)
            cat_list.append(s.data[:, 1:])

        out = np.hstack(cat_list)  # type: ignore

        return ColorCoordinates(out, cs)

    def plot(self, swatch=None, color_space=None, dim_names=("L", "a", "b")):
        import proplot as pplt

        cmap = self.to_matlab_colormap(swatch)

        c = ColorCoordinates(swatch.data, self.interpolation_space).copy()  # type: ignore
        c = convert(c, color_space or self.perceptual_space)
        d = c.data

        fig = pplt.figure(share=False, refwidth=2, refheight=1)
        layout = [[1, 2], [1, 2]]
        axs = fig.subplots(layout)  # type: ignore
        x = np.linspace(0, 1, len(d[0]))

        axs[0].format(xlabel="x", ylabel=dim_names[0])
        axs[0].scatter(x, d[0], c=cmap(x))

        for i, j, k in zip([1], [2], [1]):
            axs[i].scatter(d[j], d[k], c=cmap(x))
            axs[i].format(xlabel=dim_names[j], ylabel=dim_names[k], aspect="equal")
        return swatch

    def to_matlab_colormap(
        self, swatch, name: str = "test_cmap", register=False
    ) -> LinearSegmentedColormap:

        swatch = convert(swatch, "srgb1")
        cmap = LinearSegmentedColormap.from_list(name, swatch.data.T)
        if register:
            cm.register(cmap)
        return cmap


_REPR_PNG_SIZE = (512, 64)

_interpolator = ColorInterpolator()._interpolate


def _repr_png(color: ColorCoordinates):
    """Generate a PNG representation of the Colormap."""

    color = color.copy()
    color = convert(color, "srgb255")

    if len(color.data.shape) == 1:
        # single color
        pixels = np.zeros((_REPR_PNG_SIZE[1], _REPR_PNG_SIZE[0], 3))
        pixels[:, :] = color.data
    else:
        x0 = np.linspace(0, 1, color.data.shape[1])
        pixels = _interpolator(_IPYTHON_DISPLAY_DATA, x0, color.data)

        pixels = np.rollaxis(pixels, 0, 3)

    pixels = pixels.astype(np.uint8)
    png_bytes = io.BytesIO()
    Image.fromarray(pixels).save(png_bytes, format="png")
    return png_bytes.getvalue()


def _repr_html(color: ColorCoordinates):
    """https://github.com/matplotlib/matplotlib/blob/9b1fcf67c4228c4a2788af5bcaf0c6fde09a55bf/lib/matplotlib/colors.py"""

    color = color.copy()

    png_bytes = _repr_png(color)
    png_base64 = base64.b64encode(png_bytes).decode("ascii")

    cs = color.color_space.name

    if len(color.data.shape) == 1:
        # color.convert("srgbhex", mode="clip")
        name = f"<strong>colorio:</strong> {cs} {str(color.data)}"
    else:
        name = f"<strong>colorio:</strong> {color.data.shape[1]} samples, {cs}"

    return (
        '<div style="vertical-align: middle;">'
        f"{name}"
        "</div>"
        '<div class="cmap"><img '
        f'alt="{name} colormap" '
        f'title="{name}" '
        'style="border: 1px solid #555;" '
        f'src="data:image/png;base64,{png_base64}"></div>'
    )


with contextlib.suppress(ImportError):
    from IPython import get_ipython

    html_formatter = get_ipython().display_formatter.formatters["text/html"]  # type: ignore
    html_formatter.for_type(ColorCoordinates, _repr_html)

try:
    from colorcet.sineramp import sineramp

    _IPYTHON_DISPLAY_DATA = sineramp((_REPR_PNG_SIZE[1], _REPR_PNG_SIZE[0])) / 255.0
except ImportError:
    _IPYTHON_DISPLAY_DATA = np.tile(
        np.linspace(0, 1, _REPR_PNG_SIZE[0]), (_REPR_PNG_SIZE[1], 1)
    )


def install_colormaps():

    try:  # only install if not yet installed
        cm["xenosite"]
    except KeyError:
        _ci = ColorInterpolator("srlab2", "cam16ucs")

        _colormap = _ci.diverging_swatch(blue, red, lightness=0.5)
        _ci.to_matlab_colormap(_colormap, "xenosite_bwr", register=True)
        _ci.to_matlab_colormap(_colormap, "xenosite", register=True)

        _colormap = _ci.diverging_swatch(red, blue, lightness=0.5)
        _ci.to_matlab_colormap(_colormap, "xenosite_rwb", register=True)

        _colormap = _ci.diverging_swatch(green, purple, lightness=0.5)
        _ci.to_matlab_colormap(_colormap, "xenosite_gwp", register=True)

        _colormap = _ci.diverging_swatch(purple, orange, lightness=0.5)
        _ci.to_matlab_colormap(_colormap, "xenosite_pwo", register=True)

        _colormap = _ci.diverging_swatch(black, black, lightness=0.5)
        _ci.to_matlab_colormap(_colormap, "xenosite_kwk", register=True)
