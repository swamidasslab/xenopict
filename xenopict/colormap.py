import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from PIL import Image
import contextlib
import io
import base64

from colorio.cs import ColorCoordinates, HSV, HSL
import numpy as np
import logging

logger = logging.getLogger(__name__)


class ColorInterpolator(object):
    def __init__(
        self,
        interpolation_space: str = "oklab",
        perceptual_space: str = "cam16ucs",
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
                c.convert(self.interpolation_space)
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

    def set_lightness(self, color: ColorCoordinates, lightness: float = 20):
        color = color.copy()
        color.convert(self.perceptual_space)
        color.lightness

    def perceptually_scale_swatch(self, swatch) -> ColorCoordinates:
        swatch = swatch.copy()
        swatch.convert("srgb1", mode="clip")
        swatch.convert(self.interpolation_space)

        diff = self.color_diff(swatch)
        x_scaled = np.cumsum(diff)
        x_scaled /= x_scaled.max()

        x = np.linspace(0, 1, self.resolution)

        scaled_swatch = self._interpolate(x, x_scaled, swatch.data)
        return ColorCoordinates(scaled_swatch, self.interpolation_space)

    def color_diff(self, swatch: ColorCoordinates):
        swatch_pcs = swatch.copy()
        swatch_pcs.convert(self.perceptual_space)
        swatch_pcs = swatch_pcs.data

        diff = ((swatch_pcs - np.roll(swatch_pcs, 1, axis=1)) ** 2).sum(axis=0) ** 0.5
        diff[0] = 0.0
        return diff * swatch_pcs.shape[1]

    def diverging_swatch(
        self,
        neutral_color: ColorCoordinates,
        left_colors: list[ColorCoordinates],
        right_colors: list[ColorCoordinates],
        perceptually_scale=True,
    ) -> ColorCoordinates:

        left = list(reversed(left_colors)) + [neutral_color]
        right = [neutral_color] + right_colors
        left = self.many_color_swatch(left, perceptually_scale=perceptually_scale)
        right = self.many_color_swatch(right, perceptually_scale=perceptually_scale)
        return self.concatenate_swatches([left, right])

    def concatenate_swatches(
        self, swatches: list[ColorCoordinates]
    ) -> ColorCoordinates:
        cs = swatches[0].color_space

        cat_list = [swatches[0].data]
        for s in swatches[1:]:
            s = s.copy()
            if s.color_space != cs:
                s.convert(cs)
            cat_list.append(s.data[:, 1:])

        out = np.hstack(cat_list)  # type: ignore

        return ColorCoordinates(out, cs)

    def plot(self, swatch, color_space=None):
        import proplot as pplt

        c = ColorCoordinates(swatch.data, self.interpolation_space).copy()
        c.convert(color_space or self.perceptual_space)
        d = c.data

        fig = pplt.figure(share=False)
        layout = [[1, 2, 3], [4, 5, 6]]
        axs = fig.subplots(layout)  # type: ignore
        x = np.linspace(0, 1, len(d[0]))
        for i in range(3):
            axs[i].format(title=f"dim{i+1}")
            axs[i].plot(x, d[i], "k")

        for i, j, k in zip([3, 4, 5], [0, 1, 0], [1, 2, 2]):
            axs[i].plot(d[j], d[k], "k")
            axs[i].format(xlabel=f"dim{j+1}", ylabel=f"dim{k+1}")

    def to_matlab_colormap(
        self, swatch, name: str = "test_cmap", register=False
    ) -> LinearSegmentedColormap:

        swatch.convert("srgb1", mode="clip")
        cmap = LinearSegmentedColormap.from_list(name, swatch.data.T)
        if register:
            cm.register_cmap(name, cmap)
        return cmap


_REPR_PNG_SIZE = (512, 64)

_interpolator = ColorInterpolator()._interpolate


def _repr_png(color: ColorCoordinates):
    """Generate a PNG representation of the Colormap."""

    color = color.copy()
    color.convert("srgb255", mode="clip")

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


with contextlib.suppress(NameError):
    html_formatter = get_ipython().display_formatter.formatters["text/html"]  # ignore
    html_formatter.for_type(ColorCoordinates, _repr_html)

try:
    from colorcet.sineramp import sineramp

    _IPYTHON_DISPLAY_DATA = sineramp((_REPR_PNG_SIZE[1], _REPR_PNG_SIZE[0])) / 255.0
except Exception:
    _IPYTHON_DISPLAY_DATA = np.tile(
        np.linspace(0, 1, _REPR_PNG_SIZE[0]), (_REPR_PNG_SIZE[1], 1)
    )


reds = [
    ColorCoordinates([1.0, 0, 0], "srgb1"),
    ColorCoordinates([0.2, 0, 0.0], "srgb1"),
]

blues = [
    ColorCoordinates([0.0, 0, 1], "srgb1"),
    ColorCoordinates([0.0, 0.0, 0.2], "srgb1"),
]

greens = [
    ColorCoordinates([0.0, 0.8, 0], "srgb1"),
    ColorCoordinates([0.0, 0.1, 0.0], "srgb1"),
]

yellows = [
    ColorCoordinates([0.8, 0.8, 0], "srgb1"),
    ColorCoordinates([0.1, 0.1, 0.0], "srgb1"),
]

purples = [
    ColorCoordinates([0.5, 0.0, 1], "srgb1"),
    ColorCoordinates([0.05, 0.0, 0.1], "srgb1"),
]

pinks = [
    ColorCoordinates([1.0, 0.5, 1], "srgb1"),
    ColorCoordinates([0.1, 0.05, 0.1], "srgb1"),
]

white = ColorCoordinates([1.0, 1, 1], "srgb1")


oranges = [
    ColorCoordinates([1.0, 0.6, 0], "srgb1"),
    ColorCoordinates([0.2, 0.12, 0], "srgb1"),
]
