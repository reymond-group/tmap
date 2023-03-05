from typing import Iterable, Optional, Dict, Any
from matplotlib import pyplot as plt
from tmap.core import TMAPEmbedding
from matplotlib.axes import Axes
from tmap.helpers import set_defaults


def plot(
    tmap_embedding: TMAPEmbedding,
    ax: Optional[Axes] = None,
    draw_edges: bool = True,
    show: bool = False,
    aspect: str = "equal",
    line_kws: Dict[str, Any] = None,
    scatter_kws: Dict[str, Any] = None,
):
    if ax is None:
        fig, ax = plt.subplots()

    line_kws = {} if line_kws is None else line_kws.copy()
    scatter_kws = {} if scatter_kws is None else scatter_kws.copy()

    line_kws = set_defaults(line_kws, {"linestyle": "-"})

    if draw_edges:
        for line in tmap_embedding.get_lines():
            ax.plot([line.x1, line.x2], [line.y1, line.y2], zorder=1, **line_kws)

    ax.scatter(tmap_embedding.x, tmap_embedding.y, zorder=2, **scatter_kws)
    ax.set_aspect(aspect)

    if show:
        plt.show()

    return ax
