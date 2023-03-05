from typing import Iterable, Callable
import tmap as tm
from tmap.core import TMAPEmbedding
from .base_layout_generator import BaseLayoutGenerator
from scipy.spatial.distance import cosine as cosine_distance
from annoy import AnnoyIndex


class AnnoyLayoutGenerator(BaseLayoutGenerator):
    def __init__(
        self,
        n_trees: int = 10,
        metric: str = "angular",
        distance_function: Callable[[Iterable, Iterable], float] = cosine_distance,
        create_mst: bool = True,
        keep_knn: bool = False,
        k: int = 10,
        kc: int = 10,
        fme_iterations: int = 100,
        fme_threads: int = 4,
        fme_precision: int = 4,
        sl_repeats: int = 1,
        sl_extra_scaling_steps: int = 2,
        sl_scaling_min: float = 1,
        sl_scaling_max: float = 1,
        sl_scaling_type: tm.ScalingType = tm.ScalingType.RelativeToDrawing,
        mmm_repeats: int = 1,
        placer: tm.Placer = tm.Placer.Barycenter,
        merger: tm.Merger = tm.LocalBiconnected,
        merger_factor: float = 2,
        merger_adjustment: int = 0,
        node_size: float = 1 / 65,
    ) -> None:
        super().__init__(
            create_mst,
            keep_knn,
            k,
            kc,
            fme_iterations,
            fme_threads,
            fme_precision,
            sl_repeats,
            sl_extra_scaling_steps,
            sl_scaling_min,
            sl_scaling_max,
            sl_scaling_type,
            mmm_repeats,
            placer,
            merger,
            merger_factor,
            merger_adjustment,
            node_size,
        )

        self.n_trees = n_trees
        self.metric = metric
        self.distance_function = distance_function

    def layout(self, X: Iterable, create_mst: bool, keep_knn: bool) -> TMAPEmbedding:
        n = len(X)
        index = AnnoyIndex(len(X[0]), metric=self.metric)
        edge_list = []

        for i, v in enumerate(X):
            index.add_item(i, v)
        index.build(self.n_trees)

        for i in range(n):
            for j in index.get_nns_by_item(i, self.k):
                edge_list.append((i, j, self.distance_function(X[i], X[j])))

        return super().layout_from_edge_list(n, edge_list, create_mst, keep_knn)
