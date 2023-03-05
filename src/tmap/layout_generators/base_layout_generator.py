import _tmap as tm
from typing import Iterable, List, Tuple

from tmap.core import TMAPEmbedding


class BaseLayoutGenerator:
    def __init__(
        self,
        create_mst: bool = True,
        keep_knn: bool = False,
        k: int = 10,
        kc: int = 10,
        fme_iterations: int = 100,
        fme_threads: int = 4,
        fme_precision: int = 4,
        sl_repeats: int = 1,
        sl_extra_scaling_steps: int = 2,
        sl_scaling_min: float = 1.0,
        sl_scaling_max: float = 1.0,
        sl_scaling_type: tm.ScalingType = tm.ScalingType.RelativeToDrawing,
        mmm_repeats: int = 1,
        placer: tm.Placer = tm.Placer.Barycenter,
        merger: tm.Merger = tm.LocalBiconnected,
        merger_factor: float = 2.0,
        merger_adjustment: int = 0,
        node_size: float = 1.0 / 65.0,
    ) -> None:
        self.create_mst = create_mst
        self.keep_knn = keep_knn
        self.k = k
        self.kc = kc
        self.fme_iterations = fme_iterations
        self.fme_threads = fme_threads
        self.fme_precision = fme_precision
        self.sl_repeats = sl_repeats
        self.sl_extra_scaling_steps = sl_extra_scaling_steps
        self.sl_scaling_min = sl_scaling_min
        self.sl_scaling_max = sl_scaling_max
        self.sl_scaling_type = sl_scaling_type
        self.mmm_repeats = mmm_repeats
        self.placer = placer
        self.merger = merger
        self.merger_factor = merger_factor
        self.merger_adjustment = merger_adjustment
        self.node_size = node_size

    def layout(self, X: Iterable, create_mst: bool, keep_knn: bool) -> TMAPEmbedding:
        raise NotImplementedError()

    def layout_from_edge_list(
        self,
        n: int,
        edge_list: List[Tuple[int, int, float]],
        create_mst: bool,
        keep_knn: bool,
    ) -> TMAPEmbedding:
        return TMAPEmbedding(
            *tm.layout_from_edge_list(
                n,
                edge_list,
                self.get_layout_configuration(),
                create_mst=create_mst,
                keep_knn=keep_knn,
            )
        )

    def get_layout_configuration(self) -> tm.LayoutConfiguration:
        cfg = tm.LayoutConfiguration()
        cfg.k = self.k
        cfg.kc = self.kc
        cfg.fme_iterations = self.fme_iterations
        cfg.fme_threads = self.fme_threads
        cfg.fme_precision = self.fme_precision
        cfg.sl_repeats = self.sl_repeats
        cfg.sl_extra_scaling_steps = self.sl_extra_scaling_steps
        cfg.sl_scaling_min = self.sl_scaling_min
        cfg.sl_scaling_max = self.sl_scaling_max
        cfg.sl_scaling_type = self.sl_scaling_type
        cfg.mmm_repeats = self.mmm_repeats
        cfg.placer = self.placer
        cfg.merger = self.merger
        cfg.merger_factor = self.merger_factor
        cfg.merger_adjustment = self.merger_adjustment
        cfg.node_size = self.node_size

        return cfg
