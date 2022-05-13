from typing import Dict, Any, Iterable
from tmap.core import TMAPEmbedding
import _tmap as tm
from _tmap import Minhash, LSHForest, VectorUchar, layout_from_lsh_forest
from .base_layout_generator import BaseLayoutGenerator


class BuiltinLayoutGenerator(BaseLayoutGenerator):
    def __init__(
        self,
        n_permutations: int = 2048,
        n_trees: int = 128,
        weighted: bool = False,
        weighted_sample_size: int = 128,
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

        self.n_permutations = n_permutations
        self.n_trees = n_trees
        self.weighted = weighted
        self.weighted_sample_size = weighted_sample_size

    def layout(self, X: Iterable, create_mst: bool, keep_knn: bool) -> TMAPEmbedding:
        enc = Minhash(self.n_permutations)
        lf = LSHForest(self.n_trees)

        data = []
        for x in X:
            data.append(VectorUchar(x))

        lf.batch_add(enc.batch_from_binary_array(data))
        lf.index()

        # TODO: Implement weighted
        return TMAPEmbedding(
            *layout_from_lsh_forest(
                lf,
                super().get_layout_configuration(),
                create_mst=create_mst,
                keep_knn=keep_knn,
            )
        )
