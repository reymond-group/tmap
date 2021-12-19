import attr
from typing import Dict, Any, Iterable
from tmap.core import TMAPEmbedding
from _tmap import Minhash, LSHForest, VectorUchar, layout_from_lsh_forest
from .base_layout_generator import BaseLayoutGenerator


@attr.s(auto_attribs=True)
class BuiltinLayoutGenerator(BaseLayoutGenerator):
    n_permutations: int = 2048
    n_trees: int = 128
    weighted: bool = False
    weighted_sample_size: int = 128

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
                self.get_layout_configuration(),
                create_mst=create_mst,
                keep_knn=keep_knn,
            )
        )
