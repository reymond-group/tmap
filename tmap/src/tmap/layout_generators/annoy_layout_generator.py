# This is an optional feature
try:
    from typing import Iterable, Callable
    import attr
    import tmap as tm
    from tmap.core import TMAPEmbedding
    from .base_layout_generator import BaseLayoutGenerator
    from scipy.spatial.distance import cosine as cosine_distance
    from annoy import AnnoyIndex

    @attr.s(auto_attribs=True)
    class AnnoyLayoutGenerator(BaseLayoutGenerator):
        n_trees: int = 10
        metric: str = "angular"
        distance_function: Callable[[Iterable, Iterable], float] = cosine_distance

        def layout(
            self, X: Iterable, create_mst: bool, keep_knn: bool
        ) -> TMAPEmbedding:
            n = len(X)
            index = AnnoyIndex(len(X[0]), metric=self.metric)
            edge_list = []

            for i, v in enumerate(X):
                index.add_item(i, v)
            index.build(self.n_trees)

            for i in range(n):
                for j in index.get_nns_by_item(i, self.k):
                    edge_list.append((i, j, cosine_distance(X[i], X[j])))

            return self.layout_from_edge_list(n, edge_list, create_mst, keep_knn)

except ModuleNotFoundError:
    ...
