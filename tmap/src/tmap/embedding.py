from typing import Iterable
import _tmap as tm
from scipy.spatial.distance import cosine as cosine_distance
from annoy import AnnoyIndex
from tmap.core import TMAPEmbedding
from tmap.layout_generators import BaseLayoutGenerator, BuiltinLayoutGenerator


def embed(
    X: Iterable,
    create_mst: bool = True,
    keep_knn: bool = False,
    index: str = "builtin",
    layout_generator: BaseLayoutGenerator = BuiltinLayoutGenerator(),
) -> TMAPEmbedding:
    if index not in ["builtin", "annoy"]:
        raise ValueError('Argument index has to be one of ("builtin", "annoy").')

    return layout_generator.layout(X, create_mst, keep_knn)

    # n = len(X)

    # coords_x: tm.VectorFloat = tm.VectorFloat()
    # coords_y: tm.VectorFloat = tm.VectorFloat()
    # edge_indices_s: tm.VectorUint = tm.VectorUint()
    # edge_indices_t: tm.VectorUint = tm.VectorUint()
    # graph_properties: tm.GraphProperties = tm.GraphProperties()

    # if index == "builtin":

    # elif "annoy":
    #     annoy = AnnoyIndex(len(X[0]), metric="angular")
    #     annoy_graph = []

    #     for i, v in enumerate(X):
    #         annoy.add_item(i, v)
    #     annoy.build(10)

    #     for i in range(n):
    #         for j in annoy.get_nns_by_item(i, 10):
    #             annoy_graph.append((i, j, cosine_distance(X[i], X[j])))

    #     (
    #         coords_x,
    #         coords_y,
    #         edge_indices_s,
    #         edge_indices_t,
    #         graph_properties,
    #     ) = tm.layout_from_edge_list(
    #         n, annoy_graph, create_mst=create_mst, keep_knn=keep_knn
    #     )

    # return TMAPEmbedding(
    #     coords_x, coords_y, edge_indices_s, edge_indices_t, graph_properties
    # )
