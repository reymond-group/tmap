from typing import Iterable
from tmap.core import TMAPEmbedding
from tmap.layout_generators import BaseLayoutGenerator, BuiltinLayoutGenerator


def embed(
    X: Iterable,
    create_mst: bool = True,
    keep_knn: bool = False,
    layout_generator: BaseLayoutGenerator = BuiltinLayoutGenerator(),
) -> TMAPEmbedding:
    return layout_generator.layout(X, create_mst, keep_knn)
