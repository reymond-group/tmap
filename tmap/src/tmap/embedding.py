from typing import Iterable
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
