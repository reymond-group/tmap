from typing import List, Tuple
import tmap as tm
import numpy as np
import matplotlib.pyplot as plt


# import s_gd2
#  = [0,1,2,3,4]
# J = [1,2,3,4,0]
# X = Is_gd2.layout(, J)


def random_vectors(n: int = 10, dims=4096) -> np.ndarray:
    vecs = []

    for _ in range(n):
        vecs.append(np.random.choice(2, dims))

    return np.array(vecs)


def test_graph() -> List[Tuple[int, int, float]]:
    return [
        (0, 1, 0.1),
        (1, 2, 1.0),
        (2, 0, 2.0),
        (2, 3, 0.1),
        (3, 0, 1.0),
    ]


def main():
    data = random_vectors(10000)
    te = tm.embed(data, layout_generator=tm.layout_generators.BuiltinLayoutGenerator())

    tm.plot(
        te,
        show=True,
        line_kws={"linestyle": "--", "color": "gray"},
        scatter_kws={"s": 5},
    )


if __name__ == "__main__":
    main()
