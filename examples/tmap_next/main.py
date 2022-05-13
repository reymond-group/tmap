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

    for i in range(n):
        vec = np.random.choice(2, dims)

        if i % 2 == 0:
            vec[:2000] = 1

        vecs.append(vec)

    return np.array(vecs)


def main():
    data = random_vectors(100)
    te = tm.embed(
        data,
        # layout_generator=tm.layout_generators.AnnoyLayoutGenerator(),
        layout_generator=tm.layout_generators.BuiltinLayoutGenerator(),
        keep_knn=True,
    )

    tm.plot(
        te,
        show=True,
        line_kws={"linestyle": "--", "color": "gray"},
        scatter_kws={"s": 5},
    )


if __name__ == "__main__":
    main()
