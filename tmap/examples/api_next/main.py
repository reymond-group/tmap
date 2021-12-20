import tmap as tm
import numpy as np


def random_vectors(n: int = 10, dims=4096) -> np.ndarray:
    vecs = []

    for _ in range(n):
        vecs.append(np.random.choice(2, dims))

    return np.array(vecs)


def main():
    data = random_vectors(1000)
    te = tm.embed(data, layout_generator=tm.layout_generators.AnnoyLayoutGenerator())
    tm.plot(
        te,
        show=True,
        line_kws={"linestyle": "--", "color": "gray"},
        scatter_kws={"s": 5},
    )


if __name__ == "__main__":
    main()
