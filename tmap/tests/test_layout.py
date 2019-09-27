import pytest
import random
import tmap as tm


class TestLayout:
    def test_lf_layout(self):
        random.seed(42)
        data = []
        for _ in range(100):
            row = []
            for _ in range(10):
                row.append(random.randint(0, 20))
            data.append(tm.VectorUint(row))

        mh = tm.Minhash()
        lf = tm.LSHForest()

        lf.batch_add(mh.batch_from_sparse_binary_array(data))
        lf.index()

        x, y, s, t, gp = tm.layout_from_lsh_forest(lf)
        assert len(x) == 100
        assert len(s) == 99
