import pytest
import random
import tmap as tm


class TestLSHForest:
    def test_init(self):
        lf = tm.LSHForest()
        assert lf is not None

    def test_query(self):
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

        assert lf.size() == len(data)

        r = lf.query_linear_scan_by_id(0, 10)
        assert r[0][1] == 0
        assert r[1][1] == 26

    def test_knn_graph(self):
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

        f = tm.VectorUint()
        t = tm.VectorUint()
        w = tm.VectorFloat()

        lf.get_knn_graph(f, t, w, 10)
        assert len(f) == 1000
        assert t[0] == 0
        assert t[1] == 26
        assert t[2] == 36
        assert t[3] == 67
        assert t[4] == 33
        assert t[5] == 83
