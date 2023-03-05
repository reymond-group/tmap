import pytest
import tmap as tm


class TestMinhash:
    def test_init(self):
        mh = tm.Minhash()
        assert mh is not None

    def test_from_binary_array(self):
        mh = tm.Minhash(8)
        a = mh.from_binary_array(tm.VectorUchar([1, 1, 1, 1, 0, 1]))
        b = mh.from_binary_array(tm.VectorUchar([1, 1, 1, 1, 1, 0]))
        assert len(a) == 8
        assert mh.get_distance(a, b) == 0.125

    def test_from_sparse_binary_array(self):
        mh = tm.Minhash(8)
        a = mh.from_sparse_binary_array(tm.VectorUint([6, 22, 26, 62, 626, 226622]))
        b = mh.from_sparse_binary_array(tm.VectorUint([6, 22, 26, 62, 262, 226622]))
        assert len(a) == 8
        assert round(mh.get_distance(a, b), 2) == 0.25

    def test_from_string_array(self):
        mh = tm.Minhash(8)
        a = mh.from_string_array(["a", "b", "c", "x", "y", "z"])
        b = mh.from_string_array(["a", "b", "c", "v", "y", "z"])
        assert len(a) == 8
        assert round(mh.get_distance(a, b), 3) == 0.375

    def test_from_weight_array(self):
        mh = tm.Minhash(8, 42, 64)
        a = mh.from_weight_array(tm.VectorFloat([0.2, 0.6, 0.22, 0.26, 0.62, 0.66]))
        b = mh.from_weight_array(tm.VectorFloat([0.26, 0.6, 0.22, 0.26, 0.62, 1.0]))
        assert len(a) == 128
        # assert round(mh.get_weighted_distance(a, b), 3) == 0.094
