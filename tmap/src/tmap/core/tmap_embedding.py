from typing import Generator
from .line import Line


class TMAPEmbedding:
    def __init__(self, x, y, s, t, gp):
        self._x = x
        self._y = y
        self._s = s
        self._t = t
        self._gp = gp

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def s(self):
        return self._s

    @property
    def t(self):
        return self._t

    @property
    def graph_properties(self):
        return self._gp

    def get_lines(self) -> Generator[Line, None, None]:
        for s, t in zip(self.s, self.t):
            yield Line(self.x[s], self.x[t], self.y[s], self.y[t])
