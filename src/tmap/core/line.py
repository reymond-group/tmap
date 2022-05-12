from dataclasses import dataclass


@dataclass
class Line:
    """A simple data type to represent lines."""

    x1: float
    x2: float
    y1: float
    y2: float
