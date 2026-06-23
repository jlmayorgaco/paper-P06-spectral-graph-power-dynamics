from dataclasses import dataclass


@dataclass(frozen=True)
class Contour:
    """Rectangular contour for nonlinear eigenvalue search."""

    re_min: float
    re_max: float
    im_min: float
    im_max: float
    nodes: int

    def contains(self, s: complex) -> bool:
        return self.re_min <= s.real <= self.re_max and self.im_min <= s.imag <= self.im_max

