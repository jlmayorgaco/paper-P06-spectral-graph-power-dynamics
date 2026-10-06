def damping_ratio(s: complex) -> float:
    """Return zeta(s) = -Re(s) / |s|."""

    magnitude = abs(s)
    if magnitude == 0.0:
        return 0.0
    return -s.real / magnitude

