import pytest

from spectral_ibr.domain.services.pole_sensitivity import compute_pole_sensitivity


def test_pole_sensitivity_placeholder_is_explicit() -> None:
    with pytest.raises(NotImplementedError):
        compute_pole_sensitivity()

