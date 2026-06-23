import numpy as np

from spectral_ibr.domain.value_objects.self_energy import SelfEnergy


def test_self_energy_evaluates_resolvent() -> None:
    sigma = SelfEnergy(
        d0=np.zeros((1, 1), dtype=np.complex128),
        c=np.ones((1, 1), dtype=np.complex128),
        acc=np.array([[2.0]], dtype=np.complex128),
        b=np.ones((1, 1), dtype=np.complex128),
    )

    assert np.allclose(sigma(3.0 + 0.0j), np.array([[1.0 + 0.0j]]))

