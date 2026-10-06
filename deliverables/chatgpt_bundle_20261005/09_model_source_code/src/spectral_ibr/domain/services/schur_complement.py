import numpy as np

from spectral_ibr.domain.entities.linearized_model import LinearizedModel
from spectral_ibr.domain.entities.reduced_object import ReducedObject
from spectral_ibr.domain.value_objects.self_energy import SelfEnergy


def build_self_energy(model: LinearizedModel) -> SelfEnergy:
    ne = model.aee.shape[0]
    return SelfEnergy(
        d0=np.zeros((ne, ne), dtype=np.complex128),
        c=model.aec.astype(np.complex128),
        acc=model.acc.astype(np.complex128),
        b=model.ace.astype(np.complex128),
    )


def build_schur_reduced_object(model: LinearizedModel) -> ReducedObject:
    sigma = build_self_energy(model)
    aee = model.aee.astype(np.complex128)

    def evaluate(s: complex) -> np.ndarray:
        return s * np.eye(aee.shape[0], dtype=np.complex128) - aee - sigma(s)

    return ReducedObject(name="schur_electrical", evaluate=evaluate)

