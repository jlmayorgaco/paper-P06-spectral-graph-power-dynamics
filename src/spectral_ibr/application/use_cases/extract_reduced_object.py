from spectral_ibr.application.ports.case_repository import ICaseRepository
from spectral_ibr.domain.entities.power_system_case import PowerSystemCase
from spectral_ibr.domain.entities.reduced_object import ReducedObject
from spectral_ibr.domain.services.schur_complement import build_schur_reduced_object


def extract_reduced_object(
    case: PowerSystemCase,
    repository: ICaseRepository,
) -> ReducedObject:
    return build_schur_reduced_object(repository.load(case))

