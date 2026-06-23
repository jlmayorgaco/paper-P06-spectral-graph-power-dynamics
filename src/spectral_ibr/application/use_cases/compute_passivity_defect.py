from spectral_ibr.domain.value_objects.control_distance import ControlDistance


def compute_passivity_defect(dc: float, phi_c: float, passivity_defect: float) -> ControlDistance:
    return ControlDistance(dc=dc, phi_c=phi_c, passivity_defect=passivity_defect)

