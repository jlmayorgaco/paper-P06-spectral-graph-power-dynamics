from spectral_ibr.infrastructure.solvers.beyn_contour_solver import BeynContourSolver


def test_beyn_solver_adapter_exists() -> None:
    assert BeynContourSolver() is not None

