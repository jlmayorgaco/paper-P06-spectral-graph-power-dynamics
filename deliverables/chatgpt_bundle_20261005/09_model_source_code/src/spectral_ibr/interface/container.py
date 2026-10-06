from dataclasses import dataclass
from pathlib import Path

from spectral_ibr.infrastructure.persistence.filesystem_case_repository import (
    FilesystemCaseRepository,
)
from spectral_ibr.infrastructure.persistence.json_ledger_writer import JsonLedgerWriter
from spectral_ibr.infrastructure.solvers.beyn_contour_solver import BeynContourSolver


@dataclass(frozen=True)
class Container:
    case_repository: FilesystemCaseRepository
    nep_solver: BeynContourSolver
    ledger_writer: JsonLedgerWriter


def build_container(project_root: Path = Path.cwd()) -> Container:
    return Container(
        case_repository=FilesystemCaseRepository(project_root / "data"),
        nep_solver=BeynContourSolver(),
        ledger_writer=JsonLedgerWriter(),
    )

