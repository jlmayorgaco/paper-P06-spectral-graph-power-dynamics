from pathlib import Path

import typer

from spectral_ibr.domain.errors import SpectralIbrError
from spectral_ibr.experiments.base_gate import GateRunConfig
from spectral_ibr.experiments.registry import get_gate, list_gates

app = typer.Typer(help="Spectral IBR research CLI.")
gate_app = typer.Typer(help="Run reproducibility gates.")


@gate_app.command("list")
def list_available_gates() -> None:
    for gate_name in list_gates():
        typer.echo(gate_name)


@gate_app.command("run")
def run_gate(
    name: str,
    case: str | None = typer.Option(None, "--case"),
    profile: str | None = typer.Option(None, "--profile"),
    config: Path | None = typer.Option(None, "--config"),
    output_root: Path = typer.Option(Path("outputs"), "--output-root"),
    seed: int = typer.Option(0, "--seed"),
) -> None:
    gate = get_gate(name)
    try:
        result = gate.run(
            GateRunConfig(
                name=name,
                case_name=case,
                profile_name=profile,
                config_path=config,
                output_root=output_root,
                seed=seed,
            )
        )
    except SpectralIbrError as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc
    typer.echo(f"{result.name}: {result.verdict}")
    if result.message:
        typer.echo(result.message)


app.add_typer(gate_app, name="gate")


if __name__ == "__main__":
    app()
