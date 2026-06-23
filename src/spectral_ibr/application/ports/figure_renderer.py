from pathlib import Path
from typing import Protocol


class IFigureRenderer(Protocol):
    def render(self, payload: dict[str, object], destination: Path) -> None:
        """Render a figure artifact."""

