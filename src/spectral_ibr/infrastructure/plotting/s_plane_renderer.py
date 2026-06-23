from pathlib import Path


class SPlaneRenderer:
    def render(self, payload: dict[str, object], destination: Path) -> None:
        raise NotImplementedError("Render pole map figure here.")

