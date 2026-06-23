from pathlib import Path


class DampingRegionRenderer:
    def render(self, payload: dict[str, object], destination: Path) -> None:
        raise NotImplementedError("Render damping-region figure here.")

