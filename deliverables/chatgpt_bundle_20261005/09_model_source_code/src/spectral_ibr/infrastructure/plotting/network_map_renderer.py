from pathlib import Path


class NetworkMapRenderer:
    def render(self, payload: dict[str, object], destination: Path) -> None:
        raise NotImplementedError("Render network-map figure here.")

