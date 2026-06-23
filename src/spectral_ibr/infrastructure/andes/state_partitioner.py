from dataclasses import dataclass


@dataclass(frozen=True)
class StatePartition:
    electrical: tuple[str, ...]
    control: tuple[str, ...]


class StatePartitioner:
    def __init__(
        self,
        control_markers: tuple[str, ...] = (
            "PLL",
            "REGC",
            "REGCP",
            "REGF",
            "REEC",
            "REPC",
            "IBR",
            "WT",
        ),
    ) -> None:
        self.control_markers = tuple(marker.upper() for marker in control_markers)

    def partition(self, state_names: tuple[str, ...]) -> StatePartition:
        electrical: list[str] = []
        control: list[str] = []
        for name in state_names:
            bucket = control if self._is_control_state(name) else electrical
            bucket.append(name)
        return StatePartition(tuple(electrical), tuple(control))

    def _is_control_state(self, state_name: str) -> bool:
        normalized = state_name.upper()
        return any(marker in normalized for marker in self.control_markers)
