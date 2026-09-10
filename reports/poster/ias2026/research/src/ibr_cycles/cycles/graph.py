from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import networkx as nx
import numpy as np

from ..actions.green import ActionGreen
from .holonomy import CycleHolonomy, holonomy


@dataclass(frozen=True)
class ActionGraph:
    """Directed action-to-action feedback graph at one frequency.

    An edge b -> a exists when the block K_ab is above threshold, that is when
    action b is visible to action a through the controller-dressed grid.
    """

    s: complex
    names: tuple[str, ...]
    weights: np.ndarray
    threshold: float

    @property
    def graph(self) -> nx.DiGraph:
        g = nx.DiGraph()
        g.add_nodes_from(self.names)
        for i, target in enumerate(self.names):
            for j, source in enumerate(self.names):
                if i == j:
                    continue
                weight = float(self.weights[i, j])
                if weight > self.threshold:
                    g.add_edge(source, target, weight=weight)
        return g

    @property
    def strongly_connected_components(self) -> list[tuple[str, ...]]:
        components = nx.strongly_connected_components(self.graph)
        return sorted(
            (tuple(sorted(component)) for component in components),
            key=lambda item: (-len(item), item),
        )

    @property
    def feedback_cores(self) -> list[tuple[str, ...]]:
        """Components of size two or more, i.e. those that can close a loop."""

        return [c for c in self.strongly_connected_components if len(c) > 1]

    @property
    def max_component_size(self) -> int:
        components = self.strongly_connected_components
        return max((len(c) for c in components), default=0)

    def simple_cycles(self, *, max_length: int = 4) -> list[tuple[str, ...]]:
        cycles = [
            tuple(cycle)
            for cycle in nx.simple_cycles(self.graph)
            if 2 <= len(cycle) <= max_length
        ]
        return sorted(cycles, key=lambda item: (len(item), item))


def build_action_graph(
    green: ActionGreen, s: complex, *, relative_threshold: float = 1e-8
) -> ActionGraph:
    """Weight every ordered action pair by the norm of its Green block."""

    k = green.k(s)
    names = green.names
    n = len(names)
    weights = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            weights[i, j] = float(np.linalg.norm(green.block(k, i, j), 2))
    largest = float(weights.max()) if weights.size else 0.0
    return ActionGraph(
        s=s,
        names=names,
        weights=weights,
        threshold=relative_threshold * max(largest, 1e-300),
    )


def rank_cycles(
    green: ActionGreen,
    s: complex,
    cycles: Sequence[Sequence[str]],
) -> list[CycleHolonomy]:
    """Rank closed cycles by the magnitude of their gauge-invariant score."""

    evaluated = [holonomy(green, cycle, s) for cycle in cycles]
    return sorted(evaluated, key=lambda item: -item.magnitude)
