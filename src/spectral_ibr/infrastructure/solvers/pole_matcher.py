from spectral_ibr.domain.entities.pole import Pole


class PoleMatcher:
    def match(self, left: list[Pole], right: list[Pole]) -> list[tuple[Pole, Pole]]:
        return list(zip(left, right, strict=False))

