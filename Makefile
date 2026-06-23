.PHONY: install test lint format gate-0d figs paper-tx1

install:
	python -m pip install -e ".[dev]"

test:
	python -m pytest

lint:
	python -m ruff check src tests
	python -m mypy src

format:
	python -m ruff format src tests

gate-0d:
	spectral-ibr gate run phase0d

figs:
	spectral-ibr gate run figs

paper-tx1:
	latexmk -pdf reports/papers/tx1_controller_aware_bridge/main.tex

