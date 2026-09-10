"""Controller-mediated interaction calculus for portfolios of grid actions.

The package answers one question: is the dynamic security of a portfolio of
interventions inferable from the security of its lower-order subsets?  It does
so through the chain

    MODEL -> REDUCE -> Sigma(s) -> K(s) -> Q(s) -> CYCLES -> PROVENANCE -> REPAIR

Every step that is an exact algebraic identity is verified numerically to
machine precision; every step that is empirical carries an uncertainty budget.
"""

__version__ = "0.1.0"
