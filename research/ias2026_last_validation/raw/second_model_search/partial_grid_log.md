# Second-model discovery execution log

The initial implementation launched a full 7-scale × 2-coordinate × 16-portfolio sweep. The first completed slice was `PLL scale=0.125`: 16/16 initialized and 0 unstable. It was interrupted after that slice because the inherited PowerDynamics network rebuild/linearization cost was disproportionate to the remaining gate.

The retained discovery runner uses the same physically meaningful bandwidth grid at five scales and two representative portfolios: the proper three-target stable reference and the four-target flagship candidate. The nominal 16-portfolio second-model census remains preserved in `raw/p5/p5_simplegfldc_ieee39_portfolios.csv`; the reduced policy search is not represented as a full 16-portfolio claim.
