# Blind holdout protocol and limitation

The genuine alternative-model holdout used frozen inputs, frozen predictions,
and a post-reveal label file. The prediction hash was committed before labels
were revealed. The holdout contained four portfolios and all four were stable.

The result is therefore a single-class holdout: accuracy is 4/4, while
balanced accuracy, unstable precision/recall, MCC, and Cohen kappa are
undefined. The holdout cannot demonstrate discrimination or blocker discovery.
A second preregistered holdout with both classes is required.

The historical procedural P6 holdout is retained separately and is not used as
independent evidence. No retrospective label repair is treated as a fresh
holdout.
