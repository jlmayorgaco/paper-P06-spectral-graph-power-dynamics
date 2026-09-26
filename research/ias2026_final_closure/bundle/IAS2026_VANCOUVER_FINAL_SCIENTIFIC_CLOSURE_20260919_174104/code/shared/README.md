# Shared reproducibility scope

The final closure root is self-contained: Python scripts resolve the campaign
root using `code/python/campaign_root.py`, and Julia scripts resolve it using
`code/julia/campaign_root.jl`. The copied raw source snapshots are inputs, not
generated replacements. The package script excludes paper/poster PDFs and
asserts that the final ZIP contains no files named `FINAL_PAPER.pdf` or
`FINAL_POSTER_DRAFT.pdf`.
