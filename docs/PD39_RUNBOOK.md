# PD39 runbook

Run from the isolated worktree `C:\tmp\pd39` with the real Julia executable:

```text
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. scripts\pd39\00_model_qualification.jl
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. scripts\pd39\01_baseline.jl
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. scripts\pd39\02_single_replacement.jl
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. scripts\pd39\03_portfolio_campaign_parallel.jl
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. scripts\pd39\04_link_outage_diagnostic.jl
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. scripts\pd39\05_summarize_results.jl
```

The first three commands are the qualification, baseline, and single-replacement stages. The fourth is the exhaustive 256-portfolio × 9-scenario primary campaign; its parallel driver uses two workers by default and resumes complete checkpoints. The fifth is the independent 46-branch diagnostic. The sixth generates the summary tables and rank associations.

Before the campaign, verify:

```text
C:\Users\walla\AppData\Local\Programs\Julia\julia.exe --project=. test\runtests.jl
```

Expected scaffold checks are 7/7. The campaign must retain its generated CSV/TOML outputs and must not overwrite the qualification, baseline, or single-replacement directories.
