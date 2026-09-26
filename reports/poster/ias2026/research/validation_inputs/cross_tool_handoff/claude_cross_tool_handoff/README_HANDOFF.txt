Files for final cross-tool validation

1) ieee39_cross_tool_validation_pack.zip
   - canonical Python reference
   - pandapower static parity
   - live/equation-equivalent ANDES dynamic parity
   - frozen gates and comparison script

2) dynamic_forest_ieee39_validation.zip
   - prior network/forest/line-sensitivity experiments
   - FROZEN Phase-E holdout in holdout_lines.json
   - seed: 20260911
   - line indices: 3, 9, 15, 17, 22, 26, 31, 32, 35, 39, 42, 44
   - these are project canonical line indices, not bus numbers

Claude should extract both to a separate validation-input directory and must not overwrite frozen research outputs.
Create a separate pandapower environment, e.g. .venv/xtool-pandapower, and do not modify the frozen ANDES environments.
