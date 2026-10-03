# Poster IAS 2026

For a request targeting one poster module, read `SECTIONS.md`, `DESIGN_SYSTEM.md`, and `poster_layout.tex` first. Edit the requested file in `sections/` and keep its fixed outer width and height. Preserve the outer `minipage`, panel height, and standalone preview block. Change text, graphics, typography, and spacing inside the module as needed. Use the roles and colors in `poster_design_system.tex` for standard text and components; local diagram labels may use a tuned size. Keep the scientific claims consistent with `generated/results.tex` and `generated/claims.tex`.

Compile the target with `./compile_section.ps1 -Section <name> -StrictLayout`, inspect its PNG, then compile `full` with `-StrictLayout` and inspect the integrated poster PNG. A clean TeX log alone does not prove that TikZ elements fit inside their box.

When the user requests a change to the entire layout, adjust `poster_layout.tex` and all affected sections together. The complete page must retain aligned left and right edges, 8 mm column gaps, and 4 mm row gaps unless the user asks to change them.
