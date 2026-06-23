# Architecture

Dependencies point inward:

`infrastructure` and `interface` depend on `application`; `application` depends on
`domain`; `domain` does not depend on IO, plotting, ANDES, or the CLI.

Gates in `src/spectral_ibr/experiments` orchestrate use cases and write artifacts
to `outputs/`.

