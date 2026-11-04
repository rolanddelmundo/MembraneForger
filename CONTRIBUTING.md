# Reporting a problem or contributing

Everyone taking part follows the [code of conduct](CODE_OF_CONDUCT.md).

## A build failed or looks wrong

Open an issue at https://github.com/rolanddelmundo/MembraneForger/issues and attach, from the `--out` directory:

- `membranebuilder.log`, whose last lines name the stage that failed;
- `run_manifest.json`, which records the settings, the software versions and the checksums of the inputs;
- `audit.json` and `membrane_validation.md` if the build got that far.

Say the command you ran, your GROMACS version (`gmx --version`), and whether the input is a bundled example. Please do not
attach structures you are not allowed to share; the manifest and log are usually enough to see what happened.
`membraneforger --check` with the same options, run first, answers most input problems in seconds.

## Changing the code

Install with `pip install -e ".[figures,dev]"`, then run the same checks as CI:

```bash
ruff check membraneforger
python -m unittest discover -s membraneforger/tests -t .
```

The command-line tests need GROMACS on the `PATH`. Describe in the pull request what you changed and which tests you ran.
The [module map](membraneforger/docs/modules.md) shows where each stage lives.
