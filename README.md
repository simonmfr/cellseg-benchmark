# Benchmarking cell segmentation in brain tissue

## Structure of this repository
- `src`
  Function definitions
- `scripts`
  Scripts for metric calculations and segmentation algorithms
- `notebooks`
  Jupyter notebooks for development and analysis
- `archive`
  Symlink to raw MERSCOPE data on DSS
- `data`
  Symlink to processed data on DSS 

## Development
We are using `ruff` and the [ruff pre-commit hook](https://github.com/astral-sh/ruff-pre-commit) to check and format the code and docstrings

### Installation
Install ruff and pre-commit in your environment and install the pre-commit hooks for ruff defined in `.pre-commit-config.yml`
```
pip install ruff
pip install pre-commit
pre-commit install
# (optional: run against all files & fix any errors that are in your current codebase)
pre-commit run --all-files
```

### Basic usage
The ruff config is located in `pyproject.toml`. See the [ruff documentation of rules](https://docs.astral.sh/ruff/rules/) for all possible rules that we can enable / disable.
As we have installed the pre-commit hook, ruff formatting and liniting will run automatically for all changed files whenever you do git commit. If there are errors, you will get a detailled messaged of the offending code and the error. Fix the errors, add the changed file and try to commit again.

You can also manually run the ruff formatter and checker on all files with: 
```
ruff format
ruff check --fix
```
or
```
pre-commit run --all-files
```

If necessary, you can also temporarily disable all pre-commit hooks when committing by using the `--no-verify` flag with `git commit`.

