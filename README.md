# Spatial Transcriptomics - Brain Segmentation Benchmark

## Structure of this repository
- `scripts`
  Scripts for calculating metrics and running segmentation algorithms
- `notebooks`
  Jupyter notebooks for development and analysis
- `data`
  Raw and processed data 
  - `raw`
  - `processed`

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
As we have installed the pre-commit hook, ruff formatting and liniting will run automatically for all changed files whenever you do git commit. 

- install ruff in your environment `pip install ruff`
- run ruff formatter and checker with: 
  ```
  ruff format
  ruff check --fix
  ```
- TODO: install as precommit hook

