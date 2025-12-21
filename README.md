# README

For a more enjoyable version of `README`, see [`report/main.pdf`](https://github.com/Al3cLee/multipole_expansion/blob/main/report/main.pdf).
That file contains both theory and implementation details.

This project uses `uv` for Python environment management. Read the documentation at:

- https://docs.astral.sh/uv/

The `uv` tool is a program that _does not depend on Python_.
This means `uv` can install Python for the user and manage both Python versions
and packages.

## TLDR

Make sure `uv` and `git` are installed. Then, starting from any directory, run these commands:

```bash
git clone https://github.com/Al3cLee/multipole_expansion
cd multipole_expansion
uv sync # This may take a little while, since python will be installed
uv run run_demo.py
uv run test_higher_orders.py
```

## Requirements

- The `uv` package manager should be installed following the [install guide](https://docs.astral.sh/uv/getting-started/installation/). On MacOS, installing via `homebrew` could result in `homebrew` attempting to install the `rust` (~600MB) language. In this case, abort and refer to the official install command.
- The python version should already be specified in `pyproject.toml`, for example `requires-python = "~=3.11.0"`.
- Before running `uv sync`, there should be no `.venv` folder in the project directory. You can check this by running `ls -la`, and if you find an existing `.venv`, remove it with `rm -rf .venv`.
- The machine should be MacOS (Intel or Apple silicon are both okay), Linux x86-64, or Windows x86-64, though Linux or MacOS is recommended. Note that Linux arm-64 is not supported by `symbolica`.

## Minimal Working Example

Navigate to the project directory and install dependencies. The installation should start by installing python 3.11.x and then proceed to install symbolica.

```bash
cd multipole_expansion
uv sync
```

Once installation completes, you can run the demonstration scripts to see the multipole expansion in action.

```bash
uv run run_demo.py
uv run test_higher_orders.py
```

The first script demonstrates the basic multipole expansion up to quadrupole order and verifies that the Taylor expansion and Q tensor formulations agree. The second script pushes the implementation to higher orders, computing multipole moments up to n=7 and beyond, showing that the recursive algorithm works for arbitrary orders.

If you are on Linux or MacOS, running

```bash
./run_demo.sh
```

will run the demo and also save the terminal output to `run_demo.txt`.

## Dive Into the Code

This project folder is organized as follows. The main package `multipole_expansion` contains the core implementation, while test scripts and documentation live at the top level.

```
multipole_expansion/
├── .venv/                  # Virtual environment (created by uv)
├── pyproject.toml          # Project configuration
├── requirements.txt        # Pip-compatible dependencies
├── multipole_expansion/    # Main package
│   ├── multipole_moments.py      # Q tensor construction
│   ├── taylor_expansion.py       # Taylor series approach
│   ├── derivatives.py            # Derivatives of 1/r
│   ├── contraction.py            # Einstein summation
│   └── ...
├── run_demo.py            # Basic demonstration
├── test_higher_orders.py  # High-order verification
└── ...
```

Although `uv run <file>` is the most convenient way to run python files, when viewing code and using the language server features of your code editor (such as "go to definition"), it is helpful to make your editor aware of the locally specified environment. For this purpose, you can run `source .venv/bin/activate`. After that, the command-line prompt should look like `(multipole-expansion) ➜  multipole_expansion`, and "go to definition" will help you jump between files.
