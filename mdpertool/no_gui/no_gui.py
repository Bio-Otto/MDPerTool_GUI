"""Legacy no_gui entrypoint kept as a thin wrapper around cli_window."""

import argparse

from .cli_window import add_arguments_tu_subparsers, run_mdpertool_from_cli


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run MDPerTool energy dissipation workflow from the command line."
        )
    )
    add_arguments_tu_subparsers(parser)
    args = parser.parse_args()
    run_mdpertool_from_cli(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
