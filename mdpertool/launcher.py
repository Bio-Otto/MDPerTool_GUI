import argparse
import sys

from mdpertool.no_gui import add_arguments_tu_subparsers, run_mdpertool_from_cli


def run_mdpertool():
    parser = argparse.ArgumentParser(description="WELCOME TO MDPERTOOL ENTRY WINDOW")
    subparsers = parser.add_subparsers(dest="mode", required=True)
    subparsers.add_parser("gui", help="Graphical User Interface")
    cli_parser = subparsers.add_parser("cli", help="Command line interface")

    add_arguments_tu_subparsers(cli_parser)
    args = parser.parse_args()

    if args.mode == "cli":
        run_mdpertool_from_cli(args)
        return

    if args.mode == "gui":
        try:
            from mdpertool.ui_main import run_mdpertool as run_gui
        except (ModuleNotFoundError, ImportError) as exc:
            err_text = str(exc).lower()
            missing_pyside = (
                isinstance(exc, ModuleNotFoundError)
                and bool(exc.name)
                and exc.name.lower().startswith("pyside2")
            )
            broken_qt_runtime = (
                "pyside2" in err_text
                or "qtuitools" in err_text
                or "undefined symbol" in err_text
                or "version qt_5" in err_text
            )

            if missing_pyside or broken_qt_runtime:
                parser.error(
                    "GUI mode requires a compatible PySide2/Qt runtime. "
                    "Use a fresh conda env with python=3.9 and install mdpertool from bio-otto/conda-forge. "
                    "Do not install PySide2 via pip into a python 3.11 env."
                )
            raise

        run_gui()
        return

    parser.print_help()
    sys.exit(1)
