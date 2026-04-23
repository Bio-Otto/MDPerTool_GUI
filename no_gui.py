"""Repository-level legacy CLI entrypoint forwarding to mdpertool.no_gui."""

from mdpertool.no_gui.no_gui import main


if __name__ == "__main__":
    raise SystemExit(main())
