"""Case 4: independent center-attraction alpha sweep."""

from __future__ import annotations

if __package__ in {None, ""}:
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from cases._common import main_for_case


if __name__ == "__main__":
    raise SystemExit(main_for_case(4))
