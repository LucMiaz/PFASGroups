"""Entry point — run with:  python -m gui  or  pfasgroups-gui"""
import sys
import os

# Ensure the PFASGroups package root is on the path when run from the
# repository root (e.g. python -m gui).
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

from gui.app import create_app


def main():
    app, window = create_app(sys.argv)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
