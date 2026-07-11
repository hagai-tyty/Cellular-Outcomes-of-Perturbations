"""Preview the CellFate-Rx console UI without running the real pipeline.

    python show_ui.py

Shows the progress bar, clean log lines, tables, and plain-language verdict
using representative numbers, so you can see the new output in a few seconds.
"""
from cellfate.common.console import _demo

if __name__ == "__main__":
    _demo()
