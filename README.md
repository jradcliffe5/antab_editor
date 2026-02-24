# ANTAB Editor

This project is a small local toolkit for working with ANTAB TSYS tables.  
It gives you a GUI for interactive cleanup and editing, plus shared parsing/writing helpers.

## What’s Included

- `antab_editor.py`: Main GUI editor. Includes plotting, table editing, selection tools, smoothing, and expected TSYS fill.
- `antab_io.py`: Shared ANTAB parsing and writing helpers used by the GUI.

## Requirements

- Python 3.x
- `numpy`
- `matplotlib`
- `tkinter` (usually bundled with Python on macOS)

## Project Files

- `ek053a.antab`: Example ANTAB file.
- `sefd_values.txt`: SEFD lookup table used for expected TSYS calculation.

## GUI Features

- Plot one or more TSYS indices and select points by click or drag.
- Keep the table view synced with plot selection.
- Replace or blank values (`-99.0` is treated as missing and not plotted).
- Interpolate blanked values between valid neighbors.
- Edit station `GAIN` fields (`DPFU`, `FREQ`, `POLY`).
- Preview smoothing (Mean, Median, Gaussian) with a configurable time window.
- Apply smoothing directly to data.
- Apply expected TSYS/SEFD across the current block (and set station `DPFU` to `1`).
- Save changes in place or write to a new file.

## Quick Start

```bash
./antab_editor.py ek053a.antab
```

## Notes

- Use `-99.0` for missing values. Missing values are excluded from plotting and smoothing.
- Expected TSYS/SEFD uses `sefd_values.txt` and the station `GAIN` `FREQ`, then rewrites that station’s `DPFU` to `1`.
