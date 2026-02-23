# ANTAB Editor (Local)

Small set of tools for interactively editing ANTAB TSYS tables.

## Tools

- `antab_editor.py`: Main GUI editor (multi-plot, table view, selection tools, smoothing, expected TSYS fill).
- `antab_io.py`: Shared parsing/writing helpers used by the GUI.

## Requirements

- Python 3.x
- `numpy`
- `matplotlib`
- `tkinter` (usually bundled with Python on macOS)

## Files

- `ek053a.antab`: Example ANTAB file.
- `sefd_values.txt`: SEFD lookup table used for expected TSYS calculation.

## GUI Highlights

- Multi-index plotting and selection (rectangle or clicks)
- Table view synced with plot selection
- Replace/blank values (`-99.0` indicates missing and is not plotted)
- Smoothing preview (Mean/Median/Gaussian) with configurable time window
- Apply smoothing to data
- Apply expected TSYS from SEFD + DPFU
- Save to original file or to a new file

## Usage

```bash
./antab_editor.py ek053a.antab
```

## Notes

- Missing values should be set to `-99.0` (these are excluded from plots and smoothing).
- Expected TSYS uses `sefd_values.txt` and the `GAIN` section (`DPFU` and `FREQ`) from the ANTAB file.
