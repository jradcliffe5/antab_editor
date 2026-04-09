# ANTAB Editor

This project is a small local toolkit for working with ANTAB TSYS tables.  
It gives you a GUI for interactive cleanup and editing, and a non-interactive
batch mode for scripted or headless use — both via the same script.

## What's Included

- `antab_editor.py`: GUI editor and batch processor. Pass batch flags to skip
  the GUI entirely; omit them to open the interactive editor.
- `antab_io.py`: Shared ANTAB parsing and writing helpers.

## Requirements

- Python 3.x
- `numpy`
- `matplotlib` and `tkinter` (GUI mode only — usually bundled with Python on macOS)
- `astropy` (optional — required only for `--from-fits`)

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
- Apply expected TSYS/SEFD across the current block (sets station `DPFU` to `1`).
- Apply expected TSYS/SEFD to the current selection only.
- Undo the last change (button or `Ctrl/Cmd+Z`).
- Load station, frequency, and IF parameters from one or more FITS-IDI files.
- Save changes in place or write to a new file.

## Quick Start

```bash
# Open the GUI
./antab_editor.py ek053a.antab

# Batch mode — no display needed
./antab_editor.py ek053a.antab --expected-tsys -o out.antab
```

## Examples

### Edit an existing ANTAB file (GUI)

```bash
./antab_editor.py ek053a.antab
```

Opens the GUI with `ek053a.antab` loaded. If no path is given and a single `.antab`
file exists in the current directory, it is loaded automatically.

### Generate a blank ANTAB from FITS-IDI files (GUI)

```bash
cd /path/to/experiment
./antab_editor.py
```

With no `.antab` file present, the **Generate Blank ANTAB** dialog opens automatically.

1. Click **Add files…** and select your FITS-IDI files (or let the tool auto-detect
   them if they are in the current directory or a `data/` subdirectory).
2. Click **Load from FITS** — stations, frequency, IFs, and the observation time range
   are populated from the file headers.
3. Choose which stations and band to include, set the time interval, and click
   **Generate**. The result is saved as `blank.antab` (editable before saving).
4. The editor then opens automatically with the generated file.

### Clean up noisy TSYS values (GUI)

1. Open the file and select a station block from the left panel.
2. Click or drag on the plot to select outlier points.
3. Use **Replace Selected** to overwrite them, or **Blank Selected** to set them to
   `-99.0` (missing), then **Interpolate Blanks** to fill in from neighbors.
4. Preview smoothing with the **Smoothing** controls, then **Apply Smoothing** if
   satisfied.
5. Press `Ctrl+Z` (or `Cmd+Z`) to undo any single step if needed.
6. **Save** to overwrite in place, or **Save As** to write a new file.

### Apply expected TSYS to a selection (GUI)

1. Select the points you want to replace on the plot.
2. Click **Apply Expected\*DPFU to Selection** — values are replaced with SEFD
   estimates from `sefd_values.txt` using the station frequency.
3. To apply to the entire block instead, use **Apply Expected TSYS (whole antenna)**.

### Non-interactive / scripted use

Any batch flag triggers non-interactive mode — tkinter and matplotlib are never
imported, so no display is required:

```bash
# Apply expected TSYS to all stations
./antab_editor.py input.antab --expected-tsys -o output.antab

# Smooth with a 5-minute median window
./antab_editor.py input.antab --smooth median 5 -o output.antab

# Chain: expected TSYS → interpolate blanks → smooth (applied in that order)
./antab_editor.py input.antab --expected-tsys --interpolate --smooth gaussian 5 -o output.antab

# Process only one station
./antab_editor.py input.antab --smooth mean 10 --station EF -o output.antab

# Generate a blank ANTAB from FITS-IDI files (astropy required)
./antab_editor.py --from-fits scan*.fits --interval 1 -o blank.antab

# Generate from FITS and immediately smooth
./antab_editor.py --from-fits scan*.fits --interval 1 --smooth median 5 -o output.antab

# Override stations, polarizations, or band when generating from FITS
./antab_editor.py --from-fits scan*.fits --stations EF,WB,MC --pols RL --band 6 --interval 2 -o output.antab
```

Run `./antab_editor.py --help` for the full option list.

## Notes

- Use `-99.0` for missing values. Missing values are excluded from plotting and smoothing.
- Expected TSYS/SEFD uses `sefd_values.txt` and the station `GAIN` `FREQ`, then rewrites that station's `DPFU` to `1`.
- FITS-IDI loading requires `astropy` (`pip install astropy`). It reads observation parameters (stations, frequency, IFs) but does not read TSYS data from the FITS file.
- Only one level of undo is stored; each edit overwrites the previous snapshot.
