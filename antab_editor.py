#!/usr/bin/env python3
from __future__ import annotations
import argparse
import copy
import math
import os
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

# tkinter and matplotlib are imported lazily in _import_gui_deps(), called only
# when the GUI is actually launched.  This lets the module be used non-interactively
# (e.g. CLI batch mode) without requiring a display or those packages.
tk = ttk = filedialog = None  # populated by _import_gui_deps()
matplotlib = FigureCanvasTkAgg = NavigationToolbar2Tk = Figure = RectangleSelector = None


def _import_gui_deps() -> None:
    """Import tkinter and matplotlib into module globals."""
    global tk, ttk, filedialog
    global matplotlib, FigureCanvasTkAgg, NavigationToolbar2Tk, Figure, RectangleSelector
    try:
        import tkinter as _tk
        from tkinter import ttk as _ttk, filedialog as _fd
        tk, ttk, filedialog = _tk, _ttk, _fd
    except Exception as exc:
        print(f"tkinter is required for the GUI: {exc}")
        raise SystemExit(1)
    try:
        import matplotlib as _mpl
        _mpl.use("TkAgg")
        from matplotlib.backends.backend_tkagg import (
            FigureCanvasTkAgg as _FCA, NavigationToolbar2Tk as _NT2
        )
        from matplotlib.figure import Figure as _Figure
        from matplotlib.widgets import RectangleSelector as _RS
        matplotlib, FigureCanvasTkAgg, NavigationToolbar2Tk = _mpl, _FCA, _NT2
        Figure, RectangleSelector = _Figure, _RS
    except Exception as exc:
        print(f"matplotlib is required for the GUI: {exc}")
        raise SystemExit(1)

import numpy as np

from antab_io import parse_antab, write_antab, RawSegment, TsysSegment, TsysBlock, DataRow

MISSING_VALUE = -99.0
SEFD_FILENAME = "sefd_values.txt"


def first_antab_in_cwd() -> Optional[str]:
    for name in os.listdir("."):
        if name.lower().endswith(".antab"):
            return name
    return None


def parse_time(day: str, time: str) -> float:
    h, m, s = time.split(":")
    seconds = int(h) * 3600 + int(m) * 60 + int(s)
    return int(day) + (seconds / 86400.0)

def _parse_float_list(text: str) -> List[float]:
    _, _, rhs = text.partition("=")
    rhs = rhs.replace(",", " ")
    values = []
    for item in rhs.split():
        try:
            values.append(float(item))
        except ValueError:
            continue
    return values


def _format_float_list(values: List[float]) -> str:
    return ", ".join(f"{v:g}" for v in values)


def _band_cm_to_hz(cm: float) -> float:
    c = 299792458.0
    wavelength_m = cm / 100.0
    return c / wavelength_m


def _normalize_station_key(station: str) -> str:
    text = station.strip()
    if text == "":
        return ""
    token = text.split()[0].strip(",;")
    return token.upper()


def _load_sefd_table_from_dirs(
    search_dirs: List[Optional[str]],
) -> Tuple[Dict[str, Dict[float, float]], List[Tuple[float, float]]]:
    """Load sefd_values.txt from the first directory that contains it.

    Returns (sefd_table, sefd_bands).
    """
    sefd_table: Dict[str, Dict[float, float]] = {}
    sefd_bands: List[Tuple[float, float]] = []

    candidates = [
        os.path.join(d, SEFD_FILENAME) for d in search_dirs if d is not None
    ]
    sefd_path = next((p for p in candidates if os.path.exists(p)), None)
    if sefd_path is None:
        return sefd_table, sefd_bands

    with open(sefd_path, "r", encoding="utf-8") as f:
        lines = [l.rstrip("\n") for l in f if l.strip() != ""]
    if not lines:
        return sefd_table, sefd_bands

    header = [v.strip() for v in lines[0].split("|")[1:]]
    band_cm: List[Optional[float]] = []
    for h in header:
        try:
            band_cm.append(float(h))
        except ValueError:
            band_cm.append(None)
    band_hz: List[Optional[float]] = [
        (_band_cm_to_hz(cm) if cm is not None else None) for cm in band_cm
    ]
    sefd_bands = [
        (cm, hz)
        for cm, hz in zip(band_cm, band_hz)
        if cm is not None and hz is not None
    ]
    for line in lines[1:]:
        parts = [v.strip() for v in line.split("|")]
        if not parts:
            continue
        station = parts[0]
        if station == "":
            continue
        station_key = _normalize_station_key(station)
        sefd_table.setdefault(station_key, {})
        for hz, cell in zip(band_hz, parts[1:]):
            if hz is None or cell == "":
                continue
            try:
                sefd_table[station_key][hz] = float(cell)
            except ValueError:
                continue
    return sefd_table, sefd_bands


def _parse_gain_info(path: str) -> Dict[str, Dict[str, List[float]]]:
    gain_info: Dict[str, Dict[str, List[float]]] = {}
    if not os.path.exists(path):
        return gain_info
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        parts = line.split()
        if parts and parts[0].upper() == "GAIN":
            station = _normalize_station_key(parts[1]) if len(parts) > 1 else ""
            dpfu_vals: List[float] = []
            freq_vals: List[float] = []
            poly_vals: List[float] = []
            i += 1
            while i < len(lines) and lines[i].strip() != "/":
                l = lines[i].strip()
                key = l.split("=", 1)[0].strip().upper() if l else ""
                if key == "DPFU":
                    dpfu_vals = _parse_float_list(l)
                elif key == "FREQ":
                    freq_vals = _parse_float_list(l)
                elif key == "POLY":
                    poly_vals = _parse_float_list(l)
                i += 1
            if station:
                gain_info[station] = {"dpfu": dpfu_vals, "freq": freq_vals, "poly": poly_vals}
        else:
            i += 1
    return gain_info


class AntabGui:
    def __init__(self, root: tk.Tk, path: str = None, segments=None) -> None:
        self.root = root
        self.path = path
        self.segments = segments or []
        self.blocks: List[TsysBlock] = [
            seg.block for seg in self.segments if isinstance(seg, TsysSegment)
        ]

        self.block_index = 0 if self.blocks else -1
        self.block = self.blocks[0] if self.blocks else None

        self.xs: np.ndarray = np.array([], dtype=float)
        self.rows: List[DataRow] = []
        self.series_cache: Dict[str, np.ndarray] = {}
        self.selected_masks: Dict[str, np.ndarray] = {}
        self.max_highlight_points = 20000
        self.max_table_rows = 3000
        self.table_sync_skipped = False
        self.table_sync_rows = 0
        self.gain_info = _parse_gain_info(path) if path else {}
        self.sefd_bands: List[Tuple[float, float]] = []
        self.sefd_table = self._load_sefd_table()

        self.series_artists: Dict[str, Tuple] = {}
        self.selected_scatter = None
        self.rect_selector: Optional[RectangleSelector] = None
        self.smooth_enabled = tk.BooleanVar(value=False)
        self.smooth_window_var = tk.StringVar(value="5")
        self.smooth_method_var = tk.StringVar(value="Mean")
        self.save_write_new = tk.BooleanVar(value=False)
        self.save_out_var = tk.StringVar(value="")
        self.gain_dpfu_var = tk.StringVar(value="")
        self.gain_freq_var = tk.StringVar(value="")
        self.gain_poly_var = tk.StringVar(value="")
        self.table = None
        self.row_iids: List[str] = []
        self.iid_to_row: Dict[str, int] = {}
        self.table_selection_guard = False
        self.table_guard_job = None
        self.undo_snapshot: Optional[Dict[str, Any]] = None
        self.undo_button = None

        self._build_ui()
        if self.blocks:
            self._load_block(0)

    def _load_sefd_table(self) -> Dict[str, Dict[float, float]]:
        search_dirs = [
            os.path.dirname(os.path.abspath(self.path)) if self.path else None,
            os.path.dirname(os.path.abspath(__file__)),
            os.getcwd(),
        ]
        sefd_table, sefd_bands = _load_sefd_table_from_dirs(search_dirs)
        self.sefd_bands = sefd_bands
        return sefd_table

    def _build_ui(self) -> None:
        self.root.title(
            f"ANTAB GUI Editor — {os.path.basename(self.path)}" if self.path else "ANTAB GUI Editor"
        )
        self.root.geometry("1200x750")
        self.root.update_idletasks()
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        x = max(0, (sw - 1200) // 2)
        y = max(0, (sh - 750) // 2)
        self.root.geometry(f"1200x750+{x}+{y}")

        self.main = ttk.Frame(self.root, padding=6)
        self.main.pack(fill=tk.BOTH, expand=True)

        self.left_container = ttk.Frame(self.main, width=260)
        self.left_container.pack(side=tk.LEFT, fill=tk.Y, padx=(0, 6))
        self.left_canvas = tk.Canvas(self.left_container, highlightthickness=0)
        self.left_scroll = ttk.Scrollbar(self.left_container, orient=tk.VERTICAL, command=self.left_canvas.yview)
        self.left_canvas.configure(yscrollcommand=self.left_scroll.set)
        self.left_scroll.pack(side=tk.RIGHT, fill=tk.Y)
        self.left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        self.left = ttk.Frame(self.left_canvas)
        self.left_window = self.left_canvas.create_window((0, 0), window=self.left, anchor="nw")

        def _left_on_configure(event):
            self.left_canvas.configure(scrollregion=self.left_canvas.bbox("all"))

        def _left_canvas_resize(event):
            self.left_canvas.itemconfigure(self.left_window, width=event.width)

        self.left.bind("<Configure>", _left_on_configure)
        self.left_canvas.bind("<Configure>", _left_canvas_resize)

        def _left_mousewheel(event):
            if event.delta:
                self.left_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
            elif getattr(event, "num", None) == 4:
                self.left_canvas.yview_scroll(-1, "units")
            elif getattr(event, "num", None) == 5:
                self.left_canvas.yview_scroll(1, "units")

        self.left_canvas.bind("<MouseWheel>", _left_mousewheel)
        self.left_canvas.bind("<Button-4>", _left_mousewheel)
        self.left_canvas.bind("<Button-5>", _left_mousewheel)
        self.left.bind("<MouseWheel>", _left_mousewheel)
        self.left.bind("<Button-4>", _left_mousewheel)
        self.left.bind("<Button-5>", _left_mousewheel)

        self.right = ttk.Panedwindow(self.main, orient=tk.VERTICAL)
        self.right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        block_frame = ttk.LabelFrame(self.left, text="TSYS Block")
        block_frame.pack(fill=tk.X, pady=(0, 8))
        self.block_var = tk.StringVar()
        block_names = [f"{i}: {b.station}" for i, b in enumerate(self.blocks)]
        self.block_combo = ttk.Combobox(block_frame, textvariable=self.block_var, values=block_names, state="readonly")
        if self.blocks:
            self.block_combo.current(0)
        self.block_combo.pack(fill=tk.X, padx=6, pady=6)
        self.block_combo.bind("<<ComboboxSelected>>", self._on_block_change)

        index_frame = ttk.LabelFrame(self.left, text="Indices (multi-select)")
        index_frame.pack(fill=tk.BOTH, expand=True)

        self.index_list = tk.Listbox(index_frame, selectmode=tk.EXTENDED, exportselection=False)
        self.index_list.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(6, 0), pady=6)
        self.index_list.bind("<<ListboxSelect>>", self._on_index_select)

        scrollbar = ttk.Scrollbar(index_frame, orient=tk.VERTICAL, command=self.index_list.yview)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y, padx=(0, 6), pady=6)
        self.index_list.configure(yscrollcommand=scrollbar.set)

        btn_frame = ttk.Frame(self.left)
        btn_frame.pack(fill=tk.X, pady=(8, 0))

        ttk.Button(btn_frame, text="Select All", command=self._select_all_indices).pack(fill=tk.X, pady=2)
        ttk.Button(btn_frame, text="Clear Plot", command=self._clear_indices).pack(fill=tk.X, pady=2)

        sel_frame = ttk.LabelFrame(self.left, text="Selection")
        sel_frame.pack(fill=tk.X, pady=(8, 0))
        self.add_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(sel_frame, text="Add to selection", variable=self.add_var).pack(anchor=tk.W, padx=6, pady=4)
        ttk.Button(sel_frame, text="Clear Selection", command=self._clear_selection).pack(fill=tk.X, padx=6, pady=4)

        edit_frame = ttk.LabelFrame(self.left, text="Edit Values")
        edit_frame.pack(fill=tk.X, pady=(8, 0))
        ttk.Label(edit_frame, text="Value").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.value_entry = ttk.Entry(edit_frame)
        self.value_entry.pack(fill=tk.X, padx=6)
        ttk.Button(edit_frame, text="Apply to Selection", command=self._apply_value).pack(fill=tk.X, padx=6, pady=(6, 2))
        ttk.Button(edit_frame, text="Blank Selection", command=self._blank_value).pack(fill=tk.X, padx=6, pady=(0, 6))
        ttk.Button(edit_frame, text="Interpolate Blanks", command=self._interpolate_blanks).pack(fill=tk.X, padx=6, pady=(0, 6))
        ttk.Button(edit_frame, text="Apply Expected TSYS (Whole Antenna)", command=self._apply_expected).pack(fill=tk.X, padx=6, pady=(0, 6))
        ttk.Button(edit_frame, text="Apply Expected*DPFU to Selection", command=self._apply_expected_selection).pack(fill=tk.X, padx=6, pady=(0, 6))

        gain_frame = ttk.LabelFrame(self.left, text="GAIN Values")
        gain_frame.pack(fill=tk.X, pady=(8, 0))
        ttk.Label(gain_frame, text="DPFU").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.gain_dpfu_entry = ttk.Entry(gain_frame, textvariable=self.gain_dpfu_var)
        self.gain_dpfu_entry.pack(fill=tk.X, padx=6)
        ttk.Label(gain_frame, text="FREQ (MHz)").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.gain_freq_entry = ttk.Entry(gain_frame, textvariable=self.gain_freq_var)
        self.gain_freq_entry.pack(fill=tk.X, padx=6)
        ttk.Label(gain_frame, text="POLY").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.gain_poly_entry = ttk.Entry(gain_frame, textvariable=self.gain_poly_var)
        self.gain_poly_entry.pack(fill=tk.X, padx=6)
        ttk.Button(gain_frame, text="Apply GAIN Values", command=self._apply_gain_values).pack(fill=tk.X, padx=6, pady=(6, 6))

        smooth_frame = ttk.LabelFrame(self.left, text="Smoothing")
        smooth_frame.pack(fill=tk.X, pady=(8, 0))
        ttk.Label(smooth_frame, text="Window (min)").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.smooth_entry = ttk.Entry(smooth_frame, textvariable=self.smooth_window_var)
        self.smooth_entry.pack(fill=tk.X, padx=6)
        self.smooth_entry.bind("<Return>", lambda e: self._update_plot())
        self.smooth_entry.bind("<FocusOut>", lambda e: self._update_plot())
        ttk.Label(smooth_frame, text="Method").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.smooth_method_combo = ttk.Combobox(
            smooth_frame,
            textvariable=self.smooth_method_var,
            values=["Mean", "Median", "Gaussian"],
            state="readonly",
        )
        self.smooth_method_combo.pack(fill=tk.X, padx=6)
        self.smooth_method_combo.bind("<<ComboboxSelected>>", lambda e: self._update_plot())
        ttk.Checkbutton(smooth_frame, text="Preview", variable=self.smooth_enabled,
                        command=self._update_plot).pack(anchor=tk.W, padx=6, pady=4)
        ttk.Button(smooth_frame, text="Apply Smoothing", command=self._apply_smoothing).pack(fill=tk.X, padx=6, pady=(0, 6))

        file_frame = ttk.LabelFrame(self.left, text="File")
        file_frame.pack(fill=tk.X, pady=(8, 0))
        ttk.Button(file_frame, text="Open File…", command=self._open_file).pack(fill=tk.X, padx=6, pady=(6, 2))
        ttk.Button(file_frame, text="Generate Blank ANTAB…", command=self._open_generate_blank_dialog).pack(fill=tk.X, padx=6, pady=(0, 6))

        save_frame = ttk.LabelFrame(self.left, text="Save")
        save_frame.pack(fill=tk.X, pady=(8, 0))
        self.undo_button = ttk.Button(save_frame, text="Undo Last Change", command=self._undo_last_change)
        self.undo_button.pack(fill=tk.X, pady=(0, 6))
        self._update_undo_button_state()
        ttk.Button(save_frame, text="Save", command=self._save).pack(fill=tk.X)
        ttk.Button(save_frame, text="Generate Blank ANTAB...", command=self._open_generate_blank_dialog).pack(fill=tk.X, pady=(4, 0))
        ttk.Checkbutton(save_frame, text="Save to new file", variable=self.save_write_new,
                        command=self._toggle_save_output).pack(anchor=tk.W, padx=6, pady=(6, 0))
        ttk.Label(save_frame, text="Output path").pack(anchor=tk.W, padx=6, pady=(0, 2))
        self.save_out_entry = ttk.Entry(save_frame, textvariable=self.save_out_var)
        self.save_out_entry.pack(fill=tk.X, padx=6, pady=(0, 6))
        self.save_out_entry.configure(state="disabled")

        self.status_var = tk.StringVar(value="Ready")
        self.status = ttk.Label(self.left, textvariable=self.status_var, wraplength=240)
        self.status.pack(fill=tk.X, pady=(10, 0))

        plot_frame = ttk.Frame(self.right)
        table_frame = ttk.Frame(self.right)
        self.right.add(plot_frame, weight=3)
        self.right.add(table_frame, weight=2)

        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        toolbar.update()

        self.canvas.mpl_connect("pick_event", self._on_pick)

        table_cols = ["DOY", "UT"]
        self.table = ttk.Treeview(table_frame, columns=table_cols, show="headings", selectmode="extended")
        self.table.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.table.bind("<<TreeviewSelect>>", self._on_table_select)

        table_scroll_y = ttk.Scrollbar(table_frame, orient=tk.VERTICAL, command=self.table.yview)
        table_scroll_y.pack(side=tk.RIGHT, fill=tk.Y)
        self.table.configure(yscrollcommand=table_scroll_y.set)

        table_scroll_x = ttk.Scrollbar(table_frame, orient=tk.HORIZONTAL, command=self.table.xview)
        table_scroll_x.pack(side=tk.BOTTOM, fill=tk.X)
        self.table.configure(xscrollcommand=table_scroll_x.set)

        self.root.bind_all("<Control-z>", self._undo_last_change)
        self.root.bind_all("<Command-z>", self._undo_last_change)

    def _set_status(self, text: str) -> None:
        self.status_var.set(text)

    def _make_undo_snapshot(self) -> Dict[str, Any]:
        return {
            "segments": copy.deepcopy(self.segments),
            "gain_info": copy.deepcopy(self.gain_info),
            "block_index": self.block_index,
        }

    def _save_undo_snapshot(self, snapshot: Dict[str, Any], label: str) -> None:
        snapshot["label"] = label
        self.undo_snapshot = snapshot
        self._update_undo_button_state()

    def _update_undo_button_state(self) -> None:
        if self.undo_button is None:
            return
        state = "normal" if self.undo_snapshot is not None else "disabled"
        self.undo_button.configure(state=state)

    def _undo_last_change(self, event=None):
        if self.undo_snapshot is None:
            self._set_status("Nothing to undo.")
            return "break"

        snapshot = self.undo_snapshot
        self.undo_snapshot = None
        self._update_undo_button_state()

        self.segments = snapshot["segments"]
        self.gain_info = snapshot["gain_info"]
        self.blocks = [seg.block for seg in self.segments if isinstance(seg, TsysSegment)]
        if not self.blocks:
            self._set_status("Undo failed: no TSYS blocks in restored state.")
            return "break"

        target_index = int(snapshot.get("block_index", 0))
        target_index = max(0, min(target_index, len(self.blocks) - 1))
        block_names = [f"{i}: {b.station}" for i, b in enumerate(self.blocks)]
        self.block_combo.configure(values=block_names)
        self.block_combo.current(target_index)
        self._load_block(target_index)

        label = str(snapshot.get("label", "last change"))
        self._set_status(f"Undid {label}.")
        return "break"

    def _set_table_guard(self) -> None:
        self.table_selection_guard = True
        if self.table_guard_job is None:
            self.table_guard_job = self.root.after(0, self._clear_table_guard)

    def _clear_table_guard(self) -> None:
        self.table_selection_guard = False
        self.table_guard_job = None

    def _load_block(self, index: int) -> None:
        self.block_index = index
        self.block = self.blocks[index]
        self.index_list.delete(0, tk.END)
        for name in self.block.index:
            self.index_list.insert(tk.END, name)
        self.index_list.selection_clear(0, tk.END)
        if self.block.index:
            self.index_list.selection_set(0)
        self.series_cache.clear()
        self.selected_masks.clear()

        xs_list: List[float] = []
        self.rows = []
        for row in self.block.data_rows:
            if not isinstance(row, DataRow):
                continue
            xs_list.append(parse_time(row.day, row.time))
            self.rows.append(row)
        self.xs = np.asarray(xs_list, dtype=float)
        self.selected_masks = {name: np.zeros(len(self.rows), dtype=bool) for name in self.block.index}

        self._setup_table()
        self._populate_table()
        self._load_gain_fields()
        self._update_plot()
        self._set_status(f"Loaded TSYS {self.block.station}")
        self._toggle_save_output()

    def _on_block_change(self, event=None) -> None:
        raw = self.block_combo.get().split(":", 1)[0]
        try:
            idx = int(raw)
        except ValueError:
            return
        self._load_block(idx)

    def _selected_indices(self) -> List[str]:
        selected = self.index_list.curselection()
        return [self.index_list.get(i) for i in selected]

    def _active_indices(self) -> List[str]:
        selected = self._selected_indices()
        return selected if selected else list(self.block.index)

    def _clear_selection_masks(self) -> None:
        for mask in self.selected_masks.values():
            mask[:] = False

    def _selected_row_mask(self, indices: Optional[List[str]] = None) -> np.ndarray:
        if indices is None:
            indices = self._active_indices()
        if not self.rows:
            return np.zeros(0, dtype=bool)
        mask = np.zeros(len(self.rows), dtype=bool)
        for name in indices:
            m = self.selected_masks.get(name)
            if m is not None:
                mask |= m
        return mask

    def _select_all_indices(self) -> None:
        self.index_list.selection_set(0, tk.END)
        self._update_plot()

    def _clear_indices(self) -> None:
        self.index_list.selection_clear(0, tk.END)
        self._update_plot()

    def _on_index_select(self, event=None) -> None:
        self._update_plot()

    def _build_series(self, index_name: str) -> np.ndarray:
        if index_name in self.series_cache:
            return self.series_cache[index_name]
        idx = self.block.index.index(index_name)
        ys = np.full(len(self.rows), np.nan, dtype=float)
        for i, row in enumerate(self.rows):
            if idx >= len(row.values):
                continue
            val = row.values[idx].strip()
            if val == "":
                continue
            try:
                fval = float(val)
                if fval == MISSING_VALUE:
                    continue
                ys[i] = fval
            except ValueError:
                ys[i] = math.nan
        self.series_cache[index_name] = ys
        return ys

    def _smooth_window_days(self) -> Optional[float]:
        raw = self.smooth_window_var.get().strip()
        if raw == "":
            return None
        try:
            minutes = float(raw)
        except ValueError:
            return None
        if minutes <= 0:
            return None
        return minutes / (24.0 * 60.0)

    def _default_save_output(self) -> str:
        base, ext = os.path.splitext(self.path)
        if not ext:
            ext = ".antab"
        return f"{base}_edited{ext}"

    def _save_output_path(self) -> str:
        raw = self.save_out_var.get().strip()
        if raw == "":
            return self._default_save_output()
        return raw

    def _toggle_save_output(self) -> None:
        if self.save_write_new.get():
            if self.save_out_var.get().strip() == "":
                self.save_out_var.set(self._default_save_output())
            self.save_out_entry.configure(state="normal")
        else:
            self.save_out_entry.configure(state="disabled")

    def _smooth_method_key(self) -> str:
        raw = self.smooth_method_var.get().strip().lower()
        if raw.startswith("med"):
            return "median"
        if raw.startswith("gau"):
            return "gaussian"
        return "mean"

    def _parse_list_entry(self, raw: str, field_name: str) -> Tuple[bool, Optional[List[float]]]:
        text = raw.strip()
        if text == "":
            return True, None
        values: List[float] = []
        for token in text.replace(",", " ").split():
            try:
                values.append(float(token))
            except ValueError:
                self._set_status(f"Invalid {field_name} value: {token}")
                return False, None
        if not values:
            self._set_status(f"No {field_name} values provided.")
            return False, None
        return True, values

    def _load_gain_fields(self) -> None:
        station = _normalize_station_key(self.block.station)
        info = self.gain_info.get(station, {})
        self.gain_dpfu_var.set(_format_float_list(info.get("dpfu", [])))
        self.gain_freq_var.set(_format_float_list(info.get("freq", [])))
        self.gain_poly_var.set(_format_float_list(info.get("poly", [])))

    def _set_gain_line(self, lines: List[str], start: int, end: int, key: str, values: List[float]) -> Tuple[int, bool]:
        text = _format_float_list(values)
        for idx in range(start, end):
            raw = lines[idx]
            head = raw.strip().split("=", 1)[0].strip().upper() if raw.strip() else ""
            if head == key:
                indent = raw[:len(raw) - len(raw.lstrip())]
                new_line = f"{indent}{key} = {text}"
                if lines[idx] != new_line:
                    lines[idx] = new_line
                    return end, True
                return end, False
        lines.insert(end, f"{key} = {text}")
        return end + 1, True

    def _set_station_gain_values(
        self,
        station: str,
        dpfu_vals: Optional[List[float]] = None,
        freq_vals: Optional[List[float]] = None,
        poly_vals: Optional[List[float]] = None,
    ) -> Tuple[bool, bool]:
        station = _normalize_station_key(station)
        found = False
        changed = False
        for seg in self.segments:
            if not isinstance(seg, RawSegment):
                continue
            lines = seg.lines
            i = 0
            while i < len(lines):
                parts = lines[i].strip().split()
                if len(parts) >= 2 and parts[0].upper() == "GAIN" and _normalize_station_key(parts[1]) == station:
                    found = True
                    end = i + 1
                    while end < len(lines) and lines[end].strip() != "/":
                        end += 1
                    if dpfu_vals is not None:
                        end, did = self._set_gain_line(lines, i + 1, end, "DPFU", dpfu_vals)
                        changed = changed or did
                    if freq_vals is not None:
                        end, did = self._set_gain_line(lines, i + 1, end, "FREQ", freq_vals)
                        changed = changed or did
                    if poly_vals is not None:
                        end, did = self._set_gain_line(lines, i + 1, end, "POLY", poly_vals)
                        changed = changed or did
                    i = end
                i += 1
        if not found:
            return False, False
        info = self.gain_info.setdefault(station, {"dpfu": [], "freq": [], "poly": []})
        if dpfu_vals is not None:
            info["dpfu"] = dpfu_vals
        if freq_vals is not None:
            info["freq"] = freq_vals
        if poly_vals is not None:
            info["poly"] = poly_vals
        return True, changed

    def _smooth_series(self, index_name: str, window_days: float, method: str) -> np.ndarray:
        y = self._build_series(index_name)
        if window_days is None or window_days <= 0:
            return np.full_like(y, np.nan)
        mask = np.isfinite(y)
        if not mask.any():
            return np.full_like(y, np.nan)
        x_valid = self.xs[mask]
        y_valid = y[mask]
        n = x_valid.size
        half = window_days / 2.0
        smoothed_valid = np.full(n, np.nan, dtype=float)
        if method == "mean":
            for i in range(n):
                xi = x_valid[i]
                left = np.searchsorted(x_valid, xi - half, side="left")
                right = np.searchsorted(x_valid, xi + half, side="right")
                if right <= left:
                    continue
                smoothed_valid[i] = float(np.mean(y_valid[left:right]))
        elif method == "median":
            for i in range(n):
                xi = x_valid[i]
                left = np.searchsorted(x_valid, xi - half, side="left")
                right = np.searchsorted(x_valid, xi + half, side="right")
                if right <= left:
                    continue
                smoothed_valid[i] = float(np.median(y_valid[left:right]))
        else:  # gaussian
            sigma = window_days / 2.355
            if sigma <= 0:
                return np.full_like(y, np.nan)
            for i in range(n):
                xi = x_valid[i]
                left = np.searchsorted(x_valid, xi - half, side="left")
                right = np.searchsorted(x_valid, xi + half, side="right")
                if right <= left:
                    continue
                dt = x_valid[left:right] - xi
                w = np.exp(-0.5 * (dt / sigma) ** 2)
                wsum = float(np.sum(w))
                if wsum <= 0:
                    continue
                smoothed_valid[i] = float(np.sum(w * y_valid[left:right]) / wsum)
        smoothed = np.full_like(y, np.nan)
        smoothed[mask] = smoothed_valid
        return smoothed

    def _dpfu_for_index(self, dpfu_vals: List[float], index_name: str) -> Optional[float]:
        if not dpfu_vals:
            return None
        if len(dpfu_vals) == 1:
            return dpfu_vals[0]
        upper = index_name.upper()
        is_r = "R" in upper
        is_l = "L" in upper
        if is_r and not is_l:
            return dpfu_vals[0]
        if is_l and not is_r:
            return dpfu_vals[1] if len(dpfu_vals) > 1 else dpfu_vals[0]
        return dpfu_vals[0]

    def _expected_tsys_for_index(self, index_name: str, include_dpfu: bool = True) -> Optional[float]:
        station_key = _normalize_station_key(self.block.station)
        gain = self.gain_info.get(station_key)
        if not gain:
            return None
        dpfu_vals = gain.get("dpfu", [])
        freq_vals = gain.get("freq", [])
        if not freq_vals:
            return None
        freq_mhz = sum(freq_vals) / len(freq_vals)
        freq_hz = freq_mhz * 1e6
        sefd_map = self.sefd_table.get(station_key)
        if not sefd_map:
            return None
        nearest = min(sefd_map.keys(), key=lambda f: abs(f - freq_hz))
        sefd_jy = sefd_map[nearest]
        if not include_dpfu:
            return sefd_jy
        dpfu = self._dpfu_for_index(dpfu_vals, index_name)
        if dpfu is None:
            return None
        return dpfu * sefd_jy

    def _set_station_dpfu_to_one(self, station: str) -> bool:
        station = _normalize_station_key(station)
        dpfu_count = max(1, len(self.gain_info.get(station, {}).get("dpfu", [])))
        _, changed = self._set_station_gain_values(station, dpfu_vals=[1.0] * dpfu_count)
        self._load_gain_fields()
        return changed

    def _apply_gain_values(self) -> None:
        station = _normalize_station_key(self.block.station)
        ok, dpfu_vals = self._parse_list_entry(self.gain_dpfu_var.get(), "DPFU")
        if not ok:
            return
        ok, freq_vals = self._parse_list_entry(self.gain_freq_var.get(), "FREQ")
        if not ok:
            return
        ok, poly_vals = self._parse_list_entry(self.gain_poly_var.get(), "POLY")
        if not ok:
            return
        if dpfu_vals is None and freq_vals is None and poly_vals is None:
            self._set_status("Enter at least one GAIN field to apply.")
            return
        snapshot = self._make_undo_snapshot()
        found, changed = self._set_station_gain_values(
            station,
            dpfu_vals=dpfu_vals,
            freq_vals=freq_vals,
            poly_vals=poly_vals,
        )
        if not found:
            self._set_status(f"No GAIN block found for station {station}.")
            return
        self._load_gain_fields()
        if changed:
            self._save_undo_snapshot(snapshot, "GAIN edit")
            self._set_status(f"Updated GAIN values for {station}.")
        else:
            self._set_status(f"GAIN values for {station} already match.")

    def _apply_expected(self) -> None:
        indices = list(self.block.index)
        if not indices or not self.rows:
            self._set_status("No TSYS data available for this block.")
            return
        station = _normalize_station_key(self.block.station)
        snapshot = self._make_undo_snapshot()
        dpfu_changed = self._set_station_dpfu_to_one(station)
        rows_idx = np.arange(len(self.rows), dtype=int)
        rows_touched = set()
        applied_indices = 0
        for index_name in indices:
            expected = self._expected_tsys_for_index(index_name)
            if expected is None:
                continue
            idx = self.block.index.index(index_name)
            for row_idx in rows_idx:
                row = self.rows[int(row_idx)]
                self._ensure_value_length(row, idx)
                row.values[idx] = f"{expected:.1f}"
                rows_touched.add(int(row_idx))
            if index_name in self.series_cache:
                self.series_cache[index_name][rows_idx] = expected
            applied_indices += 1
        if applied_indices == 0:
            if dpfu_changed:
                self._save_undo_snapshot(snapshot, "expected TSYS apply (whole antenna)")
            self._set_status("No expected TSYS available (missing GAIN/SEFD data).")
            return
        self._save_undo_snapshot(snapshot, "expected TSYS apply (whole antenna)")
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Replaced all values with expected TSYS for {applied_indices} indices; set GAIN {station} DPFU to 1.")

    def _apply_expected_selection(self) -> None:
        total_selected = int(sum(mask.sum() for mask in self.selected_masks.values()))
        if total_selected == 0:
            self._set_status("No points selected.")
            return

        snapshot = None
        rows_touched = set()
        total_applied = 0
        applied_indices = 0
        for index_name, mask in self.selected_masks.items():
            if not mask.any():
                continue

            # expected(TSYS) * DPFU, using R/L assignment from index name.
            expected = self._expected_tsys_for_index(index_name, include_dpfu=True)
            if expected is None:
                continue

            rows_idx = np.where(mask)[0]
            if rows_idx.size == 0:
                continue
            if snapshot is None:
                snapshot = self._make_undo_snapshot()

            idx = self.block.index.index(index_name)
            for row_idx in rows_idx:
                row = self.rows[int(row_idx)]
                self._ensure_value_length(row, idx)
                row.values[idx] = f"{expected:.1f}"
                rows_touched.add(int(row_idx))
            if index_name in self.series_cache:
                self.series_cache[index_name][rows_idx] = expected
            total_applied += int(rows_idx.size)
            applied_indices += 1

        if total_applied == 0:
            self._set_status("No expected TSYS*DPFU available for selected points (missing GAIN/SEFD data).")
            return

        self._save_undo_snapshot(snapshot, "expected TSYS*DPFU apply (selection)")
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Applied expected TSYS*DPFU to {total_applied} selected points across {applied_indices} indices.")

    def _setup_table(self) -> None:
        if self.table is None:
            return
        cols = ["DOY", "UT"] + list(self.block.index)
        self.table.configure(columns=cols)
        for col in cols:
            self.table.heading(col, text=col)
            if col in ("DOY", "UT"):
                self.table.column(col, width=70, anchor=tk.W, stretch=False)
            else:
                self.table.column(col, width=70, anchor=tk.E, stretch=False)

    def _row_values(self, row: DataRow) -> List[str]:
        values = [row.day, row.time]
        for i in range(len(self.block.index)):
            values.append(row.values[i] if i < len(row.values) else "")
        return values

    def _populate_table(self) -> None:
        if self.table is None:
            return
        self._set_table_guard()
        for iid in self.table.get_children():
            self.table.delete(iid)
        self.row_iids = []
        self.iid_to_row = {}
        for idx, row in enumerate(self.rows):
            iid = self.table.insert("", tk.END, values=self._row_values(row))
            self.row_iids.append(iid)
            self.iid_to_row[iid] = idx

    def _refresh_table_rows(self, row_indices: Set[int]) -> None:
        if self.table is None:
            return
        for row_idx in row_indices:
            if row_idx < 0 or row_idx >= len(self.rows):
                continue
            iid = self.row_iids[row_idx]
            self.table.item(iid, values=self._row_values(self.rows[row_idx]))

    def _update_plot(self) -> None:
        self.ax.clear()
        self.series_artists.clear()

        selected = self._selected_indices()
        if not selected:
            self.ax.set_title("Select indices to plot")
            self._clear_selection_masks()
            if self.selected_scatter is not None:
                self.selected_scatter.set_offsets(np.empty((0, 2)))
            self._sync_table_selection()
            self.canvas.draw_idle()
            return

        active = set(selected)
        for name, mask in self.selected_masks.items():
            if name not in active:
                mask[:] = False

        window_days = self._smooth_window_days() if self.smooth_enabled.get() else None
        method = self._smooth_method_key()

        valid_x_mask = np.zeros(len(self.xs), dtype=bool)
        y_values = []
        for index_name in selected:
            ys = self._build_series(index_name)
            if np.isfinite(ys).any():
                valid_x_mask |= np.isfinite(ys)
                y_values.append(ys[np.isfinite(ys)])
            scatter = self.ax.scatter(self.xs, ys, s=12, picker=5, label=index_name, marker="o", alpha=0.9)
            scatter.set_gid(index_name)
            self.series_artists[index_name] = (None, scatter)
            if window_days is not None:
                smooth = self._smooth_series(index_name, window_days, method)
                if np.isfinite(smooth).any():
                    facecolors = scatter.get_facecolors()
                    color = facecolors[0] if facecolors is not None and len(facecolors) > 0 else None
                    self.ax.plot(self.xs, smooth, linestyle="--", linewidth=1.2, alpha=0.9,
                                 color=color, label=f"_{index_name}_smooth")

        self.ax.set_xlabel("Day of Year (fractional)")
        self.ax.set_ylabel("TSYS")
        self.ax.set_title(f"TSYS {self.block.station}")
        if valid_x_mask.any() and y_values:
            xs_valid = self.xs[valid_x_mask]
            y_concat = np.concatenate(y_values)
            x_min = float(np.min(xs_valid))
            x_max = float(np.max(xs_valid))
            y_min = float(np.min(y_concat))
            y_max = float(np.max(y_concat))
            x_pad = (x_max - x_min) * 0.02 if x_max > x_min else 0.5
            y_pad = (y_max - y_min) * 0.05 if y_max > y_min else 1.0
            self.ax.set_xlim(x_min - x_pad, x_max + x_pad)
            self.ax.set_ylim(y_min - y_pad, y_max + y_pad)
        self.ax.legend(loc="upper right", fontsize="small", ncol=2)

        self.selected_scatter = self.ax.scatter([], [], s=60, facecolors="none", edgecolors="red", linewidths=1.2)
        self._update_selected_scatter(redraw=False)

        if self.rect_selector is None:
            props = dict(facecolor="none", edgecolor="black", linestyle="--", linewidth=1.0, alpha=0.9)
            try:
                self.rect_selector = RectangleSelector(self.ax, self._on_rectangle,
                                                       useblit=True,
                                                       button=[1],
                                                       interactive=False,
                                                       drag_from_anywhere=True,
                                                       spancoords="data",
                                                       props=props)
            except TypeError:
                self.rect_selector = RectangleSelector(self.ax, self._on_rectangle,
                                                       useblit=True,
                                                       button=[1],
                                                       interactive=False,
                                                       drag_from_anywhere=True,
                                                       spancoords="data",
                                                       rectprops=props)
            self.rect_selector.set_active(True)
        self.canvas.draw_idle()

    def _on_rectangle(self, eclick, erelease) -> None:
        if eclick.xdata is None or erelease.xdata is None:
            return
        if eclick.ydata is None or erelease.ydata is None:
            return
        indices = self._selected_indices()
        if not indices:
            return

        x1, x2 = sorted([eclick.xdata, erelease.xdata])
        y1, y2 = sorted([eclick.ydata, erelease.ydata])

        if not self.add_var.get():
            self._clear_selection_masks()

        xs = self.xs
        mask_x = (xs >= x1) & (xs <= x2)
        for index_name in indices:
            ys = self._build_series(index_name)
            mask = mask_x & (ys >= y1) & (ys <= y2) & ~np.isnan(ys)
            if self.add_var.get():
                self.selected_masks[index_name] |= mask
            else:
                self.selected_masks[index_name][:] = mask

        self._update_selected_scatter()


    def _on_pick(self, event) -> None:
        artist = getattr(event, "artist", None)
        if artist is None:
            return
        index_name = artist.get_gid()
        if not index_name:
            return
        if not self.add_var.get():
            self._clear_selection_masks()
        mask = self.selected_masks.get(index_name)
        if mask is None:
            return
        for i in event.ind:
            idx = int(i)
            if 0 <= idx < mask.size:
                mask[idx] = ~mask[idx]
        self._update_selected_scatter()

    def _update_selected_scatter(self, redraw: bool = True) -> None:
        if self.selected_scatter is None:
            return
        indices = self._selected_indices()
        if not indices:
            indices = list(self.selected_masks.keys())

        count = int(sum(self.selected_masks[name].sum() for name in indices if name in self.selected_masks))
        if count == 0:
            offsets = np.empty((0, 2))
        elif count > self.max_highlight_points:
            offsets = np.empty((0, 2))
        else:
            xs_sel = []
            ys_sel = []
            for index_name in indices:
                mask = self.selected_masks.get(index_name)
                if mask is None or not mask.any():
                    continue
                ys = self._build_series(index_name)
                rows_arr = np.where(mask)[0]
                if rows_arr.size == 0:
                    continue
                yvals = ys[rows_arr]
                ok = ~np.isnan(yvals)
                if not ok.any():
                    continue
                xs_sel.extend(self.xs[rows_arr[ok]].tolist())
                ys_sel.extend(yvals[ok].tolist())
            offsets = np.column_stack([xs_sel, ys_sel]) if xs_sel else np.empty((0, 2))
        self.selected_scatter.set_offsets(offsets)
        rows_count, table_synced = self._sync_table_selection()
        if count == 0:
            status = "Selected points: 0"
        else:
            status = f"Selected points: {count}"
            if count > self.max_highlight_points:
                status += " (highlight skipped)"
        if rows_count and not table_synced:
            status += f" | rows: {rows_count} (table sync skipped)"
        self._set_status(status)
        if redraw:
            self.canvas.draw_idle()

    def _sync_table_selection(self) -> Tuple[int, bool]:
        if self.table is None:
            return 0, True
        indices = self._selected_indices()
        if not indices:
            indices = list(self.selected_masks.keys())
        rows_mask = self._selected_row_mask(indices)
        rows_selected = np.where(rows_mask)[0]
        rows_count = int(rows_selected.size)
        if rows_count > self.max_table_rows:
            return rows_count, False
        self._set_table_guard()
        iids = [self.row_iids[row_idx] for row_idx in rows_selected if 0 <= row_idx < len(self.row_iids)]
        self.table.selection_set(iids)
        return rows_count, True

    def _on_table_select(self, event=None) -> None:
        if self.table is None or self.table_selection_guard:
            return
        selected_iids = self.table.selection()
        row_indices = [self.iid_to_row[iid] for iid in selected_iids if iid in self.iid_to_row]
        index_names = self._selected_indices()
        if not index_names:
            self._set_status("Select indices to plot before selecting rows.")
            return
        rows_arr = np.asarray(row_indices, dtype=int)
        if rows_arr.size == 0:
            return
        if not self.add_var.get():
            self._clear_selection_masks()
        for index_name in index_names:
            mask = self.selected_masks.get(index_name)
            if mask is None:
                continue
            rows_ok = rows_arr[(rows_arr >= 0) & (rows_arr < mask.size)]
            if rows_ok.size:
                mask[rows_ok] = True
        self._update_selected_scatter()

    def _ensure_value_length(self, row: DataRow, idx: int) -> None:
        while len(row.values) <= idx:
            row.values.append("")

    def _apply_value(self) -> None:
        value = self.value_entry.get().strip()
        if value == "":
            self._set_status("Value is empty. Use 'Blank Selection' to blank values.")
            return
        total = int(sum(mask.sum() for mask in self.selected_masks.values()))
        if total == 0:
            self._set_status("No points selected.")
            return
        try:
            numeric_value = float(value)
            numeric_ok = True
        except ValueError:
            numeric_value = math.nan
            numeric_ok = False
        if numeric_ok and numeric_value == MISSING_VALUE:
            numeric_ok = False
            numeric_value = math.nan
        snapshot = self._make_undo_snapshot()
        rows_touched = set()
        for index_name, mask in self.selected_masks.items():
            if not mask.any():
                continue
            idx = self.block.index.index(index_name)
            rows_idx = np.where(mask)[0]
            for row_idx in rows_idx:
                row = self.rows[row_idx]
                self._ensure_value_length(row, idx)
                row.values[idx] = value
                rows_touched.add(int(row_idx))
            if index_name in self.series_cache:
                self.series_cache[index_name][rows_idx] = numeric_value if numeric_ok else math.nan
        self._save_undo_snapshot(snapshot, "value apply")
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Applied value to {total} points.")

    def _blank_value(self) -> None:
        total = int(sum(mask.sum() for mask in self.selected_masks.values()))
        if total == 0:
            self._set_status("No points selected.")
            return
        blank_value = f"{MISSING_VALUE:.1f}"
        snapshot = self._make_undo_snapshot()
        rows_touched = set()
        for index_name, mask in self.selected_masks.items():
            if not mask.any():
                continue
            idx = self.block.index.index(index_name)
            rows_idx = np.where(mask)[0]
            for row_idx in rows_idx:
                row = self.rows[row_idx]
                self._ensure_value_length(row, idx)
                row.values[idx] = blank_value
                rows_touched.add(int(row_idx))
            if index_name in self.series_cache:
                self.series_cache[index_name][rows_idx] = math.nan
        self._save_undo_snapshot(snapshot, "blank value")
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Blanked {total} points.")

    def _interpolate_missing_between(self, ys: np.ndarray) -> np.ndarray:
        out = ys.copy()
        finite = np.isfinite(ys)
        if np.count_nonzero(finite) < 2:
            return out
        x_valid = self.xs[finite]
        y_valid = ys[finite]
        order = np.argsort(x_valid)
        x_valid = x_valid[order]
        y_valid = y_valid[order]
        fill_mask = (~finite) & (self.xs >= x_valid[0]) & (self.xs <= x_valid[-1])
        if not fill_mask.any():
            return out
        out[fill_mask] = np.interp(self.xs[fill_mask], x_valid, y_valid)
        return out

    def _missing_entry_mask(self, idx: int) -> np.ndarray:
        mask = np.zeros(len(self.rows), dtype=bool)
        for row_idx, row in enumerate(self.rows):
            if idx >= len(row.values):
                mask[row_idx] = True
                continue
            raw = row.values[idx].strip()
            if raw == "":
                mask[row_idx] = True
                continue
            try:
                mask[row_idx] = float(raw) == MISSING_VALUE
            except ValueError:
                mask[row_idx] = False
        return mask

    def _interpolate_blanks(self) -> None:
        indices = self._selected_indices()
        if not indices:
            self._set_status("Select indices to interpolate.")
            return
        rows_mask = self._selected_row_mask(indices)
        if not rows_mask.any():
            rows_mask = np.ones(len(self.rows), dtype=bool)

        rows_touched = set()
        total_filled = 0
        snapshot = None
        for index_name in indices:
            idx = self.block.index.index(index_name)
            ys = self._build_series(index_name)
            interp = self._interpolate_missing_between(ys)
            missing_mask = self._missing_entry_mask(idx)
            fill_mask = rows_mask & missing_mask & np.isfinite(interp)
            if not fill_mask.any():
                continue
            rows_idx = np.where(fill_mask)[0]
            if snapshot is None:
                snapshot = self._make_undo_snapshot()
            for row_idx in rows_idx:
                row = self.rows[int(row_idx)]
                self._ensure_value_length(row, idx)
                row.values[idx] = f"{interp[row_idx]:.1f}"
                rows_touched.add(int(row_idx))
            if index_name in self.series_cache:
                self.series_cache[index_name][rows_idx] = interp[rows_idx]
            total_filled += int(rows_idx.size)

        if total_filled == 0:
            self._set_status("No blank points found between valid neighbors.")
            return
        self._save_undo_snapshot(snapshot, "blank interpolation")
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Interpolated {total_filled} blank points.")

    def _apply_smoothing(self) -> None:
        window_days = self._smooth_window_days()
        if window_days is None:
            self._set_status("Invalid smoothing window.")
            return
        indices = self._selected_indices()
        if not indices:
            self._set_status("Select indices to plot before applying smoothing.")
            return
        method = self._smooth_method_key()
        rows_touched = set()
        snapshot = None
        for index_name in indices:
            smooth = self._smooth_series(index_name, window_days, method)
            orig = self._build_series(index_name)
            mask = np.isfinite(orig) & np.isfinite(smooth)
            if not mask.any():
                continue
            idx = self.block.index.index(index_name)
            if snapshot is None:
                snapshot = self._make_undo_snapshot()
            for row_idx in np.where(mask)[0]:
                row = self.rows[row_idx]
                self._ensure_value_length(row, idx)
                row.values[idx] = f"{smooth[row_idx]:.1f}"
                rows_touched.add(int(row_idx))
            self.series_cache[index_name] = smooth.copy()
        if rows_touched:
            self._save_undo_snapshot(snapshot, f"smoothing ({method})")
            self._refresh_table_rows(rows_touched)
            self._update_plot()
            self._set_status(f"Applied smoothing ({method}) to {len(indices)} indices.")
        else:
            self._set_status("No points to smooth.")

    def _clear_selection(self) -> None:
        self._clear_selection_masks()
        self._update_selected_scatter()

    def _open_file(self) -> None:
        path = filedialog.askopenfilename(
            parent=self.root,
            title="Open ANTAB file",
            filetypes=[("ANTAB files", "*.antab"), ("All files", "*.*")],
        )
        if not path or not os.path.exists(path):
            return
        try:
            segments = parse_antab(path)
        except Exception as exc:
            self._set_status(f"Failed to open {path}: {exc}")
            return
        self.path = path
        self.segments = segments
        self.blocks = [seg.block for seg in segments if isinstance(seg, TsysSegment)]
        if not self.blocks:
            self._set_status(f"No TSYS blocks found in {path}.")
            return
        self.gain_info = _parse_gain_info(path)
        self.sefd_table = self._load_sefd_table()
        self.undo_snapshot = None
        self._update_undo_button_state()
        self.root.title(f"ANTAB GUI Editor — {os.path.basename(path)}")
        block_names = [f"{i}: {b.station}" for i, b in enumerate(self.blocks)]
        self.block_combo.configure(values=block_names)
        self._load_block(0)
        self._set_status(f"Opened {os.path.basename(path)}")

    def _open_generate_blank_dialog(self) -> None:
        GenerateBlankDialog(self.root, self.sefd_table, self.sefd_bands)

    def _save(self) -> None:
        if self.save_write_new.get():
            out_path = self._save_output_path()
            if os.path.abspath(out_path) == os.path.abspath(self.path):
                self._set_status("Output path matches current file. Choose a different path.")
                return
            write_antab(out_path, self.segments)
            self._set_status(f"Saved to {out_path}")
            return
        write_antab(self.path, self.segments)
        self._set_status(f"Saved to {self.path}")


def _read_fits_idi_params(path: str) -> Dict:
    """Read observation parameters from a single FITS-IDI file.

    If the file has a SYSTEM_TEMPERATURE table, antennas and time range come
    from it.  Otherwise (e.g. IDI2 with UV-only data) the active antennas are
    decoded from UV_DATA baselines and the time range from UV_DATA DATE+TIME.
    Times in the returned dict are fractional days relative to midnight of DATE-OBS.
    """
    try:
        from astropy.io import fits as _fits
    except ImportError:
        raise RuntimeError("astropy is required to read FITS files.\nInstall with: pip install astropy")

    import datetime as _dt

    with _fits.open(path, memmap=True) as hdul:
        ext_names = {h.name for h in hdul}

        uv_hdr = hdul["UV_DATA"].header
        date_obs: str = uv_hdr.get("DATE-OBS", "")
        ref_freq_hz: float = float(uv_hdr.get("REF_FREQ", 0.0))
        n_if: int = int(uv_hdr.get("NO_BAND", 1))

        freq_row = hdul["FREQUENCY"].data[0]
        band_offsets_hz = freq_row["BANDFREQ"]
        total_bw_hz = freq_row["TOTAL_BANDWIDTH"]
        if_freqs_hz = ref_freq_hz + band_offsets_hz
        center_freq_hz = float((if_freqs_hz[0] + if_freqs_hz[-1] + total_bw_hz[-1]) / 2)
        if_freqs_mhz = [float(f) / 1e6 for f in if_freqs_hz]

        ant_data = hdul["ANTENNA"].data
        ant_map: Dict[int, Tuple[str, str, str]] = {}
        for row in ant_data:
            ant_map[int(row["ANTENNA_NO"])] = (
                row["ANNAME"].strip(),
                row["POLTYA"].strip(),
                row["POLTYB"].strip(),
            )

        if "SYSTEM_TEMPERATURE" in ext_names:
            tsys_data = hdul["SYSTEM_TEMPERATURE"].data
            active_ant_nos = sorted(set(int(a) for a in tsys_data["ANTENNA_NO"]))
            time_start = float(tsys_data["TIME"].min())
            time_end   = float(tsys_data["TIME"].max())
        else:
            # Fall back to UV_DATA baselines and DATE+TIME columns.
            uv_data = hdul["UV_DATA"].data
            baselines = uv_data["BASELINE"]
            active_ant_nos = sorted(set(
                [int(b) // 256 for b in baselines] + [int(b) % 256 for b in baselines]
            ))
            # UV_DATA times: DATE is JD of midnight, TIME is fractional day offset.
            jd_midnight = (
                _dt.date.fromisoformat(date_obs).toordinal()
                - _dt.date(1858, 11, 17).toordinal()
                + 2400000.5
            ) if date_obs else 0.0
            abs_frac = (uv_data["DATE"] - jd_midnight) + uv_data["TIME"]
            time_start = float(abs_frac.min())
            time_end   = float(abs_frac.max())

        stations = [ant_map[a][0] for a in active_ant_nos if a in ant_map]
        first_ant = active_ant_nos[0] if active_ant_nos else None
        pol_a = ant_map[first_ant][1] if first_ant and first_ant in ant_map else "R"
        pol_b = ant_map[first_ant][2] if first_ant and first_ant in ant_map else "L"

    base_doy = 0
    if date_obs:
        try:
            base_doy = _dt.datetime.strptime(date_obs, "%Y-%m-%d").timetuple().tm_yday
        except ValueError:
            pass

    return {
        "date_obs": date_obs,
        "base_doy": base_doy,
        "n_if": n_if,
        "center_freq_hz": center_freq_hz,
        "if_freqs_mhz": if_freqs_mhz,
        "stations": stations,
        "pol_a": pol_a,
        "pol_b": pol_b,
        "time_start_frac": time_start,
        "time_end_frac": time_end,
    }


def _merge_fits_idi_params(paths: List[str]) -> Dict:
    """Read and merge parameters from multiple FITS-IDI files.

    Stations are the union across all files.  Frequency/IF/polarization come
    from the first file that provides them.  Time range is the global min/max.
    """
    if not paths:
        raise ValueError("No FITS files supplied.")
    merged: Optional[Dict] = None
    errors: List[str] = []
    for p in paths:
        try:
            params = _read_fits_idi_params(p)
        except Exception as exc:
            errors.append(f"{os.path.basename(p)}: {exc}")
            continue
        if merged is None:
            merged = dict(params)
            merged["stations"] = list(params["stations"])
        else:
            # Union of stations (preserve order, no duplicates)
            existing = {s.upper() for s in merged["stations"]}
            for s in params["stations"]:
                if s.upper() not in existing:
                    merged["stations"].append(s)
                    existing.add(s.upper())
            # Expand time range
            merged["time_start_frac"] = min(merged["time_start_frac"], params["time_start_frac"])
            merged["time_end_frac"]   = max(merged["time_end_frac"],   params["time_end_frac"])
            # Use first file's freq/IF/pol unless missing
            if merged["center_freq_hz"] == 0.0 and params["center_freq_hz"] != 0.0:
                merged["center_freq_hz"] = params["center_freq_hz"]
                merged["if_freqs_mhz"]   = params["if_freqs_mhz"]
            if merged["n_if"] == 0 and params["n_if"] > 0:
                merged["n_if"] = params["n_if"]
    if merged is None:
        raise RuntimeError("Could not read any FITS file:\n" + "\n".join(errors))
    merged["load_errors"] = errors
    return merged


# ---------------------------------------------------------------------------
# Batch / non-interactive operations
# These are standalone functions (not methods) so they work without the GUI.
# ---------------------------------------------------------------------------

def _block_rows_xs(block: TsysBlock):
    """Return (rows, xs) from a TsysBlock without needing an AntabGui instance."""
    rows: List[DataRow] = []
    xs_list: List[float] = []
    for row in block.data_rows:
        if isinstance(row, DataRow):
            rows.append(row)
            xs_list.append(parse_time(row.day, row.time))
    return rows, np.asarray(xs_list, dtype=float)


def _build_series_batch(rows: List[DataRow], block: TsysBlock, index_name: str) -> np.ndarray:
    idx = block.index.index(index_name)
    ys = np.full(len(rows), np.nan, dtype=float)
    for i, row in enumerate(rows):
        if idx >= len(row.values):
            continue
        val = row.values[idx].strip()
        if not val:
            continue
        try:
            fval = float(val)
            if fval != MISSING_VALUE:
                ys[i] = fval
        except ValueError:
            pass
    return ys


def _smooth_series_batch(xs: np.ndarray, ys: np.ndarray, window_days: float, method: str) -> np.ndarray:
    if window_days <= 0:
        return np.full_like(ys, np.nan)
    mask = np.isfinite(ys)
    if not mask.any():
        return np.full_like(ys, np.nan)
    x_valid = xs[mask]
    y_valid = ys[mask]
    n = x_valid.size
    half = window_days / 2.0
    smoothed_valid = np.full(n, np.nan, dtype=float)
    if method == "mean":
        for i in range(n):
            xi = x_valid[i]
            l = int(np.searchsorted(x_valid, xi - half, side="left"))
            r = int(np.searchsorted(x_valid, xi + half, side="right"))
            if r > l:
                smoothed_valid[i] = float(np.mean(y_valid[l:r]))
    elif method == "median":
        for i in range(n):
            xi = x_valid[i]
            l = int(np.searchsorted(x_valid, xi - half, side="left"))
            r = int(np.searchsorted(x_valid, xi + half, side="right"))
            if r > l:
                smoothed_valid[i] = float(np.median(y_valid[l:r]))
    else:  # gaussian
        sigma = window_days / 2.355
        if sigma <= 0:
            return np.full_like(ys, np.nan)
        for i in range(n):
            xi = x_valid[i]
            l = int(np.searchsorted(x_valid, xi - half, side="left"))
            r = int(np.searchsorted(x_valid, xi + half, side="right"))
            if r <= l:
                continue
            dt = x_valid[l:r] - xi
            w = np.exp(-0.5 * (dt / sigma) ** 2)
            wsum = float(np.sum(w))
            if wsum > 0:
                smoothed_valid[i] = float(np.sum(w * y_valid[l:r]) / wsum)
    smoothed = np.full_like(ys, np.nan)
    smoothed[mask] = smoothed_valid
    return smoothed


def _interpolate_missing_batch(xs: np.ndarray, ys: np.ndarray) -> np.ndarray:
    out = ys.copy()
    finite = np.isfinite(ys)
    if np.count_nonzero(finite) < 2:
        return out
    x_valid = xs[finite]
    y_valid = ys[finite]
    order = np.argsort(x_valid)
    x_valid = x_valid[order]
    y_valid = y_valid[order]
    fill_mask = (~finite) & (xs >= x_valid[0]) & (xs <= x_valid[-1])
    if fill_mask.any():
        out[fill_mask] = np.interp(xs[fill_mask], x_valid, y_valid)
    return out


def _missing_entry_mask_batch(rows: List[DataRow], idx: int) -> np.ndarray:
    mask = np.zeros(len(rows), dtype=bool)
    for i, row in enumerate(rows):
        if idx >= len(row.values):
            mask[i] = True
            continue
        raw = row.values[idx].strip()
        if not raw:
            mask[i] = True
            continue
        try:
            mask[i] = float(raw) == MISSING_VALUE
        except ValueError:
            pass
    return mask


def _set_gain_line_batch(lines: List[str], start: int, end: int, key: str, values: List[float]) -> Tuple[int, bool]:
    text = _format_float_list(values)
    for i in range(start, end):
        raw = lines[i]
        head = raw.strip().split("=", 1)[0].strip().upper() if raw.strip() else ""
        if head == key:
            indent = raw[: len(raw) - len(raw.lstrip())]
            new_line = f"{indent}{key} = {text}"
            if lines[i] != new_line:
                lines[i] = new_line
                return end, True
            return end, False
    lines.insert(end, f"{key} = {text}")
    return end + 1, True


def _set_dpfu_to_one(segments: List, station: str, gain_info: Dict) -> None:
    station = _normalize_station_key(station)
    dpfu_count = max(1, len(gain_info.get(station, {}).get("dpfu", [])))
    dpfu_vals = [1.0] * dpfu_count
    for seg in segments:
        if not isinstance(seg, RawSegment):
            continue
        lines = seg.lines
        i = 0
        while i < len(lines):
            parts = lines[i].strip().split()
            if len(parts) >= 2 and parts[0].upper() == "GAIN" and _normalize_station_key(parts[1]) == station:
                end = i + 1
                while end < len(lines) and lines[end].strip() != "/":
                    end += 1
                _set_gain_line_batch(lines, i + 1, end, "DPFU", dpfu_vals)
                i = end
            i += 1
    gain_info.setdefault(station, {"dpfu": [], "freq": [], "poly": []})["dpfu"] = dpfu_vals


def _expected_tsys_batch(block: TsysBlock, index_name: str, gain_info: Dict, sefd_table: Dict) -> Optional[float]:
    station_key = _normalize_station_key(block.station)
    gain = gain_info.get(station_key)
    if not gain:
        return None
    dpfu_vals = gain.get("dpfu", [])
    freq_vals = gain.get("freq", [])
    if not freq_vals:
        return None
    freq_hz = (sum(freq_vals) / len(freq_vals)) * 1e6
    sefd_map = sefd_table.get(station_key)
    if not sefd_map:
        return None
    nearest = min(sefd_map.keys(), key=lambda f: abs(f - freq_hz))
    if not dpfu_vals:
        return None
    upper = index_name.upper()
    if len(dpfu_vals) == 1:
        dpfu = dpfu_vals[0]
    elif "R" in upper and "L" not in upper:
        dpfu = dpfu_vals[0]
    elif "L" in upper and "R" not in upper:
        dpfu = dpfu_vals[1] if len(dpfu_vals) > 1 else dpfu_vals[0]
    else:
        dpfu = dpfu_vals[0]
    return dpfu * sefd_map[nearest]


def op_expected_tsys(segments: List, gain_info: Dict, sefd_table: Dict, station_filter: Optional[str] = None) -> None:
    for seg in segments:
        if not isinstance(seg, TsysSegment):
            continue
        block = seg.block
        if station_filter and _normalize_station_key(block.station) != station_filter:
            continue
        rows, _ = _block_rows_xs(block)
        applied = 0
        for index_name in block.index:
            expected = _expected_tsys_batch(block, index_name, gain_info, sefd_table)
            if expected is None:
                continue
            idx = block.index.index(index_name)
            for row in rows:
                while len(row.values) <= idx:
                    row.values.append("")
                row.values[idx] = f"{expected:.1f}"
            applied += 1
        if applied:
            _set_dpfu_to_one(segments, block.station, gain_info)
            print(f"  {block.station}: applied expected TSYS to {applied} indices, set DPFU=1")
        else:
            print(f"  {block.station}: no SEFD/GAIN data — skipped", file=sys.stderr)


def op_interpolate(segments: List, station_filter: Optional[str] = None) -> None:
    for seg in segments:
        if not isinstance(seg, TsysSegment):
            continue
        block = seg.block
        if station_filter and _normalize_station_key(block.station) != station_filter:
            continue
        rows, xs = _block_rows_xs(block)
        total = 0
        for index_name in block.index:
            idx = block.index.index(index_name)
            ys = _build_series_batch(rows, block, index_name)
            interp = _interpolate_missing_batch(xs, ys)
            fill_mask = _missing_entry_mask_batch(rows, idx) & np.isfinite(interp)
            for row_idx in np.where(fill_mask)[0]:
                row = rows[int(row_idx)]
                while len(row.values) <= idx:
                    row.values.append("")
                row.values[idx] = f"{interp[row_idx]:.1f}"
            total += int(fill_mask.sum())
        print(f"  {block.station}: interpolated {total} blank points")


def op_smooth(segments: List, method: str, window_minutes: float, station_filter: Optional[str] = None) -> None:
    window_days = window_minutes / (24 * 60)
    for seg in segments:
        if not isinstance(seg, TsysSegment):
            continue
        block = seg.block
        if station_filter and _normalize_station_key(block.station) != station_filter:
            continue
        rows, xs = _block_rows_xs(block)
        total = 0
        for index_name in block.index:
            idx = block.index.index(index_name)
            ys = _build_series_batch(rows, block, index_name)
            smooth = _smooth_series_batch(xs, ys, window_days, method)
            mask = np.isfinite(ys) & np.isfinite(smooth)
            for row_idx in np.where(mask)[0]:
                row = rows[int(row_idx)]
                while len(row.values) <= idx:
                    row.values.append("")
                row.values[idx] = f"{smooth[row_idx]:.1f}"
            total += int(mask.sum())
        print(f"  {block.station}: smoothed {total} values ({method}, {window_minutes:.1f} min)")


def _seconds_to_time_str_fn(seconds: int) -> str:
    h, m, s = seconds // 3600, (seconds % 3600) // 60, seconds % 60
    return f"{h:02d}:{m:02d}:{s:02d}"


def _gen_time_rows(
    start_doy: int, start_secs: int,
    end_doy: int, end_secs: int,
    interval_secs: int,
    max_rows: int = 50000,
) -> List[Tuple[str, str]]:
    cur = start_doy * 86400 + start_secs
    end = end_doy * 86400 + end_secs
    if cur > end:
        raise ValueError("Start time is after end time.")
    rows: List[Tuple[str, str]] = []
    while cur <= end:
        rows.append((str(cur // 86400).zfill(3), _seconds_to_time_str_fn(cur % 86400)))
        cur += interval_secs
        if len(rows) > max_rows:
            raise ValueError(f"Time range exceeds {max_rows} rows. Reduce range or increase --interval.")
    if not rows:
        raise ValueError("No time rows generated. Check time range and interval.")
    return rows


def generate_blank(
    fits_paths: List[str],
    out_path: str,
    sefd_table: Dict,
    sefd_bands: List,
    stations_override: Optional[List[str]] = None,
    pols_override: Optional[str] = None,
    nif_override: Optional[int] = None,
    if_start: int = 1,
    band_cm_override: Optional[float] = None,
    interval_min: float = 1.0,
) -> None:
    params = _merge_fits_idi_params(fits_paths)
    if params.get("load_errors"):
        for e in params["load_errors"]:
            print(f"  Warning: {e}", file=sys.stderr)

    if stations_override:
        stations = [_normalize_station_key(s) for s in stations_override]
    else:
        all_stations = [_normalize_station_key(s) for s in params["stations"]]
        in_sefd = [s for s in all_stations if s in sefd_table]
        stations = in_sefd if in_sefd else all_stations

    if band_cm_override is not None:
        freq_hz = _band_cm_to_hz(band_cm_override)
    elif sefd_bands and params["center_freq_hz"]:
        _, freq_hz = min(sefd_bands, key=lambda b: abs(b[1] - params["center_freq_hz"]))
    else:
        freq_hz = params["center_freq_hz"]
    freq_mhz = freq_hz / 1e6

    nif = nif_override if nif_override is not None else params["n_if"]
    if pols_override is not None:
        pols = list(pols_override.upper())
    else:
        pols = []
        for p in (params["pol_a"], params["pol_b"]):
            if p in ("R", "L") and p not in pols:
                pols.append(p)
        if not pols:
            pols = ["R", "L"]
    indices = [f"{pol}{i}" for pol in pols for i in range(if_start, if_start + nif)]

    base_doy = params["base_doy"]
    t0 = params["time_start_frac"]
    t1 = params["time_end_frac"]
    d0, s0 = base_doy + int(t0), int((t0 % 1) * 86400)
    d1, s1 = base_doy + int(t1), int((t1 % 1) * 86400)
    time_rows = _gen_time_rows(d0, s0, d1, s1, max(1, int(round(interval_min * 60))))

    index_str = ", ".join(f"'{name}'" for name in indices)
    lines: List[str] = []
    missing_sefd: List[str] = []

    for station in stations:
        sefd_map = sefd_table.get(station, {})
        if sefd_map:
            nearest_hz = min(sefd_map.keys(), key=lambda f: abs(f - freq_hz))
            sefd_jy: Optional[float] = sefd_map[nearest_hz]
        else:
            sefd_jy = None
            missing_sefd.append(station)
        val_str = (
            " ".join(f"{sefd_jy:.1f}" for _ in indices)
            if sefd_jy is not None
            else " ".join(f"{MISSING_VALUE:.1f}" for _ in indices)
        )
        lines += [
            f"GAIN {station}", "ELEV", "DPFU = 1", "POLY = 1", f"FREQ = {freq_mhz:g}", "/",
            f"TSYS {station}", "FT = 1", "TIMEOFF = 0", f"INDEX = {index_str}", "/",
        ]
        for doy_str, time_str in time_rows:
            lines.append(f"{doy_str} {time_str} {val_str}")
        lines.append("/")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    print(
        f"Generated {out_path}: {len(stations)} station(s), {len(indices)} indices, "
        f"{len(time_rows)} time rows, {freq_mhz:.1f} MHz."
    )
    if missing_sefd:
        print(f"  No SEFD data for: {', '.join(missing_sefd)} — used {MISSING_VALUE:.1f}.", file=sys.stderr)


class GenerateBlankDialog:
    """Dialog for generating a blank ANTAB file with estimated SEFD values.

    Parameters come from a FITS-IDI file (stations, frequency, IFs,
    polarizations, time range). TSYS values are filled with SEFD estimates
    from sefd_values.txt — no TSYS data is read from the FITS file.
    """

    MAX_ROWS = 50000

    def __init__(
        self,
        parent: tk.Tk,
        sefd_table: Dict[str, Dict[float, float]],
        sefd_bands: List[Tuple[float, float]],
    ) -> None:
        self.parent = parent
        self.sefd_table = sefd_table
        self.sefd_bands = sefd_bands  # [(cm, hz), ...]

        self.win = tk.Toplevel(parent)
        self.win.title("Generate Blank ANTAB")
        self.win.resizable(True, True)
        self.win.grab_set()
        self._build()

    def _build(self) -> None:
        outer = ttk.Frame(self.win, padding=10)
        outer.pack(fill=tk.BOTH, expand=True)

        # ---- FITS loading ----
        fits_frame = ttk.LabelFrame(outer, text="Load parameters from FITS-IDI files")
        fits_frame.grid(row=0, column=0, columnspan=3, sticky="ew", padx=4, pady=(0, 8))

        fits_inner = ttk.Frame(fits_frame)
        fits_inner.pack(fill=tk.BOTH, padx=6, pady=6)

        # File list
        list_col = ttk.Frame(fits_inner)
        list_col.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        self.fits_listbox = tk.Listbox(list_col, selectmode=tk.EXTENDED, exportselection=False,
                                       height=4, width=60)
        fits_sb = ttk.Scrollbar(list_col, orient=tk.VERTICAL, command=self.fits_listbox.yview)
        self.fits_listbox.configure(yscrollcommand=fits_sb.set)
        self.fits_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        fits_sb.pack(side=tk.LEFT, fill=tk.Y)

        # Buttons beside list
        btn_col = ttk.Frame(fits_inner)
        btn_col.pack(side=tk.LEFT, padx=(6, 0), anchor="n")
        ttk.Button(btn_col, text="Add files…", command=self._browse_fits).pack(fill=tk.X, pady=2)
        ttk.Button(btn_col, text="Remove selected", command=self._remove_fits).pack(fill=tk.X, pady=2)
        ttk.Button(btn_col, text="Load from FITS", command=self._load_from_fits).pack(fill=tk.X, pady=(10, 2))

        # Auto-detect IDI files in cwd or data/ directory
        self._auto_add_fits_files()

        frame = ttk.Frame(outer)
        frame.grid(row=1, column=0, columnspan=3, sticky="nsew")
        outer.rowconfigure(1, weight=1)
        outer.columnconfigure(0, weight=1)
        outer.columnconfigure(1, weight=1)
        outer.columnconfigure(2, weight=1)

        # ---- Stations ----
        sta_frame = ttk.LabelFrame(frame, text="Stations")
        sta_frame.grid(row=0, column=0, rowspan=2, sticky="nsew", padx=4, pady=4)

        self.sta_listbox = tk.Listbox(sta_frame, selectmode=tk.EXTENDED, exportselection=False, width=10, height=16)
        sta_sb = ttk.Scrollbar(sta_frame, orient=tk.VERTICAL, command=self.sta_listbox.yview)
        self.sta_listbox.configure(yscrollcommand=sta_sb.set)
        self.sta_listbox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(4, 0), pady=4)
        sta_sb.pack(side=tk.RIGHT, fill=tk.Y, pady=4, padx=(0, 4))

        self._all_stations = sorted(self.sefd_table.keys())
        for sta in self._all_stations:
            self.sta_listbox.insert(tk.END, sta)

        ttk.Button(sta_frame, text="Select All", command=lambda: self.sta_listbox.selection_set(0, tk.END)).pack(fill=tk.X, padx=4, pady=(0, 4))

        # ---- Band / Frequency ----
        bf_frame = ttk.LabelFrame(frame, text="Band / Frequency")
        bf_frame.grid(row=0, column=1, sticky="new", padx=4, pady=4)

        ttk.Label(bf_frame, text="Band (cm)").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.band_var = tk.StringVar()
        band_names = [f"{cm:g}" for cm, _ in self.sefd_bands]
        self.band_combo = ttk.Combobox(bf_frame, textvariable=self.band_var, values=band_names, state="readonly", width=12)
        if band_names:
            self.band_combo.current(0)
        self.band_combo.pack(fill=tk.X, padx=6)
        self.band_combo.bind("<<ComboboxSelected>>", self._on_band_change)

        ttk.Label(bf_frame, text="Center frequency (MHz)").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.freq_var = tk.StringVar()
        ttk.Entry(bf_frame, textvariable=self.freq_var, width=14).pack(fill=tk.X, padx=6, pady=(0, 6))

        if self.sefd_bands:
            self._update_freq_from_hz(self.sefd_bands[0][1])

        # ---- Index Pattern ----
        idx_frame = ttk.LabelFrame(frame, text="Index Pattern")
        idx_frame.grid(row=1, column=1, sticky="new", padx=4, pady=4)

        self.pol_l_var = tk.BooleanVar(value=True)
        self.pol_r_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(idx_frame, text="L polarization", variable=self.pol_l_var).pack(anchor=tk.W, padx=6, pady=2)
        ttk.Checkbutton(idx_frame, text="R polarization", variable=self.pol_r_var).pack(anchor=tk.W, padx=6, pady=2)

        ttk.Label(idx_frame, text="Number of IFs").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.nif_var = tk.StringVar(value="8")
        ttk.Entry(idx_frame, textvariable=self.nif_var, width=8).pack(anchor=tk.W, padx=6)

        ttk.Label(idx_frame, text="Starting IF number").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.if_start_var = tk.StringVar(value="1")
        ttk.Entry(idx_frame, textvariable=self.if_start_var, width=8).pack(anchor=tk.W, padx=6, pady=(0, 6))

        # ---- Time Range ----
        time_frame = ttk.LabelFrame(frame, text="Time Range")
        time_frame.grid(row=0, column=2, sticky="new", padx=4, pady=4)

        for label, attr, default in [
            ("Start DOY", "start_doy_var", "001"),
            ("Start time (HH:MM:SS)", "start_time_var", "00:00:00"),
            ("End DOY", "end_doy_var", "001"),
            ("End time (HH:MM:SS)", "end_time_var", "23:59:00"),
            ("Interval (minutes)", "interval_var", "5"),
        ]:
            ttk.Label(time_frame, text=label).pack(anchor=tk.W, padx=6, pady=(6, 2))
            var = tk.StringVar(value=default)
            setattr(self, attr, var)
            ttk.Entry(time_frame, textvariable=var, width=14).pack(fill=tk.X, padx=6)

        # ---- Output ----
        out_frame = ttk.LabelFrame(frame, text="Output")
        out_frame.grid(row=1, column=2, sticky="new", padx=4, pady=4)

        ttk.Label(out_frame, text="Output path").pack(anchor=tk.W, padx=6, pady=(6, 2))
        self.out_var = tk.StringVar(value="blank.antab")
        out_row = ttk.Frame(out_frame)
        out_row.pack(fill=tk.X, padx=6)
        ttk.Entry(out_row, textvariable=self.out_var).pack(side=tk.LEFT, fill=tk.X, expand=True)
        ttk.Button(out_row, text="...", width=3, command=self._browse_output).pack(side=tk.LEFT, padx=(4, 0))

        # ---- Status / Buttons ----
        self.dlg_status_var = tk.StringVar(value="")
        ttk.Label(outer, textvariable=self.dlg_status_var, wraplength=580).grid(
            row=2, column=0, columnspan=3, sticky="w", padx=4, pady=(6, 0)
        )

        btn_frame = ttk.Frame(outer)
        btn_frame.grid(row=3, column=0, columnspan=3, sticky="e", padx=4, pady=6)
        ttk.Button(btn_frame, text="Close", command=self.win.destroy).pack(side=tk.RIGHT, padx=(4, 0))
        ttk.Button(btn_frame, text="Generate", command=self._generate).pack(side=tk.RIGHT)

        frame.columnconfigure(0, weight=1)
        frame.columnconfigure(1, weight=1)
        frame.columnconfigure(2, weight=1)

    # ---- FITS loading ----

    @staticmethod
    def _is_fits_idi(name: str) -> bool:
        upper = name.upper()
        return (
            upper.endswith((".IDI", ".FITS"))
            or upper.endswith((".IDI1", ".IDI2", ".IDI3", ".IDI4"))
            or ("IDI" in upper and not upper.endswith(".PY"))
        )

    def _auto_add_fits_files(self) -> None:
        found: List[str] = []
        for name in sorted(os.listdir(".")):
            if os.path.isfile(name) and self._is_fits_idi(name):
                found.append(name)
        if not found:
            data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
            if os.path.isdir(data_dir):
                for name in sorted(os.listdir(data_dir)):
                    full = os.path.join(data_dir, name)
                    if os.path.isfile(full) and self._is_fits_idi(name):
                        found.append(full)
        for p in found:
            self.fits_listbox.insert(tk.END, p)

    def _browse_fits(self) -> None:
        paths = filedialog.askopenfilenames(
            parent=self.win,
            title="Select FITS-IDI files",
            filetypes=[
                ("FITS-IDI files", "*.IDI *.IDI1 *.IDI2 *.IDI3 *.IDI4 *.fits *.FITS"),
                ("All files", "*.*"),
            ],
        )
        existing = set(self.fits_listbox.get(0, tk.END))
        for p in paths:
            if p not in existing:
                self.fits_listbox.insert(tk.END, p)
                existing.add(p)

    def _remove_fits(self) -> None:
        for i in reversed(self.fits_listbox.curselection()):
            self.fits_listbox.delete(i)

    def _load_from_fits(self) -> None:
        paths = list(self.fits_listbox.get(0, tk.END))
        if not paths:
            self.dlg_status_var.set("Add at least one FITS-IDI file to the list first.")
            return
        missing = [p for p in paths if not os.path.exists(p)]
        if missing:
            self.dlg_status_var.set(f"File(s) not found: {', '.join(missing)}")
            return
        try:
            params = _merge_fits_idi_params(paths)
        except Exception as exc:
            self.dlg_status_var.set(f"Error reading FITS: {exc}")
            return

        # --- Stations ---
        fits_stations_upper = {s.upper() for s in params["stations"]}
        self.sta_listbox.selection_clear(0, tk.END)
        for i, sta in enumerate(self._all_stations):
            if sta.upper() in fits_stations_upper:
                self.sta_listbox.selection_set(i)

        # --- Frequency: pick nearest band in SEFD table ---
        center_hz = params["center_freq_hz"]
        if self.sefd_bands:
            nearest_cm, nearest_hz = min(self.sefd_bands, key=lambda bnd: abs(bnd[1] - center_hz))
            band_names = [f"{cm:g}" for cm, _ in self.sefd_bands]
            key = f"{nearest_cm:g}"
            if key in band_names:
                self.band_combo.current(band_names.index(key))
        self.freq_var.set(f"{center_hz / 1e6:.2f}")

        # --- Index pattern ---
        self.nif_var.set(str(params["n_if"]))
        self.if_start_var.set("1")
        self.pol_l_var.set("L" in (params["pol_a"], params["pol_b"]))
        self.pol_r_var.set("R" in (params["pol_a"], params["pol_b"]))

        # --- Time range ---
        base_doy = params["base_doy"]
        t0 = params["time_start_frac"]
        t1 = params["time_end_frac"]
        d0 = base_doy + int(t0);  s0 = int((t0 % 1) * 86400)
        d1 = base_doy + int(t1);  s1 = int((t1 % 1) * 86400)
        self.start_doy_var.set(str(d0).zfill(3))
        self.start_time_var.set(self._seconds_to_time_str(s0))
        self.end_doy_var.set(str(d1).zfill(3))
        self.end_time_var.set(self._seconds_to_time_str(s1))

        n_files = len(paths)
        n_fits = len(params["stations"])
        n_sefd = sum(1 for s in params["stations"] if s.upper() in {k.upper() for k in self.sefd_table})
        msg = (
            f"Merged {n_files} file(s): {n_fits} station(s) ({n_sefd} in SEFD table), "
            f"{params['n_if']} IFs, center freq {center_hz / 1e6:.1f} MHz, "
            f"DATE-OBS {params['date_obs']}."
        )
        if params.get("load_errors"):
            msg += "  Warnings: " + "; ".join(params["load_errors"])
        self.dlg_status_var.set(msg)

    # ---- helpers ----

    def _on_band_change(self, event=None) -> None:
        raw = self.band_var.get().strip()
        try:
            cm = float(raw)
        except ValueError:
            return
        self._update_freq_from_hz(_band_cm_to_hz(cm))

    def _update_freq_from_hz(self, hz: float) -> None:
        self.freq_var.set(f"{hz / 1e6:.2f}")

    def _browse_output(self) -> None:
        path = filedialog.asksaveasfilename(
            parent=self.win,
            title="Save blank ANTAB as",
            defaultextension=".antab",
            filetypes=[("ANTAB files", "*.antab"), ("All files", "*.*")],
        )
        if path:
            self.out_var.set(path)

    def _get_selected_stations(self) -> List[str]:
        return [self.sta_listbox.get(i) for i in self.sta_listbox.curselection()]

    def _build_index_list(self) -> List[str]:
        try:
            nif = int(self.nif_var.get().strip())
            if_start = int(self.if_start_var.get().strip())
        except ValueError:
            return []
        if nif <= 0:
            return []
        pols = []
        if self.pol_l_var.get():
            pols.append("L")
        if self.pol_r_var.get():
            pols.append("R")
        return [f"{pol}{i}" for pol in pols for i in range(if_start, if_start + nif)]

    def _parse_time_seconds(self, time_str: str) -> Optional[int]:
        parts = time_str.strip().split(":")
        if len(parts) != 3:
            return None
        try:
            h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
            return h * 3600 + m * 60 + s
        except ValueError:
            return None

    @staticmethod
    def _seconds_to_time_str(seconds: int) -> str:
        h = seconds // 3600
        m = (seconds % 3600) // 60
        s = seconds % 60
        return f"{h:02d}:{m:02d}:{s:02d}"

    def _generate_time_rows(self) -> Optional[List[Tuple[str, str]]]:
        try:
            start_doy = int(self.start_doy_var.get().strip())
            end_doy = int(self.end_doy_var.get().strip())
        except ValueError:
            self.dlg_status_var.set("Invalid DOY value.")
            return None
        start_secs = self._parse_time_seconds(self.start_time_var.get())
        end_secs = self._parse_time_seconds(self.end_time_var.get())
        if start_secs is None:
            self.dlg_status_var.set("Invalid start time. Use HH:MM:SS.")
            return None
        if end_secs is None:
            self.dlg_status_var.set("Invalid end time. Use HH:MM:SS.")
            return None
        try:
            interval_min = float(self.interval_var.get().strip())
            if interval_min <= 0:
                raise ValueError
        except ValueError:
            self.dlg_status_var.set("Interval must be a positive number.")
            return None
        interval_secs = max(1, int(round(interval_min * 60)))
        cur_total = start_doy * 86400 + start_secs
        end_total = end_doy * 86400 + end_secs
        if cur_total > end_total:
            self.dlg_status_var.set("Start time is after end time.")
            return None
        rows: List[Tuple[str, str]] = []
        while cur_total <= end_total:
            doy = cur_total // 86400
            secs = cur_total % 86400
            rows.append((str(doy).zfill(3), self._seconds_to_time_str(secs)))
            cur_total += interval_secs
            if len(rows) > self.MAX_ROWS:
                self.dlg_status_var.set(
                    f"Time range exceeds {self.MAX_ROWS} rows. Reduce the range or increase the interval."
                )
                return None
        if not rows:
            self.dlg_status_var.set("No rows generated. Check time range and interval.")
            return None
        return rows

    def _generate(self) -> None:
        stations = self._get_selected_stations()
        if not stations:
            self.dlg_status_var.set("Select at least one station.")
            return
        indices = self._build_index_list()
        if not indices:
            self.dlg_status_var.set("Select at least one polarization and set number of IFs > 0.")
            return
        try:
            freq_mhz = float(self.freq_var.get().strip())
        except ValueError:
            self.dlg_status_var.set("Invalid frequency.")
            return
        freq_hz = freq_mhz * 1e6
        time_rows = self._generate_time_rows()
        if time_rows is None:
            return
        out_path = self.out_var.get().strip()
        if not out_path:
            self.dlg_status_var.set("Enter output path.")
            return

        index_str = ", ".join(f"'{name}'" for name in indices)
        lines: List[str] = []
        missing_sefd: List[str] = []

        for station in stations:
            sefd_map = self.sefd_table.get(station, {})
            if sefd_map:
                nearest_hz = min(sefd_map.keys(), key=lambda f: abs(f - freq_hz))
                sefd_jy: Optional[float] = sefd_map[nearest_hz]
            else:
                sefd_jy = None
                missing_sefd.append(station)

            lines += [
                f"GAIN {station}",
                "ELEV",
                "DPFU = 1",
                "POLY = 1",
                f"FREQ = {freq_mhz:g}",
                "/",
                f"TSYS {station}",
                "FT = 1",
                "TIMEOFF = 0",
                f"INDEX = {index_str}",
                "/",
            ]
            val_str = (
                " ".join(f"{sefd_jy:.1f}" for _ in indices)
                if sefd_jy is not None
                else " ".join(f"{MISSING_VALUE:.1f}" for _ in indices)
            )
            for doy_str, time_str in time_rows:
                lines.append(f"{doy_str} {time_str} {val_str}")
            lines.append("/")

        try:
            with open(out_path, "w", encoding="utf-8") as f:
                f.write("\n".join(lines) + "\n")
        except OSError as exc:
            self.dlg_status_var.set(f"Error writing file: {exc}")
            return

        msg = (
            f"Saved {out_path}: {len(stations)} station(s), {len(indices)} indices, "
            f"{len(time_rows)} time rows."
        )
        if missing_sefd:
            msg += f" No SEFD for: {', '.join(missing_sefd)} (used {MISSING_VALUE:.1f})."
        self.dlg_status_var.set(msg)


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(
        description="ANTAB editor — GUI by default, or non-interactive with batch flags.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Non-interactive examples (no display required):
  antab_editor.py input.antab --expected-tsys -o output.antab
  antab_editor.py input.antab --smooth median 5 -o output.antab
  antab_editor.py input.antab --expected-tsys --interpolate --smooth gaussian 5 -o output.antab
  antab_editor.py --from-fits scan*.fits --interval 1 --expected-tsys -o output.antab
""",
    )
    parser.add_argument("path", nargs="?", help=".antab file path")
    # ---- FITS generation ----
    parser.add_argument("--from-fits", nargs="+", metavar="FILE",
                        help="Generate a blank ANTAB from one or more FITS-IDI files (non-interactive)")
    parser.add_argument("--stations", metavar="STA,STA,...",
                        help="Comma-separated station list for --from-fits (default: all in FITS)")
    parser.add_argument("--pols", metavar="POLS",
                        help="Polarizations for --from-fits, e.g. RL, R, L (default: from FITS)")
    parser.add_argument("--nif", type=int, metavar="N",
                        help="Number of IFs for --from-fits (default: from FITS)")
    parser.add_argument("--if-start", type=int, default=1, metavar="N",
                        help="First IF number for --from-fits (default: 1)")
    parser.add_argument("--band", type=float, metavar="CM",
                        help="Wavelength in cm for frequency lookup in --from-fits")
    parser.add_argument("--interval", type=float, default=1.0, metavar="MIN",
                        help="Time interval in minutes for --from-fits (default: 1)")
    # ---- batch operations ----
    parser.add_argument("--expected-tsys", action="store_true",
                        help="Replace TSYS values with SEFD * DPFU; sets DPFU=1 (non-interactive)")
    parser.add_argument("--smooth", nargs=2, metavar=("METHOD", "MIN"),
                        help="Smooth data: METHOD is mean/median/gaussian, MIN is window in minutes")
    parser.add_argument("--interpolate", action="store_true",
                        help="Interpolate blank (-99.0) values between valid neighbours")
    parser.add_argument("--station", metavar="STA",
                        help="Restrict batch operations to this station only")
    parser.add_argument("-o", "--output", metavar="FILE",
                        help="Output path for batch mode (default: overwrite input)")
    args = parser.parse_args(argv)

    batch_mode = bool(args.from_fits or args.expected_tsys or args.smooth or args.interpolate)

    # ------------------------------------------------------------------ #
    #  Batch / non-interactive path                                        #
    # ------------------------------------------------------------------ #
    if batch_mode:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        sefd_table, sefd_bands = _load_sefd_table_from_dirs([script_dir, os.getcwd()])
        antab_path = args.path
        out_path = args.output

        if args.from_fits:
            generated = out_path or "blank.antab"
            stations_override = [s.strip() for s in args.stations.split(",")] if args.stations else None
            try:
                generate_blank(
                    fits_paths=args.from_fits,
                    out_path=generated,
                    sefd_table=sefd_table,
                    sefd_bands=sefd_bands,
                    stations_override=stations_override,
                    pols_override=args.pols,
                    nif_override=args.nif,
                    if_start=args.if_start,
                    band_cm_override=args.band,
                    interval_min=args.interval,
                )
            except Exception as exc:
                print(f"Error: {exc}", file=sys.stderr)
                return 1
            if not (args.expected_tsys or args.smooth or args.interpolate):
                return 0
            antab_path = generated
            out_path = out_path or generated

        if not antab_path:
            parser.error("Provide an .antab file path, or use --from-fits.")
        if not os.path.exists(antab_path):
            print(f"File not found: {antab_path}", file=sys.stderr)
            return 1

        segments = parse_antab(antab_path)
        gain_info = _parse_gain_info(antab_path)
        out_path = out_path or antab_path
        station_filter = _normalize_station_key(args.station) if args.station else None

        if args.expected_tsys:
            print("Applying expected TSYS...")
            op_expected_tsys(segments, gain_info, sefd_table, station_filter)

        if args.interpolate:
            print("Interpolating blanks...")
            op_interpolate(segments, station_filter)

        if args.smooth:
            method_raw, window_raw = args.smooth
            method = method_raw.lower()
            if method not in ("mean", "median", "gaussian"):
                print(f"Unknown smooth method '{method_raw}'. Choose: mean, median, gaussian.", file=sys.stderr)
                return 1
            try:
                window_min = float(window_raw)
            except ValueError:
                print(f"Invalid smooth window '{window_raw}' — must be a number (minutes).", file=sys.stderr)
                return 1
            print(f"Smoothing ({method}, {window_min:.1f} min)...")
            op_smooth(segments, method, window_min, station_filter)

        write_antab(out_path, segments)
        print(f"Saved: {out_path}")
        return 0

    # ------------------------------------------------------------------ #
    #  GUI path                                                            #
    # ------------------------------------------------------------------ #
    _import_gui_deps()

    path = args.path

    root = tk.Tk()
    if path and os.path.exists(path):
        segments = parse_antab(path)
        AntabGui(root, path, segments)
    else:
        # No file given — open the editor in empty state.
        # Use "Open File…" or "Generate Blank ANTAB…" from inside the GUI.
        AntabGui(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
