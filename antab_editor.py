#!/usr/bin/env python3
import argparse
import math
import os
import sys
from typing import Dict, List, Optional, Set, Tuple

try:
    import tkinter as tk
    from tkinter import ttk
except Exception as exc:
    print("tkinter is required for this tool.")
    print(f"Import error: {exc}")
    raise SystemExit(1)

try:
    import matplotlib
    matplotlib.use("TkAgg")
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    from matplotlib.widgets import RectangleSelector
except Exception as exc:
    print("matplotlib is required for this tool.")
    print(f"Import error: {exc}")
    raise SystemExit(1)

import numpy as np

from antab_io import parse_antab, write_antab, TsysSegment, TsysBlock, DataRow

MISSING_VALUE = -99.0
SEFD_PATH = os.path.join(os.getcwd(), "sefd_values.txt")


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

def _band_cm_to_hz(cm: float) -> float:
    c = 299792458.0
    wavelength_m = cm / 100.0
    return c / wavelength_m


class AntabGui:
    def __init__(self, root: tk.Tk, path: str, segments) -> None:
        self.root = root
        self.path = path
        self.segments = segments
        self.blocks: List[TsysBlock] = [seg.block for seg in segments if isinstance(seg, TsysSegment)]
        if not self.blocks:
            raise RuntimeError("No TSYS blocks found in file.")

        self.block_index = 0
        self.block = self.blocks[self.block_index]

        self.xs: np.ndarray = np.array([], dtype=float)
        self.rows: List[DataRow] = []
        self.series_cache: Dict[str, np.ndarray] = {}
        self.selected_masks: Dict[str, np.ndarray] = {}
        self.max_highlight_points = 20000
        self.max_table_rows = 3000
        self.table_sync_skipped = False
        self.table_sync_rows = 0
        self.gain_info = self._parse_gain_info(path)
        self.sefd_table = self._load_sefd_table()

        self.series_artists: Dict[str, Tuple] = {}
        self.selected_scatter = None
        self.rect_selector: Optional[RectangleSelector] = None
        self.smooth_enabled = tk.BooleanVar(value=False)
        self.smooth_window_var = tk.StringVar(value="5")
        self.smooth_method_var = tk.StringVar(value="Mean")
        self.save_write_new = tk.BooleanVar(value=False)
        self.save_out_var = tk.StringVar(value="")
        self.table = None
        self.row_iids: List[str] = []
        self.iid_to_row: Dict[str, int] = {}
        self.table_selection_guard = False
        self.table_guard_job = None

        self._build_ui()
        self._load_block(0)

    def _parse_gain_info(self, path: str) -> Dict[str, Dict[str, List[float]]]:
        gain_info: Dict[str, Dict[str, List[float]]] = {}
        if not os.path.exists(path):
            return gain_info
        with open(path, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if line.startswith("GAIN"):
                parts = line.split()
                station = parts[1].upper() if len(parts) > 1 else ""
                dpfu_vals: List[float] = []
                freq_vals: List[float] = []
                i += 1
                while i < len(lines) and lines[i].strip() != "/":
                    l = lines[i].strip()
                    if l.startswith("DPFU"):
                        dpfu_vals = _parse_float_list(l)
                    elif l.startswith("FREQ"):
                        freq_vals = _parse_float_list(l)
                    i += 1
                if station:
                    gain_info[station] = {"dpfu": dpfu_vals, "freq": freq_vals}
            else:
                i += 1
        return gain_info

    def _load_sefd_table(self) -> Dict[str, Dict[float, float]]:
        sefd_table: Dict[str, Dict[float, float]] = {}
        if not os.path.exists(SEFD_PATH):
            return sefd_table
        with open(SEFD_PATH, "r", encoding="utf-8") as f:
            lines = [l.rstrip("\n") for l in f if l.strip() != ""]
        if not lines:
            return sefd_table
        header = [v.strip() for v in lines[0].split("|")[1:]]
        band_cm = []
        for h in header:
            try:
                band_cm.append(float(h))
            except ValueError:
                band_cm.append(None)
        band_hz = [(_band_cm_to_hz(cm) if cm is not None else None) for cm in band_cm]
        for line in lines[1:]:
            parts = [v.strip() for v in line.split("|")]
            if not parts:
                continue
            station = parts[0]
            if station == "":
                continue
            station_key = station.lower()
            sefd_table.setdefault(station_key, {})
            for hz, cell in zip(band_hz, parts[1:]):
                if hz is None:
                    continue
                if cell == "":
                    continue
                try:
                    sefd_table[station_key][hz] = float(cell)
                except ValueError:
                    continue
        return sefd_table

    def _build_ui(self) -> None:
        self.root.title("ANTAB GUI Editor")
        self.root.geometry("1200x750")

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
        ttk.Button(edit_frame, text="Apply Expected TSYS/SEFD", command=self._apply_expected).pack(fill=tk.X, padx=6, pady=(0, 6))

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

        save_frame = ttk.LabelFrame(self.left, text="Save")
        save_frame.pack(fill=tk.X, pady=(8, 0))
        ttk.Button(save_frame, text="Save", command=self._save).pack(fill=tk.X)
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

    def _set_status(self, text: str) -> None:
        self.status_var.set(text)

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

    def _expected_tsys_for_index(self, index_name: str) -> Optional[float]:
        station = self.block.station.upper()
        gain = self.gain_info.get(station)
        if not gain:
            return None
        dpfu_vals = gain.get("dpfu", [])
        freq_vals = gain.get("freq", [])
        if not freq_vals:
            return None
        freq_mhz = sum(freq_vals) / len(freq_vals)
        freq_hz = freq_mhz * 1e6
        sefd_map = self.sefd_table.get(station.lower())
        if not sefd_map:
            return None
        nearest = min(sefd_map.keys(), key=lambda f: abs(f - freq_hz))
        sefd_jy = sefd_map[nearest]
        if not dpfu_vals:
            return None
        if len(dpfu_vals) == 1:
            dpfu = dpfu_vals[0]
        else:
            is_r = "R" in index_name
            is_l = "L" in index_name
            if is_r and not is_l:
                dpfu = dpfu_vals[0]
            elif is_l and not is_r:
                dpfu = dpfu_vals[1] if len(dpfu_vals) > 1 else dpfu_vals[0]
            else:
                dpfu = dpfu_vals[0]
        return dpfu * sefd_jy

    def _apply_expected(self) -> None:
        indices = self._selected_indices()
        if not indices:
            self._set_status("Select indices to apply expected TSYS.")
            return
        rows_mask = self._selected_row_mask(indices)
        if not rows_mask.any():
            rows_mask = np.ones(len(self.rows), dtype=bool)
        rows_idx = np.where(rows_mask)[0]
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
            self._set_status("No expected TSYS available (missing GAIN/SEFD data).")
            return
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Applied expected TSYS to {applied_indices} indices.")

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

        for index_name in selected:
            ys = self._build_series(index_name)
            line, = self.ax.plot(self.xs, ys, linewidth=1.2, alpha=0.8, label=index_name)
            scatter = self.ax.scatter(self.xs, ys, s=22, picker=5)
            scatter.set_gid(index_name)
            self.series_artists[index_name] = (line, scatter)
            if window_days is not None:
                smooth = self._smooth_series(index_name, window_days, method)
                if np.isfinite(smooth).any():
                    self.ax.plot(self.xs, smooth, linestyle="--", linewidth=1.2, alpha=0.9,
                                 color=line.get_color(), label=f"_{index_name}_smooth")

        self.ax.set_xlabel("Day of Year (fractional)")
        self.ax.set_ylabel("TSYS")
        self.ax.set_title(f"TSYS {self.block.station}")
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
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Applied value to {total} points.")

    def _blank_value(self) -> None:
        total = int(sum(mask.sum() for mask in self.selected_masks.values()))
        if total == 0:
            self._set_status("No points selected.")
            return
        blank_value = f"{MISSING_VALUE:.1f}"
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
        self._refresh_table_rows(rows_touched)
        self._update_plot()
        self._set_status(f"Blanked {total} points.")

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
        for index_name in indices:
            smooth = self._smooth_series(index_name, window_days, method)
            orig = self._build_series(index_name)
            mask = np.isfinite(orig) & np.isfinite(smooth)
            if not mask.any():
                continue
            idx = self.block.index.index(index_name)
            for row_idx in np.where(mask)[0]:
                row = self.rows[row_idx]
                self._ensure_value_length(row, idx)
                row.values[idx] = f"{smooth[row_idx]:.1f}"
                rows_touched.add(int(row_idx))
            self.series_cache[index_name] = smooth.copy()
        if rows_touched:
            self._refresh_table_rows(rows_touched)
            self._update_plot()
            self._set_status(f"Applied smoothing ({method}) to {len(indices)} indices.")
        else:
            self._set_status("No points to smooth.")

    def _clear_selection(self) -> None:
        self._clear_selection_masks()
        self._update_selected_scatter()

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


def main(argv: List[str]) -> int:
    parser = argparse.ArgumentParser(description="GUI ANTAB editor with multi-plot selection")
    parser.add_argument("path", nargs="?", help=".antab file path")
    args = parser.parse_args(argv)

    path = args.path or first_antab_in_cwd()
    if not path:
        print("No .antab file found. Provide a path.")
        return 1
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return 1

    segments = parse_antab(path)
    root = tk.Tk()
    try:
        AntabGui(root, path, segments)
    except RuntimeError as exc:
        print(str(exc))
        return 1
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
