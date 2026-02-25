#!/usr/bin/env python3
import argparse
import copy
import math
import os
import sys
from typing import Any, Dict, List, Optional, Set, Tuple

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

    def _load_sefd_table(self) -> Dict[str, Dict[float, float]]:
        sefd_table: Dict[str, Dict[float, float]] = {}
        antab_dir = os.path.dirname(os.path.abspath(self.path))
        script_dir = os.path.dirname(os.path.abspath(__file__))
        candidates = [
            os.path.join(antab_dir, SEFD_FILENAME),
            os.path.join(script_dir, SEFD_FILENAME),
            os.path.join(os.getcwd(), SEFD_FILENAME),
        ]
        sefd_path = next((p for p in candidates if os.path.exists(p)), candidates[0])
        if not os.path.exists(sefd_path):
            return sefd_table
        with open(sefd_path, "r", encoding="utf-8") as f:
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
            station_key = _normalize_station_key(station)
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

        save_frame = ttk.LabelFrame(self.left, text="Save")
        save_frame.pack(fill=tk.X, pady=(8, 0))
        self.undo_button = ttk.Button(save_frame, text="Undo Last Change", command=self._undo_last_change)
        self.undo_button.pack(fill=tk.X, pady=(0, 6))
        self._update_undo_button_state()
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
