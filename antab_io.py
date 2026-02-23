#!/usr/bin/env python3
import cmd
import os
import re
import sys
from dataclasses import dataclass, field
from typing import List, Optional, Tuple, Union

INDEX_RE = re.compile(r"'([^']+)'")
DATA_RE = re.compile(r"^\s*(\d{1,3})\s+(\d{2}:\d{2}:\d{2})\s+(.*)$")


@dataclass
class DataRow:
    day: str
    time: str
    values: List[str]


@dataclass
class RawRow:
    line: str


Row = Union[DataRow, RawRow]


@dataclass
class TsysBlock:
    station: str
    header_lines: List[str]
    index: List[str]
    data_rows: List[Row]
    data_end_line: str

    def find_row(self, day: str, time: str) -> Optional[DataRow]:
        day = day.zfill(3)
        for row in self.data_rows:
            if isinstance(row, DataRow) and row.day == day and row.time == time:
                return row
        return None

    def list_days(self) -> List[str]:
        days = []
        for row in self.data_rows:
            if isinstance(row, DataRow):
                days.append(row.day)
        return sorted(set(days))

    def list_times_for_day(self, day: str) -> List[str]:
        day = day.zfill(3)
        times = []
        for row in self.data_rows:
            if isinstance(row, DataRow) and row.day == day:
                times.append(row.time)
        return sorted(times)


@dataclass
class RawSegment:
    lines: List[str]


@dataclass
class TsysSegment:
    block: TsysBlock


Segment = Union[RawSegment, TsysSegment]


def parse_index(lines: List[str]) -> List[str]:
    for line in lines:
        if line.strip().startswith("INDEX"):
            return INDEX_RE.findall(line)
    return []


def parse_data_line(line: str) -> Row:
    m = DATA_RE.match(line)
    if not m:
        return RawRow(line=line)
    day = m.group(1).zfill(3)
    time = m.group(2)
    values = m.group(3).split()
    return DataRow(day=day, time=time, values=values)


def parse_antab(path: str) -> List[Segment]:
    with open(path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()

    segments: List[Segment] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.strip().startswith("TSYS "):
            station = line.strip().split(maxsplit=1)[1] if len(line.strip().split()) > 1 else ""
            header_lines: List[str] = []
            while i < len(lines):
                header_lines.append(lines[i])
                if lines[i].strip() == "/":
                    i += 1
                    break
                i += 1
            index = parse_index(header_lines)
            data_rows: List[Row] = []
            data_end_line = "/"
            while i < len(lines):
                if lines[i].strip() == "/":
                    data_end_line = lines[i]
                    i += 1
                    break
                data_rows.append(parse_data_line(lines[i]))
                i += 1
            segments.append(TsysSegment(block=TsysBlock(
                station=station,
                header_lines=header_lines,
                index=index,
                data_rows=data_rows,
                data_end_line=data_end_line,
            )))
        else:
            buf: List[str] = []
            while i < len(lines) and not lines[i].strip().startswith("TSYS "):
                buf.append(lines[i])
                i += 1
            segments.append(RawSegment(lines=buf))
    return segments


def write_antab(path: str, segments: List[Segment]) -> None:
    out_lines: List[str] = []
    for seg in segments:
        if isinstance(seg, RawSegment):
            out_lines.extend(seg.lines)
            continue
        block = seg.block
        out_lines.extend(block.header_lines)
        for row in block.data_rows:
            if isinstance(row, RawRow):
                out_lines.append(row.line)
            else:
                values = " ".join(row.values)
                out_lines.append(f"{row.day} {row.time} {values}".rstrip())
        out_lines.append(block.data_end_line)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")


def first_antab_in_cwd() -> Optional[str]:
    for name in os.listdir("."):
        if name.lower().endswith(".antab"):
            return name
    return None


class AntabShell(cmd.Cmd):
    intro = "Antab editor. Type 'help' for commands."
    prompt = "antab> "

    def __init__(self, path: str, segments: List[Segment]):
        super().__init__()
        self.path = path
        self.segments = segments
        self.blocks: List[TsysBlock] = [seg.block for seg in segments if isinstance(seg, TsysSegment)]
        self.block_index = 0 if self.blocks else -1

    def _current_block(self) -> Optional[TsysBlock]:
        if self.block_index < 0 or self.block_index >= len(self.blocks):
            return None
        return self.blocks[self.block_index]

    def do_blocks(self, arg: str) -> None:
        """List TSYS blocks."""
        if not self.blocks:
            print("No TSYS blocks found.")
            return
        for i, b in enumerate(self.blocks):
            active = "*" if i == self.block_index else " "
            print(f"{active} {i}: TSYS {b.station} | index={len(b.index)} | rows={len(b.data_rows)}")

    def do_use(self, arg: str) -> None:
        """Use a TSYS block by index. Usage: use <n>"""
        try:
            idx = int(arg.strip())
        except ValueError:
            print("Usage: use <n>")
            return
        if idx < 0 or idx >= len(self.blocks):
            print("Invalid block index.")
            return
        self.block_index = idx
        b = self.blocks[idx]
        print(f"Using TSYS {b.station}")

    def do_index(self, arg: str) -> None:
        """Show index to column mapping. Values start at column 3 (after day/time)."""
        block = self._current_block()
        if not block:
            print("No TSYS block selected.")
            return
        if not block.index:
            print("No INDEX line found in header.")
            return
        for i, name in enumerate(block.index):
            col = i + 3
            print(f"col {col}: {name}")

    def do_days(self, arg: str) -> None:
        """List day-of-year values present."""
        block = self._current_block()
        if not block:
            print("No TSYS block selected.")
            return
        for day in block.list_days():
            print(day)

    def do_list(self, arg: str) -> None:
        """List times for a day. Usage: list <DOY>"""
        block = self._current_block()
        if not block:
            print("No TSYS block selected.")
            return
        day = arg.strip()
        if not day:
            print("Usage: list <DOY>")
            return
        times = block.list_times_for_day(day)
        for t in times:
            print(t)

    def do_show(self, arg: str) -> None:
        """Show a row. Usage: show <DOY> <HH:MM:SS>"""
        block = self._current_block()
        if not block:
            print("No TSYS block selected.")
            return
        parts = arg.split()
        if len(parts) != 2:
            print("Usage: show <DOY> <HH:MM:SS>")
            return
        day, time = parts
        row = block.find_row(day, time)
        if not row:
            print("Row not found.")
            return
        if not block.index:
            print(f"{row.day} {row.time} " + " ".join(row.values))
            return
        for i, name in enumerate(block.index):
            val = row.values[i] if i < len(row.values) else ""
            col = i + 3
            print(f"col {col} {name}: {val}")

    def do_set(self, arg: str) -> None:
        """Set a value by index. Usage: set <DOY> <HH:MM:SS> <INDEX> <VALUE>"""
        block = self._current_block()
        if not block:
            print("No TSYS block selected.")
            return
        parts = arg.split()
        if len(parts) < 4:
            print("Usage: set <DOY> <HH:MM:SS> <INDEX> <VALUE>")
            return
        day, time, idx_name = parts[0], parts[1], parts[2]
        value = " ".join(parts[3:])
        row = block.find_row(day, time)
        if not row:
            print("Row not found.")
            return
        if not block.index:
            print("No INDEX defined; use setcol instead.")
            return
        try:
            idx = block.index.index(idx_name)
        except ValueError:
            print("INDEX not found.")
            return
        while len(row.values) <= idx:
            row.values.append("")
        row.values[idx] = value
        print(f"Set {idx_name} at {row.day} {row.time} to {value}")

    def do_setcol(self, arg: str) -> None:
        """Set a value by column. Usage: setcol <DOY> <HH:MM:SS> <COL> <VALUE>"""
        block = self._current_block()
        if not block:
            print("No TSYS block selected.")
            return
        parts = arg.split()
        if len(parts) < 4:
            print("Usage: setcol <DOY> <HH:MM:SS> <COL> <VALUE>")
            return
        day, time, col_s = parts[0], parts[1], parts[2]
        value = " ".join(parts[3:])
        try:
            col = int(col_s)
        except ValueError:
            print("COL must be an integer.")
            return
        if col < 3:
            print("COL must be >= 3 (1=DOY, 2=UT time).")
            return
        row = block.find_row(day, time)
        if not row:
            print("Row not found.")
            return
        idx = col - 3
        while len(row.values) <= idx:
            row.values.append("")
        row.values[idx] = value
        name = block.index[idx] if block.index and idx < len(block.index) else f"col {col}"
        print(f"Set {name} at {row.day} {row.time} to {value}")

    def do_save(self, arg: str) -> None:
        """Save file. Usage: save [path]"""
        path = arg.strip() or self.path
        write_antab(path, self.segments)
        print(f"Saved to {path}")

    def do_quit(self, arg: str) -> bool:
        """Quit."""
        return True

    def do_exit(self, arg: str) -> bool:
        """Quit."""
        return True

    def do_EOF(self, arg: str) -> bool:
        print("")
        return True


def main(argv: List[str]) -> int:
    path = argv[1] if len(argv) > 1 else first_antab_in_cwd()
    if not path:
        print("No .antab file found and none provided.")
        return 1
    if not os.path.exists(path):
        print(f"File not found: {path}")
        return 1
    segments = parse_antab(path)
    shell = AntabShell(path, segments)
    if shell.blocks:
        b = shell.blocks[shell.block_index]
        print(f"Loaded {path}. Active block: TSYS {b.station}")
        if b.index:
            print("Values start at column 3 (after DOY and UT time). Type 'index' to list mapping.")
    else:
        print(f"Loaded {path}. No TSYS blocks found.")
    shell.cmdloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
