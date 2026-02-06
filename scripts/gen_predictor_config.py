#!/usr/bin/env python3
"""Generate build/predictor.mk from YAML config file."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path


def parse_value(text: str):
    text = text.strip()
    if not text:
        return ""
    if text[0] in ('"', "'") and text[-1] == text[0]:
        return text[1:-1]
    if text.lower() in ("true", "false"):
        return text.lower() == "true"
    try:
        return int(text, 10)
    except ValueError:
        return text


def parse_simple_yaml(path: Path) -> dict:
    data: dict = {}
    current_list_key = None

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.split("#", 1)[0].rstrip()
        if not line.strip():
            continue

        stripped = line.lstrip()
        if current_list_key and stripped.startswith("- "):
            item_text = stripped[2:].strip()
            data[current_list_key].append(parse_value(item_text))
            continue

        current_list_key = None
        if ":" not in line:
            raise ValueError(f"Unsupported line in {path}: {raw_line!r}")

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()

        if not value:
            data[key] = []
            current_list_key = key
            continue

        if value.startswith("[") and value.endswith("]"):
            parsed = ast.literal_eval(value)
            if not isinstance(parsed, list):
                raise ValueError(f"{key} must be a list in {path}")
            data[key] = parsed
            continue

        data[key] = parse_value(value)

    return data


def render_predictor_type(cfg: dict) -> str:
    if "predictor_type" in cfg:
        predictor_type = str(cfg["predictor_type"]).strip()
        if not predictor_type:
            raise ValueError("predictor_type cannot be empty")
        return predictor_type

    predictor = str(cfg.get("predictor", "tage")).strip()
    if not predictor:
        raise ValueError("predictor cannot be empty")

    args = cfg.get("template_args", [])
    if args is None:
        args = []
    if not isinstance(args, list):
        raise ValueError("template_args must be a YAML list")
    args_str = ",".join(str(item).strip() for item in args)
    return f"{predictor}<{args_str}>"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input YAML file")
    parser.add_argument("--output", required=True, help="Output Makefile fragment")
    args = parser.parse_args()

    cfg_path = Path(args.input)
    out_path = Path(args.output)
    cfg = parse_simple_yaml(cfg_path)
    predictor_type = render_predictor_type(cfg)
    trace = str(cfg.get("trace", "./gcc_test_trace.gz")).strip()
    trace_name = str(cfg.get("trace_name", "test")).strip()
    warmup = int(cfg.get("warmup", 1_000_000))
    measure = int(cfg.get("measure", 40_000_000))
    if not trace:
        raise ValueError("trace cannot be empty")
    if not trace_name:
        raise ValueError("trace_name cannot be empty")
    if warmup < 0 or measure < 0:
        raise ValueError("warmup and measure must be non-negative")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        (
            "# Auto-generated from "
            + str(cfg_path)
            + "\n"
            + f"PREDICTOR_TYPE := {predictor_type}\n"
            + f"TRACE := {trace}\n"
            + f"TRACE_NAME := {trace_name}\n"
            + f"WARMUP := {warmup}\n"
            + f"MEASURE := {measure}\n"
        ),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
