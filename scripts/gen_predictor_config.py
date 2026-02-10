#!/usr/bin/env python3
"""Generate build/predictor.mk from YAML config file."""

from __future__ import annotations

import argparse
import ast
from pathlib import Path
from typing import Iterable


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


def get_bool(cfg: dict, key: str, default: bool) -> bool:
    value = cfg.get(key, default)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        lowered = value.strip().lower()
        if lowered in ("true", "yes", "1"):
            return True
        if lowered in ("false", "no", "0"):
            return False
    raise ValueError(f"{key} must be a boolean")


def get_list_of_str(cfg: dict, key: str) -> list[str]:
    value = cfg.get(key, [])
    if value is None:
        return []
    if not isinstance(value, list):
        raise ValueError(f"{key} must be a YAML list")
    return [str(item).strip() for item in value if str(item).strip()]


def join_flags(flags: Iterable[str]) -> str:
    return " ".join(flag for flag in flags if flag)


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
    build_type = str(cfg.get("build_type", "release")).strip().lower()
    warnings_as_errors = get_bool(cfg, "warnings_as_errors", True)
    native_arch = get_bool(cfg, "native_arch", False)
    defines = get_list_of_str(cfg, "defines")
    include_dirs = get_list_of_str(cfg, "extra_include_dirs")
    sanitize_cfg = cfg.get("sanitize", "none")
    extra_common_flags = str(cfg.get("extra_common_flags", "")).strip()
    extra_cbp_flags = str(cfg.get("extra_cbp_flags", "")).strip()

    if build_type not in ("release", "debug", "relwithdebinfo"):
        raise ValueError("build_type must be one of: release, debug, relwithdebinfo")

    sanitize_list: list[str]
    if isinstance(sanitize_cfg, str):
        sanitize_list = [sanitize_cfg.strip().lower()] if sanitize_cfg.strip() else ["none"]
    elif isinstance(sanitize_cfg, list):
        sanitize_list = [str(item).strip().lower() for item in sanitize_cfg if str(item).strip()]
    else:
        raise ValueError("sanitize must be a string or list")

    if not sanitize_list:
        sanitize_list = ["none"]
    allowed_sanitize = {"none", "address", "undefined"}
    for san in sanitize_list:
        if san not in allowed_sanitize:
            raise ValueError("sanitize supports only: none, address, undefined")
    if "none" in sanitize_list and len(sanitize_list) > 1:
        raise ValueError("sanitize cannot combine 'none' with other sanitizers")
    if not trace:
        raise ValueError("trace cannot be empty")
    if not trace_name:
        raise ValueError("trace_name cannot be empty")
    if warmup < 0 or measure < 0:
        raise ValueError("warmup and measure must be non-negative")

# some extra flags are needed for the predictor and/or CBP builds, 
# so we generate them here based on the config
    generated_common_flags: list[str] = []
    if build_type == "debug":
        generated_common_flags.extend(["-O0", "-g"])
    elif build_type == "relwithdebinfo":
        generated_common_flags.extend(["-O2", "-g"])
    if native_arch:
        generated_common_flags.append("-march=native")
    if "address" in sanitize_list:
        generated_common_flags.extend(["-fsanitize=address", "-fno-omit-frame-pointer"])
    if "undefined" in sanitize_list:
        generated_common_flags.extend(["-fsanitize=undefined", "-fno-omit-frame-pointer"])
    if extra_common_flags:
        generated_common_flags.append(extra_common_flags)

    generated_cbp_flags: list[str] = []
    generated_cbp_flags.extend([f"-D{macro}" for macro in defines])
    generated_cbp_flags.extend([f"-I{inc}" for inc in include_dirs])
    if extra_cbp_flags:
        generated_cbp_flags.append(extra_cbp_flags)

    cbp_warn_flags = [
        "-Wall",
        "-Wextra",
        "-pedantic",
        "-Wold-style-cast",
        "-Wno-deprecated-declarations",
        "-Wno-mismatched-tags",
    ]
    if warnings_as_errors:
        cbp_warn_flags.insert(4, "-Werror")

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
            + f"CBP_WARN_FLAGS := {join_flags(cbp_warn_flags)}\n"
            + f"EXTRA_COMMON_FLAGS := {join_flags(generated_common_flags)}\n"
            + f"EXTRA_CBP_FLAGS := {join_flags(generated_cbp_flags)}\n"
        ),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
