#!/usr/bin/env python3
"""Generate the ring's old-time displacement fields for its initial velocity."""

from __future__ import annotations

import math
import re
import sys
from pathlib import Path


def read_mesh_counts(path: Path) -> tuple[int, int, int]:
    text = path.read_text()
    matches = re.findall(
        r"hex\s+\([^)]*\)\s+\((\d+)\s+(\d+)\s+(\d+)\)", text
    )
    if len(matches) != 1:
        raise RuntimeError(f"Expected one hex block in {path}, found {len(matches)}")
    return tuple(int(value) for value in matches[0])


def read_delta_t(path: Path) -> float:
    text = path.read_text()
    matches = re.findall(r"^\s*deltaT\s+([^;]+);", text, re.MULTILINE)
    if len(matches) != 1:
        raise RuntimeError(f"Expected one deltaT entry in {path}, found {len(matches)}")
    return float(matches[0])


def polygon_centre(points: list[tuple[float, float]]) -> tuple[float, float]:
    twice_area = 0.0
    x_moment = 0.0
    y_moment = 0.0
    for first, second in zip(points, points[1:] + points[:1]):
        cross = first[0] * second[1] - second[0] * first[1]
        twice_area += cross
        x_moment += (first[0] + second[0]) * cross
        y_moment += (first[1] + second[1]) * cross
    return x_moment / (3.0 * twice_area), y_moment / (3.0 * twice_area)


def cell_centres(counts: tuple[int, int, int]) -> list[tuple[float, float]]:
    radial_cells, angular_cells, axial_cells = counts
    inner_radius = 0.95
    outer_radius = 1.05
    quarter_turn = 0.5 * math.pi
    centres = []

    for _ in range(axial_cells):
        for angular_cell in range(angular_cells):
            theta0 = quarter_turn * angular_cell / angular_cells
            theta1 = quarter_turn * (angular_cell + 1) / angular_cells
            for radial_cell in range(radial_cells):
                radius0 = inner_radius + (
                    (outer_radius - inner_radius) * radial_cell / radial_cells
                )
                radius1 = inner_radius + (
                    (outer_radius - inner_radius) * (radial_cell + 1) / radial_cells
                )
                points = [
                    (radius0 * math.cos(theta0), radius0 * math.sin(theta0)),
                    (radius1 * math.cos(theta0), radius1 * math.sin(theta0)),
                    (radius1 * math.cos(theta1), radius1 * math.sin(theta1)),
                    (radius0 * math.cos(theta1), radius0 * math.sin(theta1)),
                ]
                centres.append(polygon_centre(points))
    return centres


def displacements(
    centres: list[tuple[float, float]], delta_t: float, old_time_level: int
) -> list[tuple[float, float, float]]:
    mean_radius = 1.0
    velocity_amplitude = 1.5e-3
    values = []
    for x, y in centres:
        radius = math.hypot(x, y)
        theta = math.atan2(y, x)
        radial_velocity = velocity_amplitude * math.cos(2.0 * theta)
        tangential_velocity = velocity_amplitude * math.sin(2.0 * theta) * (
            -0.5 + 1.5 * (radius - mean_radius) / mean_radius
        )
        scale = -old_time_level * delta_t
        values.append(
            (
                scale
                * (
                    radial_velocity * math.cos(theta)
                    - tangential_velocity * math.sin(theta)
                ),
                scale
                * (
                    radial_velocity * math.sin(theta)
                    + tangential_velocity * math.cos(theta)
                ),
                0.0,
            )
        )
    return values


def write_internal_field(path: Path, values: list[tuple[float, float, float]]) -> None:
    text = path.read_text()
    start = re.search(r"^internalField\s+", text, re.MULTILINE)
    boundary = re.search(r"^boundaryField\s*$", text, re.MULTILINE)
    if start is None or boundary is None or start.start() >= boundary.start():
        raise RuntimeError(f"Could not locate internalField in {path}")

    lines = ["internalField   nonuniform List<vector>", str(len(values)), "("]
    lines.extend(f"({x:.16g} {y:.16g} {z:.16g})" for x, y, z in values)
    lines.extend((")", ";", ""))
    path.write_text(text[: start.start()] + "\n".join(lines) + text[boundary.start() :])


def main() -> None:
    if len(sys.argv) != 4:
        raise SystemExit(
            "Usage: setInitialVelocity.py BLOCK_MESH_DICT CONTROL_DICT FIELD_DIR"
        )
    block_mesh_dict, control_dict, field_dir = map(Path, sys.argv[1:])
    centres = cell_centres(read_mesh_counts(block_mesh_dict))
    delta_t = read_delta_t(control_dict)
    write_internal_field(field_dir / "D_0", displacements(centres, delta_t, 1))
    write_internal_field(field_dir / "D_0_0", displacements(centres, delta_t, 2))


if __name__ == "__main__":
    main()
