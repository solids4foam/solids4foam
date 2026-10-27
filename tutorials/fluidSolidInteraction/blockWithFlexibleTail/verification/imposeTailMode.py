#!/usr/bin/env python3

"""Impose an Euler--Bernoulli tail mode on the fluid mesh-motion fields."""

import argparse
import math
import re
from pathlib import Path


BETAS = {1: 1.875104068711961, 2: 4.694091132974174}
TAIL_ROOT = 0.055
TAIL_LENGTH = 0.040
TAIL_CENTRE = 0.060


def list_body(path):
    text = path.read_text()
    match = re.search(r"\n\s*\d+\s*\n\s*\(\s*\n", text)
    if not match:
        raise ValueError(f"{path}: cannot find OpenFOAM list")
    return text[match.end():text.rfind("\n)")]


def read_points(path):
    pattern = re.compile(
        r"^\s*\(\s*([^\s()]+)\s+([^\s()]+)\s+([^\s()]+)\s*\)\s*$",
        re.MULTILINE,
    )
    return [tuple(map(float, match.groups())) for match in pattern.finditer(
        list_body(path)
    )]


def read_faces(path):
    pattern = re.compile(r"^\s*\d+\s*\(([^()]*)\)\s*$", re.MULTILINE)
    return [
        [int(label) for label in match.group(1).split()]
        for match in pattern.finditer(list_body(path))
    ]


def tail_faces(case):
    mesh = case/"constant"/"fluid"/"polyMesh"
    boundary = (mesh/"boundary").read_text()
    patch = re.search(r"\btail\s*\{(.*?)\}", boundary, re.DOTALL)
    if not patch:
        raise ValueError(f"{mesh/'boundary'}: cannot find tail patch")
    count = int(re.search(r"nFaces\s+(\d+)", patch.group(1)).group(1))
    start = int(re.search(r"startFace\s+(\d+)", patch.group(1)).group(1))
    return read_faces(mesh/"faces")[start:start + count]


def mode_displacement(point, mode, amplitude):
    beta = BETAS[mode]
    sigma = (math.cosh(beta) + math.cos(beta))/(
        math.sinh(beta) + math.sin(beta)
    )
    position = min(max((point[0] - TAIL_ROOT)/TAIL_LENGTH, 0.0), 1.0)

    def shape(value):
        argument = beta*value
        return (
            math.cosh(argument) - math.cos(argument)
            - sigma*(math.sinh(argument) - math.sin(argument))
        )

    normalisation = shape(1.0)
    argument = beta*position
    derivative = beta*(
        math.sinh(argument) + math.sin(argument)
        - sigma*(math.cosh(argument) - math.cos(argument))
    )/normalisation
    transverse = amplitude*shape(position)/normalisation
    axial = -(point[1] - TAIL_CENTRE)*amplitude*derivative/TAIL_LENGTH
    return axial, transverse, 0.0


def replace_tail_value(path, values):
    text = path.read_text()
    match = re.search(r"(\n\s*tail\s*\{)(.*?)(\n\s*\})", text, re.DOTALL)
    if not match:
        raise ValueError(f"{path}: cannot find tail boundary condition")
    value = "value           nonuniform List<vector>\n"
    value += f"{len(values)}\n(\n"
    value += "\n".join(
        f"({component[0]:.12g} {component[1]:.12g} {component[2]:.12g})"
        for component in values
    )
    value += "\n);"
    body, substitutions = re.subn(
        r"value\s+(?:uniform\s+\([^;]+\)|"
        r"nonuniform\s+List<vector>\s+\d+\s*\(.*?\))\s*;",
        value,
        match.group(2),
        count=1,
        flags=re.DOTALL,
    )
    if substitutions != 1:
        raise ValueError(f"{path}: cannot replace tail value")
    path.write_text(text[:match.start()] + match.group(1) + body
                    + match.group(3) + text[match.end():])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, default=Path("."))
    parser.add_argument(
        "--time",
        default="0",
        help="time directory containing the D and pointD fields (default: 0)",
    )
    parser.add_argument("--mode", type=int, choices=BETAS, required=True)
    parser.add_argument("--amplitude", type=float, required=True)
    args = parser.parse_args()

    case = args.case.resolve()
    points = read_points(case/"constant"/"fluid"/"polyMesh"/"points")
    faces = tail_faces(case)

    patch_points = []
    for face in faces:
        for point in face:
            if point not in patch_points:
                patch_points.append(point)

    point_values = [
        mode_displacement(points[point], args.mode, args.amplitude)
        for point in patch_points
    ]
    face_values = [
        mode_displacement(
            tuple(sum(points[point][axis] for point in face)/len(face)
                  for axis in range(3)),
            args.mode,
            args.amplitude,
        )
        for face in faces
    ]

    time_directory = case/args.time/"fluid"
    replace_tail_value(time_directory/"pointD", point_values)
    replace_tail_value(time_directory/"D", face_values)
    print(
        f"Imposed mode {args.mode} with {args.amplitude:g} m tip amplitude "
        f"at time {args.time} on {len(patch_points)} points and "
        f"{len(faces)} faces"
    )


if __name__ == "__main__":
    main()
