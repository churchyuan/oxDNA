#!/usr/bin/env python3
"""Compare full CPU/CUDA MD runs, fixed frames and GPU restarts.

Usage: python3 test/Forces/cylindrical_cuda_smoke.py build/bin/oxDNA --precision mixed
Repeat with float, and with double on a CUDA_DOUBLE=ON build.
"""
import argparse
from pathlib import Path
import tempfile

from cylindrical_smoke import final_state, run


def translation_error(first, second):
    a, b = final_state(first), final_state(second)
    return max(abs(a[i] - b[i]) for i in (0, 1, 2, 9, 10, 11))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--precision", choices=("mixed", "float", "double"), default="mixed")
    args = parser.parse_args()
    executable = args.executable.resolve()
    # The existing CUDA backend stores its timestep in float even in double
    # builds; allow that small trajectory difference independently of force
    # comparison tolerances in CylindricalCudaForces.cu.
    tolerance = {"mixed": 2e-4, "float": 2e-3, "double": 2e-5}[args.precision]
    gpu = {"backend": "CUDA", "precision": args.precision}
    cases = [
        ("cylinder_side", "cylinder_confinement", (3.8, 0, 0), 2.),
        ("cylinder_cap", "cylinder_confinement", (0, 0, 4.8), 2.),
        ("cylinder_corner", "cylinder_confinement", (3.8, 0, 4.8), 4.),
        ("pore_side", "nanopore", (3.8, 0, 0), 2.),
        ("pore_face", "nanopore", (6, 0, 5.2), 2.),
        ("pore_rim", "nanopore", (3.85, 0, 5.2), 1.125),
        ("pore_reservoir", "nanopore", (0, 0, 6), 0.),
    ]
    with tempfile.TemporaryDirectory(prefix="oxdna-cylindrical-cuda-") as temporary:
        root = Path(temporary)
        for name, kind, point, expected_energy in cases:
            _, cpu_state, _ = run(executable, root/(name+"_cpu"), kind, point)
            drift, gpu_state, initial_energy = run(executable, root/(name+"_gpu"), kind, point, **gpu)
            error = translation_error(cpu_state, gpu_state)
            assert abs(initial_energy-expected_energy) < tolerance, (name, "initial energy", initial_energy)
            assert error < tolerance, (name, "CPU/CUDA trajectory", error, tolerance)
            print("{} [{}]: initial U={:.9f}, NVE drift={:.3g}, CPU/CUDA translation error={:.3g}".format(
                name, args.precision, initial_energy, drift, error))
        for name, kind, _, _ in (cases[0], cases[3]):
            center, axis = (48, -22, 17), (0.6, 0, 0.8)
            point = (center[0]+0.8*3.8, center[1], center[2]-0.6*3.8)
            geometry = {"center": center, "axis": axis}
            _, full, _ = run(executable, root/(name+"_full"), kind, point, **geometry, **gpu)
            _, unshifted, _ = run(executable, root/(name+"_unshifted"), kind, point,
                                  fix_diffusion=0, **geometry, **gpu)
            _, first, _ = run(executable, root/(name+"_first"), kind, point, steps=50, **geometry, **gpu)
            _, resumed, _ = run(executable, root/(name+"_resumed"), kind, point,
                                steps=950, restart=first, **geometry, **gpu)
            assert translation_error(full, unshifted) < tolerance, (name, "CUDA fix_diffusion")
            assert translation_error(full, resumed) < tolerance, (name, "CUDA restart")
        print("Passed full CPU/CUDA MD comparisons, tilted absolute geometry, fix_diffusion and restart checks ({}).".format(args.precision))


if __name__ == "__main__":
    main()
