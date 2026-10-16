#!/usr/bin/env python3
"""End-to-end MD/MC, NVE convergence, restart and fix_diffusion checks.

Usage: python3 test/Forces/cylindrical_smoke.py build/bin/oxDNA
All generated files are kept in a temporary directory, removed on completion.
"""
import math
from pathlib import Path
import subprocess
import sys
import tempfile


def run(executable, folder, kind, point, *, dt=0.001, steps=1000,
        center=(0, 0, 0), axis=(0, 0, 1), fix_diffusion=1,
        restart=None, mc=False, backend="CPU", precision=None):
    folder.mkdir()
    (folder / "particle.top").write_text("1 1\n1 A -1 -1\n")
    if restart is None:
        (folder / "initial.dat").write_text(
            "t = 0\nb = 40 40 40\nE = 0 0 0\n" +
            " ".join(map(str, (*point, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0.1, 0.2, 0.3))) + "\n")
    else:
        (folder / "initial.dat").write_bytes(restart.read_bytes())
    vector = lambda v: ",".join(map(str, v))
    length_key = "thickness" if kind == "nanopore" else "height"
    (folder / "force.txt").write_text(
        "{\ntype=" + kind + "\nparticle=all\nstiff=100\nradius=4\n" + length_key +
        "=10\nparticle_radius=0.4\ncenter=" + vector(center) + "\ndir=" + vector(axis) +
        "\ngroup_name=geometry\n}\n")
    simulation = ("sim_type=MC\nensemble=nvt\ndelta_translation=0.05\ndelta_rotation=0.05\n"
                  "check_energy_every=10\ncheck_energy_threshold=1e-8\n") if mc else (
        "sim_type=MD\nthermostat=no\nrefresh_vel=0\ndt=" + str(dt) + "\n")
    (folder / "input").write_text(
        "interaction_type=DNA2\nbackend=" + backend + "\n" +
        ("backend_precision=" + precision + "\n" if precision else "") +
        "T=0.1\nsalt_concentration=0.15\n" + simulation +
        "steps=" + str(steps) + "\nseed=1907\nrestart_step_counter=" + ("0" if restart else "1") +
        "\nfix_diffusion=" + str(fix_diffusion) + "\nverlet_skin=0.2\nconf_file=initial.dat\n"
        "topology=particle.top\nexternal_forces=1\nexternal_forces_file=force.txt\n"
        "trajectory_file=trajectory.dat\nlastconf_file=last.dat\nenergy_file=energy.dat\n"
        "print_conf_interval=100\nprint_energy_every=100\ntime_scale=linear\nno_stdout_energy=1\n"
        "data_output_1={\nname=conserved.dat\nprint_every=10\n"
        "col_1={\ntype=step\n}\ncol_2={\ntype=kinetic_energy\nprecision=12\n}\n"
        "col_3={\ntype=force_energy\nprint_group=geometry\nprecision=12\n}\n}\n")
    result = subprocess.run([str(executable), "input"], cwd=str(folder),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if result.returncode:
        raise AssertionError(folder.name + ": simulation failed\n" + result.stderr[-5000:])
    rows = [[float(x) for x in line.split()] for line in (folder / "conserved.dat").read_text().splitlines()]
    assert rows and all(math.isfinite(x) for row in rows for x in row), folder.name + ": nonfinite output"
    energies = [row[1] + row[2] for row in rows]
    drift = max(abs(x - energies[0]) for x in energies)
    if not mc:
        assert drift < 0.003 * (1 + abs(energies[0])), (folder.name, drift)
    return drift, (folder / "last.dat"), rows[0][2]


def final_state(path):
    return [float(x) for x in path.read_text().splitlines()[3].split()]


def main():
    executable = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="oxdna-cylindrical-") as temporary:
        root = Path(temporary)
        cases = [
            ("cylinder_side", "cylinder_confinement", (3.8, 0, 0), 2.),
            ("cylinder_cap", "cylinder_confinement", (0, 0, 4.8), 2.),
            ("cylinder_corner", "cylinder_confinement", (3.8, 0, 4.8), 4.),
            ("pore_side", "nanopore", (3.8, 0, 0), 2.),
            ("pore_face", "nanopore", (6, 0, 5.2), 2.),
            ("pore_rim", "nanopore", (3.85, 0, 5.2), 1.125),
            ("pore_reservoir", "nanopore", (0, 0, 6), 0.),
        ]
        for name, kind, point, expected_energy in cases:
            drift, _, initial_energy = run(executable, root / name, kind, point)
            assert abs(initial_energy - expected_energy) < 1e-9, (name, initial_energy)
            print("{}: initial U={:.6f}, max |delta(K+U)|={:.3g}".format(name, initial_energy, drift))
        _, first_step, _ = run(executable, root / "first_step", "cylinder_confinement", (3.8, 0, 0), steps=1)
        state = final_state(first_step)
        assert abs(state[0] - 3.79999) < 1e-12, ("initial Verlet position", state[0])
        assert abs(state[9] - (-0.0199995)) < 1e-12, ("initial Verlet velocity", state[9])
        for name, kind, point, _ in (cases[0], cases[3]):
            run(executable, root / (name + "_mc"), kind, point, mc=True)
        coarse, _, _ = run(executable, root / "coarse", "cylinder_confinement", (3.8, 0, 4.8), dt=0.004, steps=500)
        fine, _, _ = run(executable, root / "fine", "cylinder_confinement", (3.8, 0, 4.8), dt=0.002, steps=1000)
        assert fine < 0.6*coarse + 1e-9, ("NVE timestep convergence", coarse, fine)
        for name, kind, _, _ in (cases[0], cases[3]):
            center, axis = (48, -22, 17), (0.6, 0, 0.8)
            point = (center[0] + 0.8*3.8, center[1], center[2] - 0.6*3.8)
            _, shifted, _ = run(executable, root / (name + "_shifted"), kind, point, center=center, axis=axis)
            _, unshifted, _ = run(executable, root / (name + "_unshifted"), kind, point, center=center, axis=axis, fix_diffusion=0)
            assert max(abs(a-b) for a, b in zip(final_state(shifted), final_state(unshifted))) < 1e-9, "fix_diffusion changed physics"
            _, first, _ = run(executable, root / (name + "_first"), kind, point, center=center, axis=axis, steps=50)
            _, resumed, _ = run(executable, root / (name + "_resumed"), kind, point, center=center, axis=axis, steps=950, restart=first)
            assert max(abs(a-b) for a, b in zip(final_state(shifted), final_state(resumed))) < 1e-9, "restart changed physics"
        print("Passed MD/MC integration, step convergence (drift ratio {:.3f}), rotated geometry, fix_diffusion and restart checks.".format(fine/coarse))


if __name__ == "__main__":
    main()
