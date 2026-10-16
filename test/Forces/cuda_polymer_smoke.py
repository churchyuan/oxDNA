#!/usr/bin/env python3
"""DNA2 confinement MD and RNA2 CUDA dispatch regression on existing fixtures."""
import argparse
import math
from pathlib import Path
import subprocess
import tempfile


def execute(executable, folder, settings):
    folder.mkdir()
    (folder/"input").write_text(settings)
    result = subprocess.run([str(executable), "input"], cwd=str(folder),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    log = folder/"simulation.log"
    assert result.returncode == 0, (folder.name, result.returncode, result.stderr[-5000:],
                                    log.read_text()[-5000:] if log.exists() else "no log")


def common(root, backend, precision):
    return ("sim_type=MD\nbackend=" + backend + "\n" +
            ("backend_precision=" + precision + "\n" if backend == "CUDA" else "") +
            "T=0.1\nsalt_concentration=0.15\ndt=0.0005\nthermostat=no\nseed=1907\n"
            "restart_step_counter=1\nfix_diffusion=1\nverlet_skin=0.1\nrefresh_vel=0\n"
            "trajectory_file=trajectory.dat\nlastconf_file=last.dat\nenergy_file=energy.dat\n"
            "log_file=simulation.log\nprint_conf_interval=100\nprint_energy_every=100\n"
            "time_scale=linear\nno_stdout_energy=1\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("executable", type=Path)
    parser.add_argument("--precision", choices=("mixed", "float", "double"), default="mixed")
    args = parser.parse_args()
    executable = args.executable.resolve()
    root = Path(__file__).resolve().parents[2]
    with tempfile.TemporaryDirectory(prefix="oxdna-polymer-cuda-") as temporary:
        work = Path(temporary)
        for kind, length in (("nanopore", "thickness"), ("cylinder_confinement", "height")):
            folder = work/kind
            # The legacy configuration was produced by MC and has zero momenta.
            settings = common(root, "CUDA", args.precision).replace("refresh_vel=0", "refresh_vel=1") + (
                "interaction_type=DNA2\nsteps=1000\nuse_average_seq=false\nseq_dep_file=" +
                str(root/"oxDNA2_sequence_dependent_parameters.txt") + "\nconf_file=" +
                str(root/"examples/NANOPORE/last_conf_MC.dat") + "\ntopology=" +
                str(root/"examples/NANOPORE/input_conf.top") + "\nexternal_forces=1\n"
                "external_forces_file=force.txt\ndata_output_1={\nname=external.dat\nprint_every=100\n"
                "col_1={\ntype=force_energy\n}\n}\n")
            folder.mkdir()
            (folder/"force.txt").write_text(
                "{\ntype=" + kind + "\nparticle=all\nstiff=100\nradius=4.3\n" + length +
                "=10\nparticle_radius=0.2\ncenter=0,0,0\ndir=0,0,1\n}\n")
            # execute() creates its folder, so use the prepared force file via an absolute path.
            prepared = folder/"force.txt"
            run_folder = folder/"run"
            execute(executable, run_folder, settings.replace("external_forces_file=force.txt", "external_forces_file="+str(prepared)))
            rows = (run_folder/"last.dat").read_text().splitlines()[3:]
            assert len(rows) == 120
            assert all(math.isfinite(float(x)) for row in rows for x in row.split()), kind
            energies = [float(x) for x in (run_folder/"external.dat").read_text().split()]
            assert energies and all(math.isfinite(x) and x >= 0 for x in energies)
            print("DNA2 {} [{}]: 120 nucleotides, 1000 MD steps, finite trajectory/energy, initial external U={:.6f}".format(kind, args.precision, energies[0]))
        # Exercise all RNA launch signatures repaired during CUDA compilation:
        # neighbour matrix, single-buffer edges and accumulated-buffer edges.
        for average in (True, False):
            base = ("interaction_type=RNA2\nsteps=5\nexternal_forces=0\nconf_file=" +
                    str(root/"test/RNA/FORCE_FIELD/init.dat") + "\ntopology=" +
                    str(root/"test/RNA/FORCE_FIELD/init.top") + "\nuse_average_seq=" + str(average).lower() +
                    "\nseq_dep_file=" + str(root/"rna_sequence_dependent_parameters.txt") + "\n")
            cpu = work/("rna_cpu_"+str(average))
            execute(executable, cpu, common(root, "CPU", args.precision)+base)
            # Compare actual dynamics: the existing force_and_torque observable
            # reads host force fields that the CUDA backend does not synchronise.
            reference = [list(map(float, line.split())) for line in (cpu/"last.dat").read_text().splitlines()[3:]]
            assert len(reference) == 16
            modes = (("matrix", 0, 1),) if args.precision == "double" else (
                ("matrix", 0, 1), ("edge1", 1, 1), ("edge2", 1, 2))
            # CUDABaseBackend explicitly rejects use_edge with CUDA_DOUBLE.
            for label, edge, buffers in modes:
                folder = work/("rna_"+str(average)+"_"+label)
                execute(executable, folder, common(root, "CUDA", args.precision)+base+
                        "use_edge="+str(edge)+"\nedge_n_forces="+str(buffers)+"\n")
                values = [list(map(float, line.split())) for line in (folder/"last.dat").read_text().splitlines()[3:]]
                assert len(values) == len(reference)
                error = max(abs(value-expected)/(1+abs(expected)) for row, expected_row in zip(values, reference)
                            for value, expected in zip(row, expected_row))
                assert math.isfinite(error) and error < 1e-3, ("RNA CPU/CUDA dynamics", average, label, error)
                print("RNA2 average={} {} [{}]: CPU/CUDA state scaled error={:.3g}".format(average, label, args.precision, error))
        print("Passed coupled DNA2 confinement MD and RNA2 CUDA dispatch regression.")


if __name__ == "__main__":
    main()
