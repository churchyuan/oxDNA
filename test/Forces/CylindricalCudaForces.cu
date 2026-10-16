#include "../../src/Particles/BaseParticle.h"
#include "../../src/CUDA/CUDAForces.h"
#include "../../src/CUDA/cuda_utils/CUDABox.h"
#include "../../src/CUDA/Backends/MD_CUDABackend.h"
#include "../../src/Utilities/Logger.h"
#include "../../src/Utilities/Timings.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

// Exercise the production dispatcher linked from oxdna_common, including
// force tags, CPU parameter transfer and accumulation of multiple forces.
__global__ void set_external_forces(c_number4 *, GPU_quat *, CUDA_trap *, c_number4 *, c_number4 *, llint, int, CUDABox *);

class ExternalForceBackendProbe: public MD_CUDABackend {
public:
	void initialise_symbols(int n) {
		_particles.resize(n, nullptr);
		_dt = 0.001;
		_sqr_verlet_skin = 0.04;
		// Device constants have translation-unit linkage without CUDA RDC.
		// Initialise them through the same library method as the real backend.
		_init_CUDA_MD_symbols();
	}
};

__global__ void geometry_energies(c_number4 *positions, CUDA_trap *traps, c_number *energies, int n) {
	const int i = blockIdx.x*blockDim.x + threadIdx.x;
	if(i >= n) return;
	c_number energy = 0;
	for(int j = 0; j < 2; ++j) {
		const auto trap = traps[j*n + i];
		const auto pos = positions[i];
		if(trap.type == CUDA_NANOPORE) energy += CylindricalGeometry::nanopore(trap.nanopore.geometry, pos.x, pos.y, pos.z).energy;
		if(trap.type == CUDA_CYLINDER_CONFINEMENT) energy += CylindricalGeometry::cylinder(trap.cylinder_confinement.geometry, pos.x, pos.y, pos.z).energy;
	}
	energies[i] = energy;
}

int main() {
	int device_count = 0;
	const auto status = cudaGetDeviceCount(&device_count);
	if(status != cudaSuccess || device_count == 0) {
		std::cerr << "SKIPPED: CUDA device unavailable: " << cudaGetErrorString(status) << '\n';
		return 77;
	}
	const int n = 4096;
	Logger::init();
	Logger::instance()->disable_log("cylindrical CUDA tests");
	TimingManager::init();
	std::unique_ptr<ExternalForceBackendProbe> backend(new ExternalForceBackendProbe());
	backend->initialise_symbols(n);
	Nanopore pore;
	pore._stiff = 7.; pore._radius = 2.3; pore._half_thickness = 3.1; pore._particle_radius = 0.4;
	pore._direction = LR_vector(2., -3., 4.); pore._direction.normalize();
	pore._pos0 = LR_vector(1.2, -0.7, 0.3);
	CylinderConfinement cylinder;
	cylinder._stiff = pore._stiff; cylinder._radius = pore._radius;
	cylinder._half_height = pore._half_thickness; cylinder._particle_radius = pore._particle_radius;
	cylinder._direction = pore._direction; cylinder._pos0 = pore._pos0;
	std::vector<CUDA_trap> traps(2*n);
	std::vector<c_number4> positions(n), forces(n), torques(n);
	std::vector<c_number> energies(n);
	std::vector<LR_vector> expected_forces(n, LR_vector(0., 0., 0.));
	std::vector<number> expected_energies(n, 0.);
	std::mt19937 rng(1907);
	std::uniform_real_distribution<double> random(-7., 7.);
	for(auto &trap : traps) trap.type = CUDA_TRAP_NO_FORCE;
	for(int i = 0; i < n; ++i) {
		positions[i] = {c_number(random(rng)), c_number(random(rng)), c_number(random(rng)), 0.f};
		LR_vector pos(positions[i].x, positions[i].y, positions[i].z);
		if(i % 4 == 1 || i % 4 == 3) {
			init_Nanopore_from_CPU(&traps[i].nanopore, &pore);
			expected_forces[i] += pore.value(0, pos);
			expected_energies[i] += pore.potential(0, pos);
		}
		if(i % 4 == 2 || i % 4 == 3) {
			init_CylinderConfinement_from_CPU(&traps[n+i].cylinder_confinement, &cylinder);
			expected_forces[i] += cylinder.value(0, pos);
			expected_energies[i] += cylinder.potential(0, pos);
		}
	}
	c_number4 *d_positions, *d_forces, *d_torques;
	c_number *d_energies;
	CUDA_trap *d_traps;
	CUDA_SAFE_CALL(cudaMalloc(reinterpret_cast<void **>(&d_positions), n*sizeof(c_number4)));
	CUDA_SAFE_CALL(cudaMalloc(reinterpret_cast<void **>(&d_forces), n*sizeof(c_number4)));
	CUDA_SAFE_CALL(cudaMalloc(reinterpret_cast<void **>(&d_torques), n*sizeof(c_number4)));
	CUDA_SAFE_CALL(cudaMalloc(reinterpret_cast<void **>(&d_energies), n*sizeof(c_number)));
	CUDA_SAFE_CALL(cudaMalloc(reinterpret_cast<void **>(&d_traps), 2*n*sizeof(CUDA_trap)));
	CUDA_SAFE_CALL(cudaMemcpy(d_positions, positions.data(), n*sizeof(c_number4), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_traps, traps.data(), 2*n*sizeof(CUDA_trap), cudaMemcpyHostToDevice));
	set_external_forces<<<(n+127)/128, 128>>>(d_positions, nullptr, d_traps, d_forces, d_torques, 0, 2, nullptr);
	CUDA_SAFE_CALL(cudaGetLastError());
	geometry_energies<<<(n+127)/128, 128>>>(d_positions, d_traps, d_energies, n);
	CUDA_SAFE_CALL(cudaGetLastError());
	CUDA_SAFE_CALL(cudaMemcpy(forces.data(), d_forces, n*sizeof(c_number4), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(torques.data(), d_torques, n*sizeof(c_number4), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(energies.data(), d_energies, n*sizeof(c_number), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaFree(d_positions)); CUDA_SAFE_CALL(cudaFree(d_forces));
	CUDA_SAFE_CALL(cudaFree(d_torques)); CUDA_SAFE_CALL(cudaFree(d_energies)); CUDA_SAFE_CALL(cudaFree(d_traps));
	const double tolerance = sizeof(c_number) == sizeof(float) ? 2.e-4 : 2.e-10;
	for(int i = 0; i < n; ++i) {
		const double actual[] = {forces[i].x, forces[i].y, forces[i].z, energies[i], torques[i].x, torques[i].y, torques[i].z};
		const double expected[] = {expected_forces[i].x, expected_forces[i].y, expected_forces[i].z, expected_energies[i], 0., 0., 0.};
		for(int j = 0; j < 7; ++j) {
			if(!std::isfinite(actual[j]) || std::abs(actual[j]-expected[j]) > tolerance*(1.+std::abs(expected[j]))) {
				std::cerr << "CUDA mismatch at particle " << i << ", component " << j << ": " << actual[j] << " vs " << expected[j] << '\n';
				return 1;
			}
		}
	}
	std::cout << "Passed " << n << " CPU/CUDA force and energy comparisons (including empty slots and summed forces).\n";
	backend.reset();
	TimingManager::clear();
	return 0;
}
