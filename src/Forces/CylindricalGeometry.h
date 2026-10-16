#ifndef OXDNA_CYLINDRICAL_GEOMETRY_H_
#define OXDNA_CYLINDRICAL_GEOMETRY_H_

#include <cmath>

#ifdef __CUDACC__
#define OXDNA_CYL_HD __host__ __device__
#else
#define OXDNA_CYL_HD
#endif

/** Pure geometry shared by the CPU and CUDA external forces (no box wrapping). */
namespace CylindricalGeometry {

using std::sqrt;
using std::fabs;

template<typename T> struct Parameters {
	T stiff, radius, half_length, particle_radius;
	T cx, cy, cz;
	T nx, ny, nz;
};

template<typename T, typename U>
OXDNA_CYL_HD inline Parameters<T> convert(const Parameters<U> &p) {
	return {T(p.stiff), T(p.radius), T(p.half_length), T(p.particle_radius),
		T(p.cx), T(p.cy), T(p.cz), T(p.nx), T(p.ny), T(p.nz)};
}

template<typename T> struct Evaluation {
	T energy, fx, fy, fz;
};

template<typename T> struct LocalCoordinates {
	T z, rho, px, py, pz, nx, ny, nz;
};

template<typename T>
OXDNA_CYL_HD inline LocalCoordinates<T> local_coordinates(const Parameters<T> &p, T x, T y, T z) {
	// Renormalise after precision conversion as well: CUDA float and double
	// backends must both project onto a unit axis.
	const T inv_norm = T(1) / sqrt(p.nx*p.nx + p.ny*p.ny + p.nz*p.nz);
	const T nx = p.nx * inv_norm, ny = p.ny * inv_norm, nz = p.nz * inv_norm;
	const T rx = x - p.cx, ry = y - p.cy, rz = z - p.cz;
	const T axial = nx*rx + ny*ry + nz*rz;
	const T px = rx - axial*nx, py = ry - axial*ny, pz = rz - axial*nz;
	return {axial, sqrt(px*px + py*py + pz*pz), px, py, pz, nx, ny, nz};
}

template<typename T>
OXDNA_CYL_HD inline T sign(T x) {
	// In particular, do not select an arbitrary face at the membrane midplane.
	return x > T(0) ? T(1) : (x < T(0) ? T(-1) : T(0));
}

/**
 * Membrane solid: rho >= R AND |z| <= h. Let s be its signed distance,
 * positive in the pore/reservoirs. U = k/2 max(a-s, 0)^2.
 * The Euclidean distance outside the rim gives the correct spherical-particle
 * clearance there, rather than artificially extending the cylindrical wall.
 */
template<typename T>
OXDNA_CYL_HD inline Evaluation<T> nanopore(const Parameters<T> &p, T x, T y, T z) {
	const auto c = local_coordinates(p, x, y, z);
	const T qr = p.radius - c.rho;
	const T qz = fabs(c.z) - p.half_length;
	T distance, grad_r = T(0), grad_z = T(0);
	if(qr > T(0) || qz > T(0)) {
		const T outside_r = qr > T(0) ? qr : T(0);
		const T outside_z = qz > T(0) ? qz : T(0);
		distance = sqrt(outside_r*outside_r + outside_z*outside_z);
		grad_r = -outside_r / distance;
		grad_z = sign(c.z) * outside_z / distance;
	}
	else {
		// Inside the solid, the nearer pore wall or membrane face is the exit.
		// At equal distances the signed distance is not differentiable; use
		// the mean of the two limiting gradients, preserving reflection symmetry.
		if(qr > qz) {
			distance = qr;
			grad_r = T(-1);
		}
		else if(qz > qr) {
			distance = qz;
			grad_z = sign(c.z);
		}
		else {
			distance = qr;
			grad_r = T(-0.5);
			grad_z = T(0.5) * sign(c.z);
		}
	}
	const T penetration = p.particle_radius - distance;
	if(penetration <= T(0)) return {T(0), T(0), T(0), T(0)};
	const T magnitude = p.stiff * penetration;
	const T radial = c.rho > T(0) ? magnitude * grad_r / c.rho : T(0);
	const T axial = magnitude * grad_z;
	return {T(0.5)*magnitude*penetration,
		radial*c.px + axial*c.nx, radial*c.py + axial*c.ny, radial*c.pz + axial*c.nz};
}

/**
 * Accessible centre region: rho <= R-a AND |z| <= h-a.
 * U = k/2 [max(rho-(R-a),0)^2 + max(|z|-(h-a),0)^2].
 * This squared distance to a convex closed cylinder has a continuous force,
 * including at the cap/sidewall junction: both violations contribute.
 */
template<typename T>
OXDNA_CYL_HD inline Evaluation<T> cylinder(const Parameters<T> &p, T x, T y, T z) {
	const auto c = local_coordinates(p, x, y, z);
	const T dr = c.rho - (p.radius - p.particle_radius);
	const T dz = fabs(c.z) - (p.half_length - p.particle_radius);
	const T radial_violation = dr > T(0) ? dr : T(0);
	const T axial_violation = dz > T(0) ? dz : T(0);
	const T radial = c.rho > T(0) ? -p.stiff * radial_violation / c.rho : T(0);
	const T axial = -p.stiff * axial_violation * sign(c.z);
	return {T(0.5)*p.stiff*(radial_violation*radial_violation + axial_violation*axial_violation),
		radial*c.px + axial*c.nx, radial*c.py + axial*c.ny, radial*c.pz + axial*c.nz};
}

} // namespace CylindricalGeometry

#undef OXDNA_CYL_HD
#endif
