#include "Nanopore.h"
#include "CylindricalForceUtils.h"

std::tuple<std::vector<int>, std::string> Nanopore::init(input_file &inp) {
	BaseForce::init(inp);
	CylindricalForceUtils::read_geometry(inp, "Nanopore", "thickness",
		*this, _radius, _half_thickness, _particle_radius);
	std::string particles;
	getInputString(&inp, "particle", particles, 1);
	auto ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles, "Nanopore");
	return std::make_tuple(ids, Utils::sformat(
		"Nanopore(stiff=%g, radius=%g, thickness=%g, particle_radius=%g, center=%g,%g,%g, dir=%g,%g,%g)",
		_stiff, _radius, 2.*_half_thickness, _particle_radius,
		_pos0.x, _pos0.y, _pos0.z, _direction.x, _direction.y, _direction.z));
}

CylindricalGeometry::Parameters<number> Nanopore::geometry() const {
	return CylindricalForceUtils::parameters(*this, _radius, _half_thickness, _particle_radius);
}

LR_vector Nanopore::value(llint step, LR_vector &pos) {
	const auto result = CylindricalGeometry::nanopore(geometry(), pos.x, pos.y, pos.z);
	return LR_vector(result.fx, result.fy, result.fz);
}

number Nanopore::potential(llint step, LR_vector &pos) {
	return CylindricalGeometry::nanopore(geometry(), pos.x, pos.y, pos.z).energy;
}
