#include "CylinderConfinement.h"
#include "CylindricalForceUtils.h"

std::tuple<std::vector<int>, std::string> CylinderConfinement::init(input_file &inp) {
	BaseForce::init(inp);
	CylindricalForceUtils::read_geometry(inp, "CylinderConfinement", "height",
		*this, _radius, _half_height, _particle_radius);
	if(_particle_radius >= _half_height) {
		throw oxDNAException("CylinderConfinement: require particle_radius < height/2 (an open centre region)");
	}
	std::string particles;
	getInputString(&inp, "particle", particles, 1);
	auto ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles, "CylinderConfinement");
	return std::make_tuple(ids, Utils::sformat(
		"CylinderConfinement(stiff=%g, radius=%g, height=%g, particle_radius=%g, center=%g,%g,%g, dir=%g,%g,%g)",
		_stiff, _radius, 2.*_half_height, _particle_radius,
		_pos0.x, _pos0.y, _pos0.z, _direction.x, _direction.y, _direction.z));
}

CylindricalGeometry::Parameters<number> CylinderConfinement::geometry() const {
	return CylindricalForceUtils::parameters(*this, _radius, _half_height, _particle_radius);
}

LR_vector CylinderConfinement::value(llint step, LR_vector &pos) {
	const auto result = CylindricalGeometry::cylinder(geometry(), pos.x, pos.y, pos.z);
	return LR_vector(result.fx, result.fy, result.fz);
}

number CylinderConfinement::potential(llint step, LR_vector &pos) {
	return CylindricalGeometry::cylinder(geometry(), pos.x, pos.y, pos.z).energy;
}
