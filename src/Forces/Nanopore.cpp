#include "Nanopore.h"

#include "../Particles/BaseParticle.h"

Nanopore::Nanopore() :
				BaseForce(),
				_radius(1.0),
				_half_thickness(0.5) {
}

Nanopore::~Nanopore() {
}

static LR_vector _parse_vec3(const std::string &s, const char *key) {
	double tmpf[3];
	int tmpi = sscanf(s.c_str(), "%lf,%lf,%lf", tmpf, tmpf + 1, tmpf + 2);
	if(tmpi != 3) {
		throw oxDNAException("Could not parse %s %s in external forces file. Aborting", key, s.c_str());
	}
	return LR_vector((number) tmpf[0], (number) tmpf[1], (number) tmpf[2]);
}

std::tuple<std::vector<int>, std::string> Nanopore::init(input_file &inp) {
	BaseForce::init(inp);

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	getInputNumber(&inp, "stiff", &_stiff, 1);
	getInputNumber(&inp, "radius", &_radius, 1);
	number thickness = 0.;
	getInputNumber(&inp, "thickness", &thickness, 1);
	_half_thickness = (number) 0.5 * thickness;

	std::string strdir;
	getInputString(&inp, "dir", strdir, 1);
	_direction = _parse_vec3(strdir, "dir");
	_direction.normalize();

	std::string strcenter;
	getInputString(&inp, "center", strcenter, 0);
	if(!strcenter.empty()) {
		_pos0 = _parse_vec3(strcenter, "center");
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, "Nanopore");
	std::string description = Utils::sformat("Nanopore with stiff = %g, radius = %g, thickness = %g, center = %g,%g,%g, dir = %g,%g,%g",
			_stiff, _radius, thickness, _pos0.x, _pos0.y, _pos0.z, _direction.x, _direction.y, _direction.z);
	return std::make_tuple(particle_ids, description);
}

LR_vector Nanopore::value(llint step, LR_vector &pos) {
	LR_vector r = pos - _pos0;
	number z = _direction * r;
	number absz = (number) fabs((double) z);
	if(absz >= _half_thickness) return LR_vector(0., 0., 0.);

	LR_vector r_perp = r - z * _direction;
	number rho2 = r_perp.norm();
	number rho = (number) sqrt((double) rho2);
	if(rho <= _radius) return LR_vector(0., 0., 0.);

	number d_plane = _half_thickness - absz;
	number d_cyl = rho - _radius;

	if(d_cyl < d_plane) {
		if(rho <= (number) 0.f) return LR_vector(0., 0., 0.);
		return -(_stiff * d_cyl / rho) * r_perp;
	}
	else {
		number sgn = (z >= (number) 0.f) ? (number) 1.f : (number) -1.f;
		return (_stiff * d_plane * sgn) * _direction;
	}
}

number Nanopore::potential(llint step, LR_vector &pos) {
	LR_vector r = pos - _pos0;
	number z = _direction * r;
	number absz = (number) fabs((double) z);
	if(absz >= _half_thickness) return (number) 0.f;

	LR_vector r_perp = r - z * _direction;
	number rho = (number) sqrt((double) r_perp.norm());
	if(rho <= _radius) return (number) 0.f;

	number d_plane = _half_thickness - absz;
	number d_cyl = rho - _radius;

	number d = (d_cyl < d_plane) ? d_cyl : d_plane;
	return (number) (0.5 * _stiff * SQR(d));
}
