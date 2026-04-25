#include "TwistConstraintForce.h"

#include "../Particles/BaseParticle.h"
#include "../Utilities/oxDNAException.h"

#include <cctype>

namespace {
LR_vector _parse_vector(input_file &inp, const std::string &key, bool required, const std::string &force_name) {
	std::string raw;
	int found = getInputString(&inp, key.c_str(), raw, required ? 1 : 0);
	if(!required && found != KEY_FOUND) {
		return LR_vector((number) 0., (number) 0., (number) 0.);
	}

	double vals[3];
	int parsed = sscanf(raw.c_str(), "%lf, %lf, %lf", vals, vals + 1, vals + 2);
	if(parsed != 3) {
		throw oxDNAException("Could not parse %s `%s' for %s. Aborting", key.c_str(), raw.c_str(), force_name.c_str());
	}

	return LR_vector((number) vals[0], (number) vals[1], (number) vals[2]);
}

LR_vector _parse_vector_alias(input_file &inp, const std::string &key, const std::string &alias, const std::string &force_name, const LR_vector &default_value, bool *found_key = nullptr) {
	std::string raw;
	int found = getInputString(&inp, key.c_str(), raw, 0);
	if(found != KEY_FOUND) {
		found = getInputString(&inp, alias.c_str(), raw, 0);
	}

	if(found_key != nullptr) {
		*found_key = (found == KEY_FOUND);
	}

	if(found != KEY_FOUND) {
		return default_value;
	}

	double vals[3];
	int parsed = sscanf(raw.c_str(), "%lf, %lf, %lf", vals, vals + 1, vals + 2);
	if(parsed != 3) {
		throw oxDNAException("Could not parse %s/%s `%s' for %s. Aborting", key.c_str(), alias.c_str(), raw.c_str(), force_name.c_str());
	}

	return LR_vector((number) vals[0], (number) vals[1], (number) vals[2]);
}

void _normalise_or_throw(LR_vector &v, const std::string &what, const std::string &force_name) {
	if(v.norm() <= (number) 0.) {
		throw oxDNAException("%s for %s cannot be the zero vector", what.c_str(), force_name.c_str());
	}
	v.normalize();
}

std::string _lower(std::string value) {
	for(char &c : value) {
		c = (char) std::tolower((unsigned char) c);
	}
	return value;
}
}

TwistConstraintForce::TwistConstraintForce() : BaseForce() {
	_mode = MODE_ORBIT;
	_body_axis_selector = AXIS_V3;
	_body_ref_axis_selector = AXIS_V1;
	_body_ref_axis_auto = true;

	_center = LR_vector((number) 0., (number) 0., (number) 0.);
	_orbit_axis = LR_vector((number) 0., (number) 0., (number) 1.);
	_move_dir = LR_vector((number) 0., (number) 0., (number) 0.);
	_force_mask = LR_vector((number) 1., (number) 1., (number) 1.);

	_target_axis0 = LR_vector((number) 0., (number) 0., (number) 1.);
	_target_ref0 = LR_vector((number) 1., (number) 0., (number) 0.);
	_has_target_axis = false;
	_has_target_ref = false;

	_orbit_rate = 0.;
	_orbit_phase0 = 0.;
	_move_rate = 0.;
	_axis_stiff = 0.;
	_twist_stiff = 0.;
	_twist_rate = 0.;
	_twist_phase0 = 0.;
}

std::tuple<std::vector<int>, std::string> TwistConstraintForce::init(input_file &inp) {
	BaseForce::init(inp);

	const std::string force_name = "TwistConstraintForce";

	std::string particles_string;
	getInputString(&inp, "particle", particles_string, 1);

	std::string mode_string;
	getInputString(&inp, "mode", mode_string, 1);
	mode_string = _lower(mode_string);
	if(mode_string == "orbit") {
		_mode = MODE_ORBIT;
	}
	else if(mode_string == "frozen") {
		_mode = MODE_FROZEN;
	}
	else {
		throw oxDNAException("Invalid mode `%s' for %s. Supported values are orbit and frozen", mode_string.c_str(), force_name.c_str());
	}

	if(getInputNumber(&inp, "force_stiff", &_stiff, 0) != KEY_FOUND) {
		getInputNumber(&inp, "stiff", &_stiff, 1);
	}

	getInputNumber(&inp, "move_rate", &_move_rate, 0);
	_move_dir = _parse_vector_alias(inp, "move_dir", "dir", force_name, LR_vector((number) 0., (number) 0., (number) 0.));
	if(_move_dir.norm() > (number) 0.) {
		_move_dir.normalize();
	}
	else if(std::fabs(_move_rate) > (number) 0.) {
		throw oxDNAException("move_rate is non-zero for %s but move_dir was not provided", force_name.c_str());
	}

	_pos0 = _parse_vector(inp, "pos0", true, force_name);
	_center = _parse_vector(inp, "center", true, force_name);
	_orbit_axis = _parse_vector_alias(inp, "orbit_axis", "axis", force_name, LR_vector((number) 0., (number) 0., (number) 1.));
	_normalise_or_throw(_orbit_axis, "orbit axis", force_name);

	if(getInputNumber(&inp, "orbit_rate", &_orbit_rate, 0) != KEY_FOUND) {
		getInputNumber(&inp, "rate", &_orbit_rate, 0);
	}
	if(getInputNumber(&inp, "phase0", &_orbit_phase0, 0) != KEY_FOUND) {
		getInputNumber(&inp, "base", &_orbit_phase0, 0);
	}

	_force_mask = _parse_vector_alias(inp, "force_mask", "mask", force_name, LR_vector((number) 1., (number) 1., (number) 1.));

	std::string body_axis_string;
	getInputString(&inp, "body_axis", body_axis_string, 0);
	body_axis_string = _lower(body_axis_string);
	if(body_axis_string == "" || body_axis_string == "v3") {
		_body_axis_selector = AXIS_V3;
	}
	else if(body_axis_string == "v1") {
		_body_axis_selector = AXIS_V1;
	}
	else if(body_axis_string == "v2") {
		_body_axis_selector = AXIS_V2;
	}
	else {
		throw oxDNAException("Invalid body_axis `%s' for %s. Supported values are v1, v2 and v3", body_axis_string.c_str(), force_name.c_str());
	}

	std::string body_ref_axis_string;
	getInputString(&inp, "body_ref_axis", body_ref_axis_string, 0);
	body_ref_axis_string = _lower(body_ref_axis_string);
	if(body_ref_axis_string == "" || body_ref_axis_string == "auto") {
		_body_ref_axis_auto = true;
		if(_body_axis_selector == AXIS_V1) {
			_body_ref_axis_selector = AXIS_V2;
		}
		else if(_body_axis_selector == AXIS_V2) {
			_body_ref_axis_selector = AXIS_V1;
		}
		else {
			_body_ref_axis_selector = AXIS_V1;
		}
	}
	else {
		_body_ref_axis_auto = false;
		if(body_ref_axis_string == "v1") {
			_body_ref_axis_selector = AXIS_V1;
		}
		else if(body_ref_axis_string == "v2") {
			_body_ref_axis_selector = AXIS_V2;
		}
		else if(body_ref_axis_string == "v3") {
			_body_ref_axis_selector = AXIS_V3;
		}
		else {
			throw oxDNAException("Invalid body_ref_axis `%s' for %s. Supported values are auto, v1, v2 and v3", body_ref_axis_string.c_str(), force_name.c_str());
		}
	}
	if(_body_ref_axis_selector == _body_axis_selector) {
		throw oxDNAException("body_ref_axis cannot be the same as body_axis for %s", force_name.c_str());
	}

	bool found_target_axis = false;
	_target_axis0 = _parse_vector_alias(inp, "target_axis", "target_dir", force_name, _target_axis0, &found_target_axis);
	_has_target_axis = found_target_axis;
	if(_has_target_axis) {
		_normalise_or_throw(_target_axis0, "target_axis", force_name);
	}

	bool found_target_ref = false;
	_target_ref0 = _parse_vector_alias(inp, "target_ref", "reference_axis", force_name, _target_ref0, &found_target_ref);
	_has_target_ref = found_target_ref;
	if(_has_target_ref && !_has_target_axis) {
		throw oxDNAException("target_ref requires target_axis for %s", force_name.c_str());
	}
	if(_has_target_ref) {
		_target_ref0 -= _target_axis0 * (_target_ref0 * _target_axis0);
		_normalise_or_throw(_target_ref0, "target_ref", force_name);
	}

	getInputNumber(&inp, "axis_stiff", &_axis_stiff, 0);
	getInputNumber(&inp, "twist_stiff", &_twist_stiff, 0);
	if(_axis_stiff == 0. && _twist_stiff > 0. && _has_target_axis) {
		_axis_stiff = _twist_stiff;
	}

	bool found_twist_rate = (getInputNumber(&inp, "twist_rate", &_twist_rate, 0) == KEY_FOUND);
	getInputNumber(&inp, "twist_phase0", &_twist_phase0, 0);
	if(_mode == MODE_ORBIT && _has_target_ref && !found_twist_rate) {
		_twist_rate = _orbit_rate;
	}
	if(_mode == MODE_FROZEN) {
		if(std::fabs(_orbit_rate) > (number) 0.) {
			throw oxDNAException("mode = frozen requires orbit_rate = 0 in %s", force_name.c_str());
		}
		if(std::fabs(_twist_rate) > (number) 0.) {
			throw oxDNAException("mode = frozen requires twist_rate = 0 in %s", force_name.c_str());
		}
	}

	auto particle_ids = Utils::get_particles_from_string(CONFIG_INFO->particles(), particles_string, force_name);
	std::string description = Utils::sformat("TwistConstraintForce (mode=%s, force_stiff=%g, move_rate=%g, orbit_rate=%g, axis_stiff=%g, twist_stiff=%g)",
			mode_string.c_str(), _stiff, _move_rate, _orbit_rate, _axis_stiff, _twist_stiff);

	return std::make_tuple(particle_ids, description);
}

LR_vector TwistConstraintForce::_rotated_vector(const LR_vector &v, const LR_vector &axis, number angle) const {
	number sintheta = sin(angle);
	number costheta = cos(angle);
	number olcos = ((number) 1.) - costheta;

	number xyo = axis.x * axis.y * olcos;
	number xzo = axis.x * axis.z * olcos;
	number yzo = axis.y * axis.z * olcos;
	number xsin = axis.x * sintheta;
	number ysin = axis.y * sintheta;
	number zsin = axis.z * sintheta;

	LR_matrix R(axis.x * axis.x * olcos + costheta, xyo - zsin, xzo + ysin,
			xyo + zsin, axis.y * axis.y * olcos + costheta, yzo - xsin,
			xzo - ysin, yzo + xsin, axis.z * axis.z * olcos + costheta);

	return R * v;
}

LR_vector TwistConstraintForce::_trap_position(llint step) const {
	number orbit_angle = _orbit_phase0;
	if(_mode == MODE_ORBIT) {
		orbit_angle += _orbit_rate * (number) step;
	}

	LR_vector translated_center = _center + _move_dir * (_move_rate * (number) step);
	return _rotated_vector(_pos0 - _center, _orbit_axis, orbit_angle) + translated_center;
}

LR_vector TwistConstraintForce::_target_axis(llint step) const {
	(void) step;
	return _target_axis0;
}

LR_vector TwistConstraintForce::_target_ref(llint step) const {
	number phase = _twist_phase0;
	if(_mode == MODE_ORBIT) {
		phase += _twist_rate * (number) step;
	}
	return _rotated_vector(_target_ref0, _target_axis0, phase);
}

LR_vector TwistConstraintForce::_selected_axis(AxisSelector selector) const {
	if(_current_particle == nullptr) {
		return LR_vector((number) 0., (number) 0., (number) 0.);
	}

	switch(selector) {
	case AXIS_V1:
		return _current_particle->orientationT.v1;
	case AXIS_V2:
		return _current_particle->orientationT.v2;
	case AXIS_V3:
	default:
		return _current_particle->orientationT.v3;
	}
}

LR_vector TwistConstraintForce::value(llint step, LR_vector &pos) {
	LR_vector postrap = _trap_position(step);
	LR_vector dr = pos - postrap;
	dr.x *= _force_mask.x;
	dr.y *= _force_mask.y;
	dr.z *= _force_mask.z;

	return -_stiff * dr;
}

LR_vector TwistConstraintForce::torque(llint step, LR_vector &pos) {
	(void) pos;
	if(_current_particle == nullptr || !_current_particle->is_rigid_body()) {
		return LR_vector((number) 0., (number) 0., (number) 0.);
	}

	LR_vector lab_torque((number) 0., (number) 0., (number) 0.);

	if(_has_target_axis && _axis_stiff > (number) 0.) {
		LR_vector body_axis = _selected_axis(_body_axis_selector);
		LR_vector target_axis = _target_axis(step);
		lab_torque += _axis_stiff * body_axis.cross(target_axis);
	}

	if(_has_target_ref && _twist_stiff > (number) 0.) {
		LR_vector body_ref = _selected_axis(_body_ref_axis_selector);
		LR_vector target_ref = _target_ref(step);
		lab_torque += _twist_stiff * body_ref.cross(target_ref);
	}

	return _current_particle->orientationT * lab_torque;
}

number TwistConstraintForce::potential(llint step, LR_vector &pos) {
	LR_vector dr = pos - _trap_position(step);
	dr.x *= _force_mask.x;
	dr.y *= _force_mask.y;
	dr.z *= _force_mask.z;

	number U = (number) (0.5 * _stiff * (dr * dr));

	if(_current_particle == nullptr || !_current_particle->is_rigid_body()) {
		return U;
	}

	if(_has_target_axis && _axis_stiff > (number) 0.) {
		U += _axis_stiff * (((number) 1.) - (_selected_axis(_body_axis_selector) * _target_axis(step)));
	}
	if(_has_target_ref && _twist_stiff > (number) 0.) {
		U += _twist_stiff * (((number) 1.) - (_selected_axis(_body_ref_axis_selector) * _target_ref(step)));
	}

	return U;
}
