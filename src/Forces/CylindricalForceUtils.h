#ifndef OXDNA_CYLINDRICAL_FORCE_UTILS_H_
#define OXDNA_CYLINDRICAL_FORCE_UTILS_H_

#include "BaseForce.h"
#include "CylindricalGeometry.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <limits>

namespace CylindricalForceUtils {

inline number parse_number(const std::string &text, const char *key, const char *name) {
	char *end = nullptr;
	errno = 0;
	const double parsed = std::strtod(text.c_str(), &end);
	const bool has_number = end != text.c_str();
	while(*end && std::isspace(static_cast<unsigned char>(*end))) ++end;
	// The project builds with -ffast-math, which can optimise away isfinite().
	// Inspect the exponent bits so nan/inf are still rejected in release builds.
	static_assert(sizeof(double) == sizeof(std::uint64_t) && std::numeric_limits<double>::is_iec559,
		"Cylindrical force input validation requires IEEE-754 doubles");
	std::uint64_t bits;
	std::memcpy(&bits, &parsed, sizeof(bits));
	const bool finite = (bits & UINT64_C(0x7ff0000000000000)) != UINT64_C(0x7ff0000000000000);
	if(!has_number || *end || errno == ERANGE || !finite ||
		parsed > std::numeric_limits<number>::max() || parsed < -std::numeric_limits<number>::max()) {
		throw oxDNAException("%s: %s must be a finite number, got '%s'", name, key, text.c_str());
	}
	return static_cast<number>(parsed);
}

inline number read_number(input_file &inp, const char *key, const char *name, bool mandatory, number fallback = 0.) {
	std::string text;
	if(getInputString(&inp, key, text, mandatory) == KEY_NOT_FOUND) return fallback;
	return parse_number(text, key, name);
}

inline LR_vector read_vector(input_file &inp, const char *key, const char *name, bool mandatory, LR_vector fallback) {
	std::string text;
	if(getInputString(&inp, key, text, mandatory) == KEY_NOT_FOUND) return fallback;
	const auto first = text.find(',');
	const auto second = first == std::string::npos ? std::string::npos : text.find(',', first + 1);
	if(second == std::string::npos || text.find(',', second + 1) != std::string::npos) {
		throw oxDNAException("%s: %s must contain exactly three comma-separated numbers", name, key);
	}
	return LR_vector(parse_number(text.substr(0, first), key, name),
		parse_number(text.substr(first + 1, second - first - 1), key, name),
		parse_number(text.substr(second + 1), key, name));
}

inline void read_geometry(input_file &inp, const char *name, const char *length_key,
	BaseForce &force, number &radius, number &half_length, number &particle_radius) {
	// GPU barostat energies omit external fields, and the legacy MC volume
	// move does not consistently accumulate them. Do not silently sample an
	// incorrect NPT distribution with these position-dependent boundaries.
	if(CONFIG_INFO->sim_input != nullptr) {
		bool use_barostat = false;
		getInputBool(CONFIG_INFO->sim_input, "use_barostat", &use_barostat, 0);
		std::string ensemble;
		getInputString(CONFIG_INFO->sim_input, "ensemble", ensemble, 0);
		std::transform(ensemble.begin(), ensemble.end(), ensemble.begin(),
			[](unsigned char c) { return static_cast<char>(std::tolower(c)); });
		if(use_barostat || ensemble == "npt") {
			throw oxDNAException("%s requires a fixed simulation box (NVT/NVE): this codebase does not consistently include external-field energy in volume-move acceptance", name);
		}
	}
	force._stiff = read_number(inp, "stiff", name, true);
	radius = read_number(inp, "radius", name, true);
	const number length = read_number(inp, length_key, name, true);
	half_length = number(0.5) * length;
	particle_radius = read_number(inp, "particle_radius", name, false);
	if(force._stiff <= 0. || radius <= 0. || half_length <= 0.) {
		throw oxDNAException("%s: stiff, radius and %s must be strictly positive", name, length_key);
	}
	if(particle_radius < 0. || particle_radius >= radius) {
		throw oxDNAException("%s: require 0 <= particle_radius < radius (an open centre region)", name);
	}
	force._pos0 = read_vector(inp, "center", name, false, LR_vector(0., 0., 0.));
	force._direction = read_vector(inp, "dir", name, true, LR_vector(0., 0., 1.));
	// Scale before normalising to avoid overflow/underflow for a finite input axis.
	const number scale = std::max(std::fabs(force._direction.x),
		std::max(std::fabs(force._direction.y), std::fabs(force._direction.z)));
	if(scale == 0.) throw oxDNAException("%s: dir must be a nonzero vector", name);
	force._direction /= scale;
	force._direction.normalize();
}

inline CylindricalGeometry::Parameters<number> parameters(const BaseForce &force,
	number radius, number half_length, number particle_radius) {
	return {force._stiff, radius, half_length, particle_radius,
		force._pos0.x, force._pos0.y, force._pos0.z,
		force._direction.x, force._direction.y, force._direction.z};
}

} // namespace CylindricalForceUtils
#endif
