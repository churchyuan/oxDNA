#ifndef NANOPORE_H_
#define NANOPORE_H_

#include "BaseForce.h"
#include "CylindricalGeometry.h"

/**
 * A membrane of full thickness `thickness` pierced by a cylindrical pore.
 * All lengths describe physical surfaces. `particle_radius` (default 0) adds
 * spherical centre clearance. Uses absolute coordinates, without PBC wrapping.
 * CPU MC/MD and CUDA MD share the same conservative harmonic penalty.
 */
class Nanopore: public BaseForce {
public:
	number _radius = 1.;
	number _half_thickness = 0.5;
	number _particle_radius = 0.;

public:
	Nanopore() = default;
	~Nanopore() override = default;

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;
	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;
	CylindricalGeometry::Parameters<number> geometry() const;
};

#endif
