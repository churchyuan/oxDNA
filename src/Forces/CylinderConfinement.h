#ifndef OXDNA_CYLINDER_CONFINEMENT_H_
#define OXDNA_CYLINDER_CONFINEMENT_H_

#include "BaseForce.h"
#include "CylindricalGeometry.h"

/**
 * Harmonic confinement inside a finite cylinder with two planar end caps.
 * `height` is the full physical height; particle_radius defaults to zero.
 * Forces act on particle centres in absolute coordinates, without PBC wrapping.
 * This external force is independent of the particle-particle interaction.
 */
class CylinderConfinement: public BaseForce {
public:
	number _radius = 1.;
	number _half_height = 0.5;
	number _particle_radius = 0.;

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;
	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;
	CylindricalGeometry::Parameters<number> geometry() const;
};

#endif
