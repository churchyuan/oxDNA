#ifndef NANOPORE_H_
#define NANOPORE_H_

#include "BaseForce.h"

class Nanopore: public BaseForce {
public:
	number _radius;
	number _half_thickness;

public:
	Nanopore();
	virtual ~Nanopore();

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;
	LR_vector value(llint step, LR_vector &pos) override;
	number potential(llint step, LR_vector &pos) override;
};

#endif
