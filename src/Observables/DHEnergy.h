/*
 * DHEnergy.h
 *
 *  Created on: Mar 24, 2026
 */

#ifndef DHENERGY_H_
#define DHENERGY_H_

#include "BaseObservable.h"

class DHEnergy: public BaseObservable {
protected:
	int _term_id;
	bool _per_particle;

public:
	DHEnergy();
	virtual ~DHEnergy();

	virtual void get_settings(input_file &my_inp, input_file &sim_inp);
	std::string get_output_string(llint curr_step);
};

#endif /* DHENERGY_H_ */
