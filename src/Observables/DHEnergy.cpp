/*
 * DHEnergy.cpp
 *
 *  Created on: Mar 24, 2026
 */

#include "DHEnergy.h"

#include "../Utilities/oxDNAException.h"

DHEnergy::DHEnergy() :
				_term_id(7),
				_per_particle(true) {
}

DHEnergy::~DHEnergy() {
}

void DHEnergy::get_settings(input_file &my_inp, input_file &sim_inp) {
	BaseObservable::get_settings(my_inp, sim_inp);
	getInputInt(&my_inp, "term", &_term_id, 0);
	getInputBool(&my_inp, "per_particle", &_per_particle, 0);
}

std::string DHEnergy::get_output_string(llint curr_step) {
	number energy = (number) 0.f;
	try {
		energy = _config_info->interaction->get_system_energy_term(_term_id, _config_info->particles(), _config_info->lists);
	}
	catch(oxDNAException &e) {
		energy = (number) 0.f;
	}
	if(_per_particle) {
		energy /= _config_info->N();
	}
	return Utils::sformat(_number_formatter, energy);
}
