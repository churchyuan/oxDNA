#ifndef TWISTCONSTRAINTFORCE_H_
#define TWISTCONSTRAINTFORCE_H_

#include "BaseForce.h"

class TwistConstraintForce: public BaseForce {
public:
	enum Mode {
		MODE_ORBIT,
		MODE_FROZEN
	};

	enum AxisSelector {
		AXIS_V1,
		AXIS_V2,
		AXIS_V3
	};

	TwistConstraintForce();
	virtual ~TwistConstraintForce() {
	}

	std::tuple<std::vector<int>, std::string> init(input_file &inp) override;

	virtual LR_vector value(llint step, LR_vector &pos) override;
	virtual LR_vector torque(llint step, LR_vector &pos) override;
	virtual number potential(llint step, LR_vector &pos) override;

private:
	Mode _mode;
	AxisSelector _body_axis_selector;
	AxisSelector _body_ref_axis_selector;
	bool _body_ref_axis_auto;

	LR_vector _center;
	LR_vector _orbit_axis;
	LR_vector _move_dir;
	LR_vector _force_mask;

	LR_vector _target_axis0;
	LR_vector _target_ref0;
	bool _has_target_axis;
	bool _has_target_ref;

	number _orbit_rate;
	number _orbit_phase0;
	number _move_rate;
	number _axis_stiff;
	number _twist_stiff;
	number _twist_rate;
	number _twist_phase0;

	LR_vector _rotated_vector(const LR_vector &v, const LR_vector &axis, number angle) const;
	LR_vector _trap_position(llint step) const;
	LR_vector _target_axis(llint step) const;
	LR_vector _target_ref(llint step) const;
	LR_vector _selected_axis(AxisSelector selector) const;
};

#endif // TWISTCONSTRAINTFORCE_H_
