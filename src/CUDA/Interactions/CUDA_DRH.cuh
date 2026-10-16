/*
 * CUDA_DRH.cuh
 *
 *  Minimal CUDA helpers for the DRH interaction.
 */

#include "../../Interactions/drh_model.h"
#include "../cuda_utils/CUDA_lr_common.cuh"

__constant__ int MD_N[1];
__constant__ int MD_n_forces[1];
__constant__ int MD_DRH_dna_grooving[1];
__constant__ float MD_DRH_hb_multi[1];
__constant__ float MD_DRH_F1_A[2];
__constant__ float MD_DRH_F1_RC[2];
__constant__ float MD_DRH_F1_R0[2];
__constant__ float MD_DRH_F1_BLOW[2];
__constant__ float MD_DRH_F1_BHIGH[2];
__constant__ float MD_DRH_F1_RLOW[2];
__constant__ float MD_DRH_F1_RHIGH[2];
__constant__ float MD_DRH_F1_RCLOW[2];
__constant__ float MD_DRH_F1_RCHIGH[2];
__constant__ float MD_DRH_F1_EPS[50];
__constant__ float MD_DRH_F1_SHIFT[50];
__constant__ float MD_DRH_F2_K[2];
__constant__ float MD_DRH_F2_RC[2];
__constant__ float MD_DRH_F2_R0[2];
__constant__ float MD_DRH_F2_BLOW[2];
__constant__ float MD_DRH_F2_RLOW[2];
__constant__ float MD_DRH_F2_RCLOW[2];
__constant__ float MD_DRH_F2_BHIGH[2];
__constant__ float MD_DRH_F2_RCHIGH[2];
__constant__ float MD_DRH_F2_RHIGH[2];
__constant__ float MD_DRH_F4_THETA_A[13];
__constant__ float MD_DRH_F4_THETA_B[13];
__constant__ float MD_DRH_F4_THETA_T0[13];
__constant__ float MD_DRH_F4_THETA_TS[13];
__constant__ float MD_DRH_F4_THETA_TC[13];
__constant__ float MD_DRH_F5_PHI_A[4];
__constant__ float MD_DRH_F5_PHI_B[4];
__constant__ float MD_DRH_F5_PHI_XC[4];
__constant__ float MD_DRH_F5_PHI_XS[4];

struct DRH_RNAModel {
	float RNA_POS_BASE;
	float RNA_POS_STACK;
	float RNA_POS_BACK_a1;
	float RNA_POS_BACK_a2;
	float RNA_POS_BACK_a3;
};

__constant__ DRH_RNAModel MD_DRH_rnamodel[1];

__forceinline__ __device__ int _drh_f1_index(int type, int n3, int n5) {
	return 25 * type + n3 * 5 + n5;
}

__forceinline__ __device__ void _drh_map_pair_types(int &n3, int &n5, bool is_n3_dna, bool is_n5_dna) {
	bool au_bp = ((!is_n3_dna && n3 == 3) && (n5 == 0)) || ((!is_n5_dna && n5 == 3) && (n3 == 0));
	if(au_bp) {
		n3 = N_A;
		n5 = N_A;
	}

	if((n3 == 1 && n5 == 2) || (n5 == 1 && n3 == 2)) {
		if((!is_n3_dna && n3 == 1) || (!is_n5_dna && n5 == 1)) {
			n3 = N_G;
			n5 = N_C;
		}
		else {
			n3 = N_C;
			n5 = N_G;
		}
	}
}

__forceinline__ __device__ c_number _drh_f1(c_number r, int type, int n3, int n5, bool is_n3_dna, bool is_n5_dna) {
	_drh_map_pair_types(n3, n5, is_n3_dna, is_n5_dna);

	c_number val = (c_number) 0.f;
	if(r < MD_DRH_F1_RCHIGH[type]) {
		const int eps_index = _drh_f1_index(type, n3, n5);
		if(r > MD_DRH_F1_RHIGH[type]) {
			val = MD_DRH_F1_EPS[eps_index] * MD_DRH_F1_BHIGH[type] * SQR(r - MD_DRH_F1_RCHIGH[type]);
		}
		else if(r > MD_DRH_F1_RLOW[type]) {
			c_number tmp = 1.f - expf(-(r - MD_DRH_F1_R0[type]) * MD_DRH_F1_A[type]);
			val = MD_DRH_F1_EPS[eps_index] * SQR(tmp) - MD_DRH_F1_SHIFT[eps_index];
		}
		else if(r > MD_DRH_F1_RCLOW[type]) {
			val = MD_DRH_F1_EPS[eps_index] * MD_DRH_F1_BLOW[type] * SQR(r - MD_DRH_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _drh_f1D(c_number r, int type, int n3, int n5, bool is_n3_dna, bool is_n5_dna) {
	_drh_map_pair_types(n3, n5, is_n3_dna, is_n5_dna);

	c_number val = (c_number) 0.f;
	if(r < MD_DRH_F1_RCHIGH[type]) {
		const float eps = MD_DRH_F1_EPS[_drh_f1_index(type, n3, n5)];
		if(r > MD_DRH_F1_RHIGH[type]) {
			val = 2.f * eps * MD_DRH_F1_BHIGH[type] * (r - MD_DRH_F1_RCHIGH[type]);
		}
		else if(r > MD_DRH_F1_RLOW[type]) {
			c_number tmp = expf(-(r - MD_DRH_F1_R0[type]) * MD_DRH_F1_A[type]);
			val = 2.f * eps * (1.f - tmp) * tmp * MD_DRH_F1_A[type];
		}
		else if(r > MD_DRH_F1_RCLOW[type]) {
			val = 2.f * eps * MD_DRH_F1_BLOW[type] * (r - MD_DRH_F1_RCLOW[type]);
		}
	}

	return val;
}

__forceinline__ __device__ c_number _drh_f4(c_number t, int type) {
	c_number val = (c_number) 0.f;
	t = copysignf(t - MD_DRH_F4_THETA_T0[type], (c_number) 1.f);
	if(t < MD_DRH_F4_THETA_TC[type]) {
		val = (t > MD_DRH_F4_THETA_TS[type]) ? MD_DRH_F4_THETA_B[type] * SQR(MD_DRH_F4_THETA_TC[type] - t) : 1.f - MD_DRH_F4_THETA_A[type] * SQR(t);
	}
	return val;
}

__forceinline__ __device__ c_number _drh_f4D(c_number t, int type) {
	c_number val = (c_number) 0.f;
	t -= MD_DRH_F4_THETA_T0[type];
	c_number m = copysignf((c_number) 1.f, t);
	t = copysignf(t, (c_number) 1.f);
	if(t < MD_DRH_F4_THETA_TC[type]) {
		val = (t > MD_DRH_F4_THETA_TS[type]) ? 2.f * m * MD_DRH_F4_THETA_B[type] * (t - MD_DRH_F4_THETA_TC[type]) : -2.f * m * MD_DRH_F4_THETA_A[type] * t;
	}
	return val;
}

__forceinline__ __device__ c_number _drh_f2(c_number r, int type) {
	c_number val = (c_number) 0.f;
	if(r < MD_DRH_F2_RCHIGH[type]) {
		if(r > MD_DRH_F2_RHIGH[type]) {
			val = MD_DRH_F2_K[type] * MD_DRH_F2_BHIGH[type] * SQR(r - MD_DRH_F2_RCHIGH[type]);
		}
		else if(r > MD_DRH_F2_RLOW[type]) {
			val = (MD_DRH_F2_K[type] * 0.5f) * (SQR(r - MD_DRH_F2_R0[type]) - SQR(MD_DRH_F2_RC[type] - MD_DRH_F2_R0[type]));
		}
		else if(r > MD_DRH_F2_RCLOW[type]) {
			val = MD_DRH_F2_K[type] * MD_DRH_F2_BLOW[type] * SQR(r - MD_DRH_F2_RCLOW[type]);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _drh_f2D(c_number r, int type) {
	c_number val = (c_number) 0.f;
	if(r < MD_DRH_F2_RCHIGH[type]) {
		if(r > MD_DRH_F2_RHIGH[type]) {
			val = 2.f * MD_DRH_F2_K[type] * MD_DRH_F2_BHIGH[type] * (r - MD_DRH_F2_RCHIGH[type]);
		}
		else if(r > MD_DRH_F2_RLOW[type]) {
			val = MD_DRH_F2_K[type] * (r - MD_DRH_F2_R0[type]);
		}
		else if(r > MD_DRH_F2_RCLOW[type]) {
			val = 2.f * MD_DRH_F2_K[type] * MD_DRH_F2_BLOW[type] * (r - MD_DRH_F2_RCLOW[type]);
		}
	}
	return val;
}

__forceinline__ __device__ c_number _drh_f5(c_number f, int type) {
	c_number val = (c_number) 0.f;
	if(f > MD_DRH_F5_PHI_XC[type]) {
		if(f < MD_DRH_F5_PHI_XS[type]) {
			val = MD_DRH_F5_PHI_B[type] * SQR(MD_DRH_F5_PHI_XC[type] - f);
		}
		else if(f < 0.f) {
			val = 1.f - MD_DRH_F5_PHI_A[type] * SQR(f);
		}
		else {
			val = 1.f;
		}
	}
	return val;
}

__forceinline__ __device__ c_number _drh_f5D(c_number f, int type) {
	c_number val = (c_number) 0.f;
	if(f > MD_DRH_F5_PHI_XC[type]) {
		if(f < MD_DRH_F5_PHI_XS[type]) {
			val = 2.f * MD_DRH_F5_PHI_B[type] * (f - MD_DRH_F5_PHI_XC[type]);
		}
		else if(f < 0.f) {
			val = -2.f * MD_DRH_F5_PHI_A[type] * f;
		}
	}
	return val;
}

__forceinline__ __device__ c_number4 _drh_base_site(const c_number4 &a1, bool is_dna) {
	return a1 * (is_dna ? POS_BASE : MD_DRH_rnamodel[0].RNA_POS_BASE);
}

__forceinline__ __device__ c_number4 _drh_stack_site(const c_number4 &a1, bool is_dna) {
	return a1 * (is_dna ? POS_STACK : MD_DRH_rnamodel[0].RNA_POS_STACK);
}

__forceinline__ __device__ c_number4 _drh_back_site(const c_number4 &a1, const c_number4 &a2, const c_number4 &a3, bool is_dna) {
	if(is_dna) {
		return (MD_DRH_dna_grooving[0] != 0) ? POS_MM_BACK1 * a1 + POS_MM_BACK2 * a2 : POS_BACK * a1;
	}

	return MD_DRH_rnamodel[0].RNA_POS_BACK_a1 * a1 + MD_DRH_rnamodel[0].RNA_POS_BACK_a2 * a2 + MD_DRH_rnamodel[0].RNA_POS_BACK_a3 * a3;
}

__forceinline__ __device__ void _drh_excluded_volume(const c_number4 &r, c_number4 &F, c_number sigma, c_number rstar, c_number b, c_number rc) {
	const c_number rsqr = CUDA_DOT(r, r);
	F.x = F.y = F.z = F.w = (c_number) 0.f;

	if(rsqr < SQR(rc)) {
		if(rsqr > SQR(rstar)) {
			const c_number rmod = sqrtf(rsqr);
			const c_number rrc = rmod - rc;
			const c_number fmod = 2.f * DRH_EXCL_EPS * b * rrc / rmod;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = DRH_EXCL_EPS * b * SQR(rrc);
		}
		else {
			const c_number lj_part = CUB(SQR(sigma) / rsqr);
			const c_number fmod = 24.f * DRH_EXCL_EPS * (lj_part - 2.f * SQR(lj_part)) / rsqr;
			F.x = r.x * fmod;
			F.y = r.y * fmod;
			F.z = r.z * fmod;
			F.w = 4.f * DRH_EXCL_EPS * (SQR(lj_part) - lj_part);
		}
	}
}

__device__ void _particle_particle_DRH_hb_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1, const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, bool p_is_dna, bool q_is_dna) {
	(void) ppos;
	(void) a2;
	(void) qpos;
	(void) b2;

	const int ptype = get_particle_type(ppos);
	const int qtype = get_particle_type(qpos);
	const int pbtype = get_particle_btype(ppos);
	const int qbtype = get_particle_btype(qpos);
	if(pbtype + qbtype != 3) {
		return;
	}

	const c_number hb_multi = (abs(qbtype) >= 300 && abs(pbtype) >= 300) ? MD_DRH_hb_multi[0] : 1.f;
	const c_number4 ppos_base = _drh_base_site(a1, p_is_dna);
	const c_number4 qpos_base = _drh_base_site(b1, q_is_dna);
	const c_number4 rhydro = r + qpos_base - ppos_base;
	const c_number rhydromodsqr = CUDA_DOT(rhydro, rhydro);
	if(rhydromodsqr <= SQR(DRH_HYDR_RCLOW) || rhydromodsqr >= SQR(DRH_HYDR_RCHIGH)) {
		return;
	}

	const c_number rhydromod = sqrtf(rhydromodsqr);
	const c_number4 rhydrodir = rhydro / rhydromod;
	const c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
	const c_number cost2 = -CUDA_DOT(b1, rhydrodir);
	const c_number t2 = CUDA_LRACOS(cost2);
	const c_number cost3 = CUDA_DOT(a1, rhydrodir);
	const c_number t3 = CUDA_LRACOS(cost3);
	const c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
	const c_number cost7 = -CUDA_DOT(rhydrodir, b3);
	const c_number t7 = CUDA_LRACOS(cost7);
	const c_number cost8 = CUDA_DOT(rhydrodir, a3);
	const c_number t8 = CUDA_LRACOS(cost8);

	const c_number f1 = hb_multi * _drh_f1(rhydromod, DRH_HYDR_F1, qtype, ptype, q_is_dna, p_is_dna);
	const c_number f4t1 = _drh_f4(t1, DRH_HYDR_F4_THETA1);
	const c_number f4t2 = _drh_f4(t2, DRH_HYDR_F4_THETA2);
	const c_number f4t3 = _drh_f4(t3, DRH_HYDR_F4_THETA3);
	const c_number f4t4 = _drh_f4(t4, DRH_HYDR_F4_THETA4);
	const c_number f4t7 = _drh_f4(t7, DRH_HYDR_F4_THETA7);
	const c_number f4t8 = _drh_f4(t8, DRH_HYDR_F4_THETA8);
	const c_number hb_energy = f1 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	if(hb_energy == (c_number) 0.f) {
		return;
	}

	c_number4 Ftmp = rhydrodir * (hb_energy * _drh_f1D(rhydromod, DRH_HYDR_F1, qtype, ptype, q_is_dna, p_is_dna) / f1);
	c_number4 Ttmp = make_c_number4(0.f, 0.f, 0.f, 0.f);

	const c_number f4t1D = -_drh_f4D(t1, DRH_HYDR_F4_THETA1);
	const c_number f4t2D = -_drh_f4D(t2, DRH_HYDR_F4_THETA2);
	const c_number f4t3D = _drh_f4D(t3, DRH_HYDR_F4_THETA3);
	const c_number f4t4D = _drh_f4D(t4, DRH_HYDR_F4_THETA4);
	const c_number f4t7D = -_drh_f4D(t7, DRH_HYDR_F4_THETA7);
	const c_number f4t8D = _drh_f4D(t8, DRH_HYDR_F4_THETA8);

	Ttmp -= stably_normalised(_cross(a3, b3)) * (-hb_energy * f4t4D / f4t4);
	Ttmp -= stably_normalised(_cross(a1, b1)) * (-hb_energy * f4t1D / f4t1);
	Ftmp -= stably_normalised(b1 + rhydrodir * cost2) * (hb_energy * f4t2D / (f4t2 * rhydromod));

	c_number part = -hb_energy * f4t3D / f4t3;
	Ftmp -= stably_normalised(a1 - rhydrodir * cost3) * (-part / rhydromod);
	Ttmp += stably_normalised(_cross(rhydrodir, a1)) * part;

	Ftmp -= stably_normalised(b3 + rhydrodir * cost7) * (hb_energy * f4t7D / (f4t7 * rhydromod));

	part = -hb_energy * f4t8D / f4t8;
	Ftmp -= stably_normalised(a3 - rhydrodir * cost8) * (-part / rhydromod);
	Ttmp += stably_normalised(_cross(rhydrodir, a3)) * part;

	Ttmp += _cross(ppos_base, Ftmp);
	Ftmp.w = hb_energy;
	F += Ftmp;
	T += _vectors_transpose_c_number4_product(a1, a2, a3, Ttmp);
	T.w += hb_energy;
}

__device__ void _particle_particle_DRH_excl_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1, const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, bool p_is_dna, bool q_is_dna) {
	(void) ppos;
	(void) qpos;

	const c_number4 ppos_base = _drh_base_site(a1, p_is_dna);
	const c_number4 qpos_base = _drh_base_site(b1, q_is_dna);
	const c_number4 ppos_back = _drh_back_site(a1, a2, a3, p_is_dna);
	const c_number4 qpos_back = _drh_back_site(b1, b2, b3, q_is_dna);

	c_number4 pair_force = make_c_number4(0.f, 0.f, 0.f, 0.f);
	c_number4 pair_torque = make_c_number4(0.f, 0.f, 0.f, 0.f);
	c_number4 rcenter = r + qpos_base - ppos_base;
	_drh_excluded_volume(rcenter, pair_force, DRH_EXCL_S2, DRH_EXCL_R2, DRH_EXCL_B2, DRH_EXCL_RC2);
	pair_torque += _cross(ppos_base, pair_force);
	F += pair_force;

	rcenter = r + qpos_base - ppos_back;
	_drh_excluded_volume(rcenter, pair_force, DRH_EXCL_S4, DRH_EXCL_R4, DRH_EXCL_B4, DRH_EXCL_RC4);
	pair_torque += _cross(ppos_back, pair_force);
	F += pair_force;

	rcenter = r + qpos_back - ppos_base;
	_drh_excluded_volume(rcenter, pair_force, DRH_EXCL_S3, DRH_EXCL_R3, DRH_EXCL_B3, DRH_EXCL_RC3);
	pair_torque += _cross(ppos_base, pair_force);
	F += pair_force;

	rcenter = r + qpos_back - ppos_back;
	_drh_excluded_volume(rcenter, pair_force, DRH_EXCL_S1, DRH_EXCL_R1, DRH_EXCL_B1, DRH_EXCL_RC1);
	pair_torque += _cross(ppos_back, pair_force);
	F += pair_force;

	T += _vectors_transpose_c_number4_product(a1, a2, a3, pair_torque);
}

__device__ void _particle_particle_DRH_cross_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1, const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, bool p_is_dna, bool q_is_dna) {
	(void) a2;
	(void) b2;

	const c_number4 ppos_base = _drh_base_site(a1, p_is_dna);
	const c_number4 qpos_base = _drh_base_site(b1, q_is_dna);
	const c_number4 rcstack = r + qpos_base - ppos_base;
	const c_number rcstackmodsqr = CUDA_DOT(rcstack, rcstack);
	if(rcstackmodsqr <= SQR(DRH_CRST_RCLOW) || rcstackmodsqr >= SQR(DRH_CRST_RCHIGH)) {
		return;
	}

	const c_number rcstackmod = sqrtf(rcstackmodsqr);
	const c_number4 rcstackdir = rcstack / rcstackmod;
	const c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
	const c_number cost2 = -CUDA_DOT(b1, rcstackdir);
	const c_number t2 = CUDA_LRACOS(cost2);
	const c_number cost3 = CUDA_DOT(a1, rcstackdir);
	const c_number t3 = CUDA_LRACOS(cost3);
	const c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
	const c_number cost7 = -CUDA_DOT(rcstackdir, b3);
	const c_number t7 = CUDA_LRACOS(cost7);
	const c_number cost8 = CUDA_DOT(rcstackdir, a3);
	const c_number t8 = CUDA_LRACOS(cost8);

	const c_number f2 = _drh_f2(rcstackmod, DRH_CRST_F2);
	const c_number f4t1 = _drh_f4(t1, DRH_CRST_F4_THETA1);
	const c_number f4t2 = _drh_f4(t2, DRH_CRST_F4_THETA2);
	const c_number f4t3 = _drh_f4(t3, DRH_CRST_F4_THETA3);
	const c_number f4t4 = _drh_f4(t4, DRH_CRST_F4_THETA4) + _drh_f4(PI - t4, DRH_CRST_F4_THETA4);
	const c_number f4t7 = _drh_f4(t7, DRH_CRST_F4_THETA7) + _drh_f4(PI - t7, DRH_CRST_F4_THETA7);
	const c_number f4t8 = _drh_f4(t8, DRH_CRST_F4_THETA8) + _drh_f4(PI - t8, DRH_CRST_F4_THETA8);
	const c_number cstk_energy = f2 * f4t1 * f4t2 * f4t3 * f4t4 * f4t7 * f4t8;
	if(cstk_energy == (c_number) 0.f) {
		return;
	}

	const c_number f2D = _drh_f2D(rcstackmod, DRH_CRST_F2);

	const c_number f4t1D = -_drh_f4D(t1, DRH_CRST_F4_THETA1);
	const c_number f4t2D = -_drh_f4D(t2, DRH_CRST_F4_THETA2);
	const c_number f4t3D = _drh_f4D(t3, DRH_CRST_F4_THETA3);
	const c_number f4t4D = _drh_f4D(t4, DRH_CRST_F4_THETA4) - _drh_f4D(PI - t4, DRH_CRST_F4_THETA4);
	const c_number f4t7D = -_drh_f4D(t7, DRH_CRST_F4_THETA7) + _drh_f4D(PI - t7, DRH_CRST_F4_THETA7);
	const c_number f4t8D = _drh_f4D(t8, DRH_CRST_F4_THETA8) - _drh_f4D(PI - t8, DRH_CRST_F4_THETA8);

	c_number4 Ftmp = rcstackdir * (cstk_energy * f2D / f2);
	c_number4 Ttmp = make_c_number4(0.f, 0.f, 0.f, 0.f);

	Ttmp -= stably_normalised(_cross(a1, b1)) * (-cstk_energy * f4t1D / f4t1);
	Ftmp -= stably_normalised(b1 + rcstackdir * cost2) * (cstk_energy * f4t2D / (f4t2 * rcstackmod));

	c_number part = -cstk_energy * f4t3D / f4t3;
	Ftmp -= stably_normalised(a1 - rcstackdir * cost3) * (-part / rcstackmod);
	Ttmp += stably_normalised(_cross(rcstackdir, a1)) * part;

	Ttmp -= stably_normalised(_cross(a3, b3)) * (-cstk_energy * f4t4D / f4t4);
	Ftmp -= stably_normalised(b3 + rcstackdir * cost7) * (cstk_energy * f4t7D / (f4t7 * rcstackmod));

	part = -cstk_energy * f4t8D / f4t8;
	Ftmp -= stably_normalised(a3 - rcstackdir * cost8) * (-part / rcstackmod);
	Ttmp += stably_normalised(_cross(rcstackdir, a3)) * part;

	Ttmp += _cross(ppos_base, Ftmp);
	Ftmp.w = cstk_energy;
	F += Ftmp;
	T += _vectors_transpose_c_number4_product(a1, a2, a3, Ttmp);
}

__device__ void _particle_particle_DRH_coax_interaction(const c_number4 &r, const c_number4 &ppos, const c_number4 &a1, const c_number4 &a2, const c_number4 &a3,
		const c_number4 &qpos, const c_number4 &b1, const c_number4 &b2, const c_number4 &b3, c_number4 &F, c_number4 &T, bool p_is_dna, bool q_is_dna) {
	(void) ppos;
	(void) qpos;

	const c_number4 ppos_stack = _drh_stack_site(a1, p_is_dna);
	const c_number4 qpos_stack = _drh_stack_site(b1, q_is_dna);
	const c_number4 ppos_back = _drh_back_site(a1, a2, a3, p_is_dna);
	const c_number4 qpos_back = _drh_back_site(b1, b2, b3, q_is_dna);

	const c_number4 rstack = r + qpos_stack - ppos_stack;
	const c_number rstackmodsqr = CUDA_DOT(rstack, rstack);
	if(rstackmodsqr <= SQR(DRH_CXST_RCLOW) || rstackmodsqr >= SQR(DRH_CXST_RCHIGH)) {
		return;
	}

	const c_number rstackmod = sqrtf(rstackmodsqr);
	const c_number4 rstackdir = rstack / rstackmod;
	const c_number cost5 = CUDA_DOT(a3, rstackdir);
	const c_number cost6 = -CUDA_DOT(b3, rstackdir);
	const c_number t1 = CUDA_LRACOS(-CUDA_DOT(a1, b1));
	const c_number t4 = CUDA_LRACOS(CUDA_DOT(a3, b3));
	const c_number t5 = CUDA_LRACOS(cost5);
	const c_number t6 = CUDA_LRACOS(cost6);

	const c_number4 rbackbone = r + qpos_back - ppos_back;
	const c_number rbackmod = _module(rbackbone);
	const c_number4 rbackbonedir = rbackbone / rbackmod;
	const c_number cosphi3 = CUDA_DOT(rstackdir, _cross(rbackbonedir, a1));
	const c_number cosphi4 = CUDA_DOT(rstackdir, _cross(rbackbonedir, b1));

	const c_number f2 = _drh_f2(rstackmod, DRH_CXST_F2);
	const c_number f4t1 = _drh_f4(t1, DRH_CXST_F4_THETA1) + _drh_f4(2.f * PI - t1, DRH_CXST_F4_THETA1);
	const c_number f4t4 = _drh_f4(t4, DRH_CXST_F4_THETA4);
	const c_number f4t5 = _drh_f4(t5, DRH_CXST_F4_THETA5) + _drh_f4(PI - t5, DRH_CXST_F4_THETA5);
	const c_number f4t6 = _drh_f4(t6, DRH_CXST_F4_THETA6) + _drh_f4(PI - t6, DRH_CXST_F4_THETA6);
	const c_number f5cosphi3 = _drh_f5(cosphi3, DRH_CXST_F5_PHI3);
	const c_number f5cosphi4 = _drh_f5(cosphi4, DRH_CXST_F5_PHI4);
	const c_number energy = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
	if(energy == (c_number) 0.f) {
		return;
	}

	c_number4 force = -rstackdir * _drh_f2D(rstackmod, DRH_CXST_F2) * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
	c_number4 torquep = make_c_number4(0.f, 0.f, 0.f, 0.f);

	const c_number f4t1Dsin = -_drh_f4D(t1, DRH_CXST_F4_THETA1) + _drh_f4D(2.f * PI - t1, DRH_CXST_F4_THETA1);
	const c_number f4t4Dsin = _drh_f4D(t4, DRH_CXST_F4_THETA4);
	const c_number f4t5Dsin = _drh_f4D(t5, DRH_CXST_F4_THETA5) - _drh_f4D(PI - t5, DRH_CXST_F4_THETA5);
	const c_number f4t6Dsin = -_drh_f4D(t6, DRH_CXST_F4_THETA6) + _drh_f4D(PI - t6, DRH_CXST_F4_THETA6);
	const c_number f5Dcosphi3 = _drh_f5D(cosphi3, DRH_CXST_F5_PHI3);
	const c_number f5Dcosphi4 = _drh_f5D(cosphi4, DRH_CXST_F5_PHI4);

	c_number4 dir = _cross(a1, b1);
	c_number torquemod = -f2 * f4t1Dsin * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
	torquep -= dir * torquemod;

	dir = _cross(a3, b3);
	torquemod = -f2 * f4t1 * f4t4Dsin * f4t5 * f4t6 * f5cosphi3 * f5cosphi4;
	torquep -= dir * torquemod;

	c_number fact = f2 * f4t1 * f4t4 * f4t5Dsin * f4t6 * f5cosphi3 * f5cosphi4;
	force += (a3 - rstackdir * cost5) * (fact / rstackmod);
	dir = _cross(rstackdir, a3);
	torquep -= dir * fact;

	fact = f2 * f4t1 * f4t4 * f4t5 * f4t6Dsin * f5cosphi3 * f5cosphi4;
	force += (b3 + rstackdir * cost6) * (fact / rstackmod);

	torquep -= _cross(ppos_stack, force);

	c_number4 myforce = make_c_number4(0.f, 0.f, 0.f, 0.f);
	c_number4 mytorque1 = make_c_number4(0.f, 0.f, 0.f, 0.f);

	c_number force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi4 * f5Dcosphi3;
	c_number4 rba1 = _cross(rbackbonedir, a1);
	c_number4 a1rs = _cross(a1, rstackdir);
	c_number rb_dot_a1rs = CUDA_DOT(rbackbonedir, a1rs);
	c_number rs_dot_rba1 = CUDA_DOT(rstackdir, rba1);

	c_number4 forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
	c_number4 forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod);
	myforce += forcestack + forceback;
	mytorque1 -= _cross(ppos_stack, forcestack);
	mytorque1 -= _cross(ppos_back, forceback);
	mytorque1 -= force_c * _cross(a1, _cross(rstackdir, rbackbonedir));

	force_c = f2 * f4t1 * f4t4 * f4t5 * f4t6 * f5cosphi3 * f5Dcosphi4;
	rba1 = _cross(rbackbonedir, b1);
	a1rs = _cross(b1, rstackdir);
	rb_dot_a1rs = CUDA_DOT(rbackbonedir, a1rs);
	rs_dot_rba1 = CUDA_DOT(rstackdir, rba1);
	forcestack = -(rba1 - rstackdir * rs_dot_rba1) * (force_c / rstackmod);
	forceback = -(a1rs - rbackbonedir * rb_dot_a1rs) * (force_c / rbackmod);
	myforce += forcestack + forceback;

	if(p_is_dna) {
		mytorque1 -= _cross(ppos_stack, forcestack);
		mytorque1 -= _cross(ppos_back, forceback);
	}
	else {
		mytorque1 -= _cross(ppos_stack, forcestack);
		mytorque1 -= _cross(ppos_back, forceback);
	}

	force += myforce;
	torquep += mytorque1;

	c_number4 Ftmp = make_c_number4(-force.x, -force.y, -force.z, energy);
	F += Ftmp;
	T += _vectors_transpose_c_number4_product(a1, a2, a3, torquep);
}

__global__ void drh_forces_hybrid_hb(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs, const int *number_neighs, const LR_bonds __restrict__ *bonds, const uint8_t *is_dna,
		const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	const c_number4 ppos = poss[IND];
	const LR_bonds pbonds = bonds[IND];
	const bool p_is_dna = (is_dna[IND] != 0);

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	const int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);
		if(k_index == IND || pbonds.n3 == k_index || pbonds.n5 == k_index) {
			continue;
		}

		const bool q_is_dna = (is_dna[k_index] != 0);
		if(q_is_dna == p_is_dna) {
			continue;
		}

		const c_number4 qpos = poss[k_index];
		const c_number4 r = box->minimum_image(ppos, qpos);
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[k_index], b1, b2, b3);
		_particle_particle_DRH_hb_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, p_is_dna, q_is_dna);
	}

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void drh_forces_hybrid_excl(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs, const int *number_neighs, const LR_bonds __restrict__ *bonds, const uint8_t *is_dna,
		const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	const c_number4 ppos = poss[IND];
	const LR_bonds pbonds = bonds[IND];
	const bool p_is_dna = (is_dna[IND] != 0);

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	const int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);
		if(k_index == IND || pbonds.n3 == k_index || pbonds.n5 == k_index) {
			continue;
		}

		const bool q_is_dna = (is_dna[k_index] != 0);
		if(q_is_dna == p_is_dna) {
			continue;
		}

		const c_number4 qpos = poss[k_index];
		const c_number4 r = box->minimum_image(ppos, qpos);
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[k_index], b1, b2, b3);
		_particle_particle_DRH_excl_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, p_is_dna, q_is_dna);
	}

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void drh_forces_hybrid_cross(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs, const int *number_neighs, const LR_bonds __restrict__ *bonds, const uint8_t *is_dna,
		const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	const c_number4 ppos = poss[IND];
	const LR_bonds pbonds = bonds[IND];
	const bool p_is_dna = (is_dna[IND] != 0);

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	const int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);
		if(k_index == IND || pbonds.n3 == k_index || pbonds.n5 == k_index) {
			continue;
		}

		const bool q_is_dna = (is_dna[k_index] != 0);
		if(q_is_dna == p_is_dna) {
			continue;
		}

		const c_number4 qpos = poss[k_index];
		const c_number4 r = box->minimum_image(ppos, qpos);
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[k_index], b1, b2, b3);
		_particle_particle_DRH_cross_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, p_is_dna, q_is_dna);
	}

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void drh_forces_hybrid_coax(const c_number4 __restrict__ *poss, const GPU_quat __restrict__ *orientations, c_number4 __restrict__ *forces,
		c_number4 __restrict__ *torques, const int *matrix_neighs, const int *number_neighs, const LR_bonds __restrict__ *bonds, const uint8_t *is_dna,
		const CUDABox *box) {
	if(IND >= MD_N[0]) return;

	c_number4 F = forces[IND];
	c_number4 T = torques[IND];
	const c_number4 ppos = poss[IND];
	const LR_bonds pbonds = bonds[IND];
	const bool p_is_dna = (is_dna[IND] != 0);

	c_number4 a1, a2, a3;
	get_vectors_from_quat(orientations[IND], a1, a2, a3);

	const int num_neighs = NUMBER_NEIGHBOURS(IND, number_neighs);
	for(int j = 0; j < num_neighs; j++) {
		const int k_index = NEXT_NEIGHBOUR(IND, j, matrix_neighs);
		if(k_index == IND || pbonds.n3 == k_index || pbonds.n5 == k_index) {
			continue;
		}

		const bool q_is_dna = (is_dna[k_index] != 0);
		if(q_is_dna == p_is_dna) {
			continue;
		}

		const c_number4 qpos = poss[k_index];
		const c_number4 r = box->minimum_image(ppos, qpos);
		c_number4 b1, b2, b3;
		get_vectors_from_quat(orientations[k_index], b1, b2, b3);
		_particle_particle_DRH_coax_interaction(r, ppos, a1, a2, a3, qpos, b1, b2, b3, F, T, p_is_dna, q_is_dna);
	}

	forces[IND] = F;
	torques[IND] = T;
}

__global__ void init_DRH_strand_ends(int *is_strand_end, const LR_bonds __restrict__ *bonds, int N) {
	if(IND >= N) return;

	LR_bonds pbonds = bonds[IND];
	is_strand_end[IND] = (pbonds.n3 == P_INVALID || pbonds.n5 == P_INVALID);
}
