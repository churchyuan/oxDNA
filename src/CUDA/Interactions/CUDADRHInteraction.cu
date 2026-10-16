/*
 * CUDADRHInteraction.cu
 *
 *  Created on: 26/apr/2026
 */

#include "CUDADRHInteraction.h"

#include "CUDA_DRH.cuh"
#include "../../Particles/DNANucleotide.h"
#include "../../Particles/RNANucleotide.h"
#include "../../Utilities/ConfigInfo.h"
#include "../../Utilities/oxDNAException.h"

namespace {
void copy_Model_to_DRHModel(Model &model_from, DRH_RNAModel &model_to) {
	model_to.RNA_POS_BASE = model_from.RNA_POS_BASE;
	model_to.RNA_POS_STACK = model_from.RNA_POS_STACK;
	model_to.RNA_POS_BACK_a1 = model_from.RNA_POS_BACK_a1;
	model_to.RNA_POS_BACK_a2 = model_from.RNA_POS_BACK_a2;
	model_to.RNA_POS_BACK_a3 = model_from.RNA_POS_BACK_a3;
}
}

CUDADRHInteraction::CUDADRHInteraction() {
	_edge_compatible = false;
}

CUDADRHInteraction::~CUDADRHInteraction() {
	if(_d_is_dna != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_is_dna));
	}
	if(_d_is_strand_end != nullptr) {
		CUDA_SAFE_CALL(cudaFree(_d_is_strand_end));
	}
}

void CUDADRHInteraction::get_settings(input_file &inp) {
	Logger::instance()->disable_log("CUDADRHInteraction");
	DRHInteraction::get_settings(inp);
	_dna_helper.get_settings(inp);
	_rna_helper.get_settings(inp);
	Logger::instance()->enable_log("CUDADRHInteraction");
}

void CUDADRHInteraction::get_cuda_settings(input_file &inp) {
	CUDABaseInteraction::get_cuda_settings(inp);
	_dna_helper.get_cuda_settings(inp);
	_rna_helper.get_cuda_settings(inp);
}

void CUDADRHInteraction::_init_acid_types() {
	if(_d_is_dna != nullptr) {
		return;
	}

	std::vector<uint8_t> is_dna(_N, 0);
	for(int i = 0; i < _N; i++) {
		BaseParticle *particle = CONFIG_INFO->particles()[i];
		if(dynamic_cast<DNANucleotide *>(particle) != nullptr) {
			is_dna[i] = 1;
		}
		else if(dynamic_cast<RNANucleotide *>(particle) != nullptr) {
			is_dna[i] = 0;
		}
		else {
			throw oxDNAException("CUDADRHInteraction expected DNA/RNA particles only, found an unsupported particle type at index %d", i);
		}
	}

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<uint8_t>(&_d_is_dna, sizeof(uint8_t) * _N));
	CUDA_SAFE_CALL(cudaMemcpy(_d_is_dna, is_dna.data(), sizeof(uint8_t) * _N, cudaMemcpyHostToDevice));
}

void CUDADRHInteraction::_init_strand_ends(LR_bonds *d_bonds) {
	if(_d_is_strand_end != nullptr) {
		return;
	}

	CUDA_SAFE_CALL(GpuUtils::LR_cudaMalloc<int>(&_d_is_strand_end, sizeof(int) * _N));
	init_DRH_strand_ends<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(_d_is_strand_end, d_bonds, _N);
	CUT_CHECK_ERROR("init_DRH_strand_ends error");
}

void CUDADRHInteraction::cuda_init(int N) {
	CUDABaseInteraction::cuda_init(N);
	Logger::instance()->disable_log("CUDADRHInteraction");
	DRHInteraction::init();
	Logger::instance()->enable_log("CUDADRHInteraction");

	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_N, &N, sizeof(int)));
	if(this->_use_edge) {
		CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_n_forces, &this->_n_forces, sizeof(int)));
	}
	int grooving = DNAInteraction::_grooving ? 1 : 0;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_DRH_dna_grooving, &grooving, sizeof(int)));
	float hb_multi = this->_hb_multiplier;
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_DRH_hb_multi, &hb_multi, sizeof(float)));

	c_number tmp[50];
	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				tmp[i * 25 + j * 5 + k] = this->F1_EPS[i][j][k];
			}
		}
	}
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_EPS, tmp, 50);
	for(int i = 0; i < 2; i++) {
		for(int j = 0; j < 5; j++) {
			for(int k = 0; k < 5; k++) {
				tmp[i * 25 + j * 5 + k] = this->F1_SHIFT[i][j][k];
			}
		}
	}
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_SHIFT, tmp, 50);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_A, this->F1_A, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_RC, this->F1_RC, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_R0, this->F1_R0, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_BLOW, this->F1_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_BHIGH, this->F1_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_RLOW, this->F1_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_RHIGH, this->F1_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_RCLOW, this->F1_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F1_RCHIGH, this->F1_RCHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_K, this->F2_K, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_RC, this->F2_RC, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_R0, this->F2_R0, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_BLOW, this->F2_BLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_RLOW, this->F2_RLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_RCLOW, this->F2_RCLOW, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_BHIGH, this->F2_BHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_RCHIGH, this->F2_RCHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F2_RHIGH, this->F2_RHIGH, 2);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F4_THETA_A, this->F4_THETA_A, 13);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F4_THETA_B, this->F4_THETA_B, 13);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F4_THETA_T0, this->F4_THETA_T0, 13);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F4_THETA_TS, this->F4_THETA_TS, 13);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F4_THETA_TC, this->F4_THETA_TC, 13);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F5_PHI_A, this->F5_PHI_A, 4);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F5_PHI_B, this->F5_PHI_B, 4);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F5_PHI_XC, this->F5_PHI_XC, 4);
	COPY_ARRAY_TO_CONSTANT(MD_DRH_F5_PHI_XS, this->F5_PHI_XS, 4);

	DRH_RNAModel drh_model;
	copy_Model_to_DRHModel(*(model), drh_model);
	CUDA_SAFE_CALL(cudaMemcpyToSymbol(MD_DRH_rnamodel, &drh_model, sizeof(DRH_RNAModel)));

	_init_acid_types();

	_dna_helper._use_debye_huckel = true;
	_dna_helper._use_oxDNA2_coaxial_stacking = true;
	_dna_helper._use_oxDNA2_FENE = true;
	_dna_helper._salt_concentration = DNA2Interaction::_salt_concentration;
	_dna_helper._debye_huckel_half_charged_ends = DNA2Interaction::_debye_huckel_half_charged_ends;
	_dna_helper._debye_huckel_prefactor = DNA2Interaction::_debye_huckel_prefactor;
	_dna_helper._debye_huckel_lambdafactor = DNA2Interaction::_debye_huckel_lambdafactor;
	_dna_helper.set_grooving(DNAInteraction::_grooving);
	_dna_helper.F2_K[1] = CXST_K_OXDNA2;
	_dna_helper.set_launch_cfg(_launch_cfg);
	_dna_helper.set_particle_filter(_d_is_dna, true);
	_dna_helper.cuda_init(N);

	_rna_helper._use_debye_huckel = true;
	_rna_helper._mismatch_repulsion = RNA2Interaction::_mismatch_repulsion;
	_rna_helper._salt_concentration = RNA2Interaction::_salt_concentration;
	_rna_helper._debye_huckel_half_charged_ends = RNA2Interaction::_debye_huckel_half_charged_ends;
	_rna_helper._debye_huckel_prefactor = RNA2Interaction::_debye_huckel_prefactor;
	_rna_helper._debye_huckel_lambdafactor = RNA2Interaction::_debye_huckel_lambdafactor;
	_rna_helper._RNA_HYDR_MIS = RNA2Interaction::_RNA_HYDR_MIS;
	_rna_helper._grooving = true;
	_rna_helper.set_launch_cfg(_launch_cfg);
	_rna_helper.set_particle_filter(_d_is_dna, false);
	_rna_helper.cuda_init(N);
}

void CUDADRHInteraction::_on_T_update() {
	cuda_init(_N);
}

void CUDADRHInteraction::compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) {
	_init_acid_types();
	_init_strand_ends(d_bonds);

	if(_update_st) {
		throw oxDNAException("CUDA stress tensor updates are not implemented yet for interaction_type = NA");
	}
	if(_use_edge) {
		throw oxDNAException("CUDA use_edge is not implemented yet for interaction_type = NA");
	}

	_dna_helper.set_launch_cfg(_launch_cfg);
	_rna_helper.set_launch_cfg(_launch_cfg);
	_dna_helper.compute_forces(lists, d_poss, d_orientations, d_forces, d_torques, d_bonds, d_box);
	_rna_helper.compute_forces(lists, d_poss, d_orientations, d_forces, d_torques, d_bonds, d_box);
	drh_forces_hybrid_hb<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(d_poss, d_orientations, d_forces, d_torques,
		lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_is_dna, d_box);
	CUT_CHECK_ERROR("drh_forces_hybrid_hb error");
	drh_forces_hybrid_excl<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(d_poss, d_orientations, d_forces, d_torques,
		lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_is_dna, d_box);
	CUT_CHECK_ERROR("drh_forces_hybrid_excl error");
	drh_forces_hybrid_cross<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(d_poss, d_orientations, d_forces, d_torques,
		lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_is_dna, d_box);
	CUT_CHECK_ERROR("drh_forces_hybrid_cross error");
	drh_forces_hybrid_coax<<<_launch_cfg.blocks, _launch_cfg.threads_per_block>>>(d_poss, d_orientations, d_forces, d_torques,
		lists->d_matrix_neighs, lists->d_number_neighs, d_bonds, _d_is_dna, d_box);
	CUT_CHECK_ERROR("drh_forces_hybrid_coax error");
}
