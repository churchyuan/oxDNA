/*
 * CUDADRHInteraction.h
 *
 *  Created on: 26/apr/2026
 */

#ifndef CUDADRHINTERACTION_H_
#define CUDADRHINTERACTION_H_

#include "CUDABaseInteraction.h"
#include "CUDADNAInteraction.h"
#include "CUDARNAInteraction.h"
#include "../../Interactions/DRHInteraction.h"

/**
 * @brief CUDA scaffold for the DNA/RNA hybrid interaction implemented by DRHInteraction.
 */
class CUDADRHInteraction: public CUDABaseInteraction, public DRHInteraction {
public:
	CUDADRHInteraction();
	virtual ~CUDADRHInteraction();

	uint8_t *_d_is_dna = nullptr;
	int *_d_is_strand_end = nullptr;
	CUDADNAInteraction _dna_helper;
	CUDARNAInteraction _rna_helper;

	void get_settings(input_file &inp) override;
	void get_cuda_settings(input_file &inp) override;
	void cuda_init(int N) override;
	c_number get_cuda_rcut() override {
		return this->get_rcut();
	}

	void compute_forces(CUDABaseList *lists, c_number4 *d_poss, GPU_quat *d_orientations, c_number4 *d_forces, c_number4 *d_torques, LR_bonds *d_bonds, CUDABox *d_box) override;

protected:
	void _on_T_update() override;
	void _init_acid_types();
	void _init_strand_ends(LR_bonds *d_bonds);
};

#endif /* CUDADRHINTERACTION_H_ */
