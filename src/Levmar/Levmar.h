/*
 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *

 */


#ifndef LEVMAR_H_
#define LEVMAR_H_

#include "config.h"
#include <thrust/host_vector.h>
#include "levmar_structs.h"
#include <cusp/array1d.h>

#ifdef USE_GPU
  #include <thrust/system/cuda/vector.h>
  #define memory_vector thrust::cuda::vector
  #define cusp_memory cusp::device_memory
#else
  #include <thrust/system/tbb/vector.h>
  #define memory_vector thrust::tbb::vector
  #define cusp_memory cusp::host_memory
#endif

#include "meas_func_h.h"
#include "structs_func_h.h"


class Levmar {
public:

	int k;
	int kmax;
	oof_float v;
	oof_float eps1;
	oof_float eps2;
	oof_float tau;

	int cg_it;
	oof_float cg_thresh;
	oof_float cg_thresh_rel;
	oof_float cg_thresh_abs;

	oof_float mue;

	int verbose;

	std::vector<double> residuals;

	Levmar();

	int registerOptObj(int n);
	int insertOptObj(void *obj,external_device_set_value_func_t *func,int n,int globalStart);
	int registerMeasObj(MeasurementStruct_t *m);
	int insertMeasObj(MeasurementStruct_t *m);
	int insertMeasObj(MeasurementStruct_t *m,std::vector<oof_float> &inv_cov);
	int initMemory();
	int init();
	int run();

	void printMemUsage();


	virtual ~Levmar();

	thrust::host_vector<MeasurementStruct_t> h_MeasVec;
	cusp::array1d<oof_float, cusp_memory> var;
#ifdef useHostMemory
	MeasurementStruct_t * h_ptr;
	MeasurementStruct_t * d_ptr;
	thrust::device_ptr<MeasurementStruct_t> d_MeasVec_begin;
	thrust::device_ptr<MeasurementStruct_t> d_MeasVec_end;
#else
	memory_vector<MeasurementStruct_t> d_MeasVec;
#endif

private:
	int row;
	int col;

	int nrParaOpt;
	int nrParaJac;
	int nrParaFunc;

	oof_float max_g;

	memory_vector<ValueStruct_t> d_ValueVec;
	memory_vector<oof_float> d_funcVec;
	thrust::host_vector<oof_float> h_funcVec;
	memory_vector<oof_float> d_weightVec;
	memory_vector<int> d_jac_row_ind;
	memory_vector<int> d_jac_col_ind;
	memory_vector<oof_float> d_jac_value;

	thrust::host_vector<oof_float> h_inv_covVec;

};

#endif /* LEVMAR_H_ */
