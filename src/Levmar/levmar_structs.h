/*
 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *

 */




#ifndef LEVMAR_STRUTCS_H_
#define LEVMAR_STRUTCS_H_

#include "config.h"
#include <stdio.h>

typedef void(*external_device_set_value_func_t) (void * v1,int ind,oof_float &value,bool set);
typedef void(*external_device_func_t) (void ** v,oof_float * func_result);
typedef void(*external_device_func_jac_t) (void ** v,int rowStart,oof_float * jac_value, int* jac_row_ind,int* jac_col_ind);



struct ValueStruct_t{

	void *v1;
	int ind;
	oof_float value;
	external_device_set_value_func_t* m_func;
	bool set;
	__host__ __device__
	ValueStruct_t(){
		value=0;
		v1=0;
		ind=0;
		m_func=0;
		set=false;
	}
 /*   __host__
	void operator=(oof_float x){
		value=x;
	}*/
	__host__ __device__
	void operator()()
	{
		(**m_func)(v1,ind,value,set);
	}
};
struct MeasurementStruct_t{

	void ** v;

	int robust;
	oof_float robust_para_a;
	oof_float robust_para_b;
	oof_float robust_para_c;
	oof_float * weight_value;

	oof_float * func_result;
	oof_float * h_func_result;
	oof_float * jac_value;
	int * jac_row_ind;
	int * jac_col_ind;

	int nr_row;
	int nr_para;
	int nr_nonzero_jac;
	int globalStart;
	int startJac;
	external_device_func_t* m_func;
	external_device_func_jac_t* m_func_jac;
	__host__ __device__
	MeasurementStruct_t():v(0),globalStart(0),func_result(0),jac_value(0),jac_row_ind(0),jac_col_ind(0){
		robust=0;
		robust_para_a=0;
		robust_para_b=0;
		robust_para_c=0;
	}


	__host__ __device__
	void operator()()
	{

		(**m_func)(v,func_result);
	}
	__host__ __device__
	void jac()
	{

		(**m_func_jac)(v,globalStart,jac_value,jac_row_ind,jac_col_ind);

	}


};


#endif
