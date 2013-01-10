/*
 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *
 *   OpenOF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.

 *   OpenOF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.

 *   You should have received a copy of the GNU General Public License
 *   along with OpenOF.  If not, see <http://www.gnu.org/licenses/>.
 */




#ifndef LEVMAR_STRUTCS_H_
#define LEVMAR_STRUTCS_H_

#include "config.h"
typedef void(*external_device_set_value_func_t) (void * v1,int ind,oof_float &value,bool set);
typedef void(*external_device_func_t) (void * v1,void * v2,void * v3,void * v4,void * v5,void * v6,void * v7,void * v8,void * v9,void * v10,oof_float * func_result);
typedef void(*external_device_func_jac_t) (void * v1,void * v2,void * v3,void * v4,void * v5,void * v6,void * v7,void * v8,void * v9,void * v10,int rowStart,oof_float * jac_value, int* jac_row_ind,int* jac_col_ind);



struct ValueStruct_t{

	void *v1;
	int ind;
	oof_float value;
	external_device_set_value_func_t* m_func;
	bool set;
	ValueStruct_t(){
		value=0;
		v1=0;
		ind=0;
		m_func=0;
		set=false;
	}

	void operator=(oof_float x){
		value=x;
	}
	void operator()()
	{
		(**m_func)(v1,ind,value,set);
	}
};
struct MeasurementStruct_t{

	void * v1;
	void * v2;
	void * v3;
	void * v4;
	void * v5;
	void * v6;
	void * v7;
	void * v8;
	void * v9;
	void * v10;


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
	MeasurementStruct_t():v1(0),v2(0),v3(0),v4(0),v5(0),v6(0),v7(0),v8(0),v9(0),v10(0),globalStart(0),func_result(0),jac_value(0),jac_row_ind(0),jac_col_ind(0){
		robust=0;
		robust_para_a=0;
		robust_para_b=0;
		robust_para_c=0;
	}


	__host__ __device__
	void operator()()
	{
		(**m_func)(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,func_result);
	}
	__host__ __device__
	void jac()
	{
		(**m_func_jac)(v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,globalStart,jac_value,jac_row_ind,jac_col_ind);

	}

};


#endif
