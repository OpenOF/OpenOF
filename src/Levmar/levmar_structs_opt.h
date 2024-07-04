/*
 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *

 */


#ifndef LEVMAR_STRUCTS_OPT_H_
#define LEVMAR_STRUCTS_OPT_H_

#include "config.h"
enum robustTyp {OF_NO,OF_HUBER,OF_TRUNCATED_HUBER,OF_L1,OF_TRUNCATED_L1,OF_TRUNCATED_L2,OF_PX_BUNDLE,OF_PX_BUNDLE_L1,OF_PX_BUNDLE_HUBER,OF_PX_BUNDLE_TRUNCATED_HUBER,OF_SORT_WEIGHT};
struct OptStruct_t{
	int nrPara;
	int globalStart;
};



struct MeasurementCombinations_t{
	int type;
	int robust;
	oof_float robust_para_a;
	oof_float robust_para_b;
	oof_float robust_para_c;

	std::vector<int> v;
	std::vector< oof_float > inv_cov;
	std::vector< oof_float > res;
	MeasurementCombinations_t(){
		type=0;
		robust=0;
		robust_para_a=0;
		robust_para_b=0;
		robust_para_c=0;

	}
	MeasurementCombinations_t(int _type,int _v1,int _v2,int _v3,int _v4,int _v5,int _v6,int _v7,int _v8,int _v9,int _v10){
		type=_type;
		v.push_back(_v1);
		v.push_back(_v2);
		v.push_back(_v3);
		v.push_back(_v4);
		v.push_back(_v5);
		v.push_back(_v6);
		v.push_back(_v7);
		v.push_back(_v8);
		v.push_back(_v9);
		v.push_back(_v10);

		robust=0;
		robust_para_a=0;
		robust_para_b=0;
		robust_para_c=0;

	}
};

#endif /* LEVMAR_STRUCTS_OPT_H_ */
