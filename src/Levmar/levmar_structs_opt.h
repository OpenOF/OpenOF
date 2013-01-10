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
	int v1;
	int v2;
	int v3;
	int v4;
	int v5;
	int v6;
	int v7;
	int v8;
	int v9;
	int v10;

	std::vector< oof_float > inv_cov;
	std::vector< oof_float > res;
	MeasurementCombinations_t(){
		type=0;
		robust=0;
		robust_para_a=0;
		robust_para_b=0;
		robust_para_c=0;
		v1=0;
		v2=0;
		v3=0;
		v4=0;
		v5=0;
		v6=0;
		v7=0;
		v8=0;
		v9=0;
		v10=0;

	}
	MeasurementCombinations_t(int _type,int _v1,int _v2,int _v3,int _v4,int _v5,int _v6,int _v7,int _v8,int _v9,int _v10){
		type=_type;
		v1=_v1;
		v2=_v2;
		v3=_v3;
		v4=_v4;
		v5=_v5;
		v6=_v6;
		v7=_v7;
		v8=_v8;
		v9=_v9;
		v10=_v10;
		robust=0;
		robust_para_a=0;
		robust_para_b=0;
		robust_para_c=0;

	}
};

#endif /* LEVMAR_STRUCTS_OPT_H_ */
