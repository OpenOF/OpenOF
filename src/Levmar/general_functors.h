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



#ifndef GENERAL_FUNCTORS_H_
#define GENERAL_FUNCTORS_H_
#include "structs_in.h"

template <typename T>
struct abs_functor
{
    __host__ __device__
        T operator()(const T& x) const {
    		if (x<0.0)
    			return -x ;
    		else
    			return x ;
        }
};

template <typename T>
struct square
{
    __host__ __device__
        T operator()(const T& x) const {
            return x*x ;
        }
};
template <typename T>
struct normal_distribution_functor
{
  const T mue;
  const T sigma;

  normal_distribution_functor(T _mue,T _sigma) : mue(_mue), sigma(_sigma) {}

  __host__ __device__
  T operator()(const T& x) const
  {
    return 1.0/(sigma*sqrt(2*3.14159265359))*exp(-0.5*((x-mue)/sigma)*((x-mue)/sigma));
  }
};
template <typename T>
struct is_equal_to
{

  const T value;

  is_equal_to(T _value) : value(_value) {}

  __host__ __device__
  bool operator()(const T& x) const
  {
    return x==value;
  }
};
struct wrapper_functor_with_jac
{
	__host__ __device__
	void operator()(MeasurementStruct_t& t)
	{
		t();
		t.jac();
	}
};

struct robust_functor
{
	__host__ __device__
	void operator()(MeasurementStruct_t& t)
	{

		for (int i=0;i<t.nr_row;i++){
			if (t.robust==OF_NO){
				t.weight_value[i]=1.0;
			}
			else if (t.robust==OF_HUBER){
				oof_float e=fabs(t.func_result[i]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else
					t.weight_value[i]=sqrt(2.0*t.robust_para_a*e-t.robust_para_a*t.robust_para_a)/e;
			}
			else if (t.robust==OF_HUBER){
				oof_float e=fabs(t.func_result[i]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else
					t.weight_value[i]=sqrt(2.0*t.robust_para_a*e-t.robust_para_a*t.robust_para_a)/e;
			}
			else if (t.robust==OF_TRUNCATED_HUBER){
				oof_float e=fabs(t.func_result[i]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else if (e<t.robust_para_b)
					t.weight_value[i]=sqrt(2.0*t.robust_para_a*e-t.robust_para_a*t.robust_para_a)/e;
				else
					t.weight_value[i]=sqrt(t.robust_para_a*t.robust_para_a+2.0*t.robust_para_a*(t.robust_para_b-t.robust_para_a))/e;
			}
			else if (t.robust==OF_TRUNCATED_L1){
				oof_float e=fabs(t.func_result[i]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0/sqrt(e);
				else
					t.weight_value[i]=sqrt(t.robust_para_a)/e;
			}
			else if (t.robust==OF_TRUNCATED_L2){
				oof_float e=fabs(t.func_result[i]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else
					t.weight_value[i]=t.robust_para_a/e;
			}
			else if (t.robust==OF_PX_BUNDLE){
				oof_float e=sqrt(t.func_result[0]*t.func_result[0]+t.func_result[1]*t.func_result[1]);

				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else
					t.weight_value[i]=t.robust_para_a/e;
			}
			else if (t.robust==OF_PX_BUNDLE_L1){
				oof_float e=sqrt(t.func_result[0]*t.func_result[0]+t.func_result[1]*t.func_result[1]);
				if (e<0.000001)
					t.weight_value[i]=1.0;
				else if (e<t.robust_para_a)
					t.weight_value[i]=1.0/sqrt(e);
				else
					t.weight_value[i]=sqrt(t.robust_para_a)/e;
			}
			else if (t.robust==OF_PX_BUNDLE_HUBER){
				oof_float e=sqrt(t.func_result[0]*t.func_result[0]+t.func_result[1]*t.func_result[1]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else
					t.weight_value[i]=sqrt(2.0*t.robust_para_a*e-t.robust_para_a*t.robust_para_a)/e;
			}
			else if (t.robust==OF_PX_BUNDLE_TRUNCATED_HUBER){
				oof_float e=sqrt(t.func_result[0]*t.func_result[0]+t.func_result[1]*t.func_result[1]);
				if (e<t.robust_para_a)
					t.weight_value[i]=1.0;
				else if (e<t.robust_para_b)
					t.weight_value[i]=sqrt(2.0*t.robust_para_a*e-t.robust_para_a*t.robust_para_a)/e;
				else
					t.weight_value[i]=sqrt(t.robust_para_a*t.robust_para_a+2.0*t.robust_para_a*(t.robust_para_b-t.robust_para_a))/e;
			}
			else if (t.robust==OF_SORT_WEIGHT){
				t.weight_value[i]=-t.robust_para_a;
			}
			else
				t.weight_value[i]=1.0;
		}

	}
};


struct wrapper_functor_set_get
{
	bool set;
	wrapper_functor_set_get(bool _set){
		set=_set;
	}


	__host__ __device__
	void operator()(ValueStruct_t &t)
	{
		t.set=set;
		t();
	}
};

struct wrapper_functor
{

	template <typename Entity>
	__host__ __device__
	void operator()(Entity & t)
	{
		t();
	}
};

struct get_value_functor{
	   __host__ __device__
	   oof_float operator()(const ValueStruct_t& x) const {
	        return x.value;
	   }
};

typedef thrust::tuple<oof_float&, ValueStruct_t&> SetValueTuple;

struct set_value_functor{
    __host__ __device__
        void operator()(SetValueTuple t)
        {
            thrust::get<1>(t).value = thrust::get<0>(t);
        }
};


template <typename T>
struct squared_plus_functor
{
    __host__ __device__
        T operator()(const T& x,const T& y) const {

    		return x + y*y;
        }
};
template <typename T>
struct mul_functor
{
    __host__ __device__
        T operator()(const T& x,const T& y) const {

    		return x * y;
        }
};
template <typename T>
struct inv_functor
{
    __host__ __device__
        T operator()(const T& x) const {

    		if (x<0.0000000001)
    			return 1.0/0.0000000001;
    		else
    			return 1.0/x;
        }
};
template <typename T>
struct square_functor
{
    __host__ __device__
        T operator()(const T& x) const {
            return x*x ;
        }
};


struct square_weight_functor
{
    __host__ __device__
        oof_float operator()(const thrust::tuple<oof_float,oof_float>& x) const {
    	if (thrust::get<1>(x)<0.0000000001)
    		return 0.0f;
    	else
            return thrust::get<0>(x)*thrust::get<0>(x)/thrust::get<1>(x) ;
        }
};



#endif /* GENERAL_FUNCTORS_H_ */
