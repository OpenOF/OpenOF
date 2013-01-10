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




#include <cusp/coo_matrix.h>
#include <cusp/hyb_matrix.h>
#include <cusp/multiply.h>
#include <cusp/transpose.h>
#include <cusp/io/matrix_market.h>

#include <cusp/precond/ainv.h>
#include <cusp/precond/diagonal.h>
#include <cusp/blas.h>
#include <cusp/print.h>
#include <cusp/krylov/cg.h>

#include <thrust/system_error.h>
#include <thrust/transform.h>
#include <thrust/generate.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/for_each.h>
#include <thrust/tuple.h>
#include <thrust/iterator/zip_iterator.h>


#include "Levmar.h"
#include "meas_func_src.h"
#include "structs_func_src.h"
#include "general_functors.h"

#define checkCudaErrors(err)  __checkCudaErrors (err, __FILE__, __LINE__)

inline void __checkCudaErrors(cudaError err, const char *file, const int line )
{
    if(cudaSuccess != err)
    {
        fprintf(stderr, "%s(%i) : CUDA Runtime API error %d: %s.\n",file, line, (int)err, cudaGetErrorString( err ) );
        exit(-1);
    }
}

Levmar::Levmar() {
	kmax=50;
	v=2;
	eps1=1e-10;
	eps2=1e-10;
	tau=0.001;

	cg_it=150;
	cg_thresh_rel=1e-12;
	cg_thresh_abs=1e-12;
	

	verbose=0;
	nrParaOpt=0;
	nrParaJac=0;
	nrParaFunc=0;

}

void createMatrix(memory_vector<int> &d_jac_row_ind,
		memory_vector<int> &d_jac_col_ind,
		memory_vector<oof_float> &d_jac_value,cusp::coo_matrix<int, oof_float, cusp_memory> &jac){

	thrust::copy(d_jac_row_ind.begin(),d_jac_row_ind.end(),jac.row_indices.begin());
	thrust::copy(d_jac_col_ind.begin(),d_jac_col_ind.end(),jac.column_indices.begin());
	thrust::copy(d_jac_value.begin(),d_jac_value.end(),jac.values.begin());

}

void getDiagonalPrec(cusp::coo_matrix<int, oof_float, cusp_memory> &jact,cusp::array1d<oof_float, cusp_memory> &invCov,cusp::coo_matrix<int, oof_float, cusp_memory> &M){

	M.resize(jact.num_rows,jact.num_rows,jact.num_rows);

	jact.sort_by_row();
	thrust::sequence(M.row_indices.begin(), M.row_indices.end());
	thrust::sequence(M.column_indices.begin(), M.column_indices.end());
	memory_vector<int> outkey(jact.num_rows);

	thrust::equal_to<int> binary_pred;
	thrust::plus<oof_float> binary_op;


	memory_vector<oof_float> v(jact.values.size());
	square_functor<oof_float> sq_op;
	mul_functor<oof_float> mul_op;

	thrust::transform(jact.values.begin(),jact.values.end(),v.begin(),sq_op);

	thrust::transform(
	thrust::make_permutation_iterator(invCov.begin(),jact.column_indices.begin()),
	thrust::make_permutation_iterator(invCov.begin(),jact.column_indices.end()),
	v.begin(),v.begin(),mul_op);

	inv_functor<oof_float> unary_op;
	thrust::reduce_by_key(jact.row_indices.begin(),jact.row_indices.end(),v.begin(),outkey.begin(),M.values.begin(),binary_pred,binary_op);

	thrust::transform(M.values.begin(),M.values.end(),M.values.begin(),unary_op);
}
oof_float getDiagonalMax(cusp::coo_matrix<int, oof_float, cusp_memory> &jact){

	jact.sort_by_row();
	memory_vector<int> outkey(jact.num_rows);
	memory_vector<oof_float> diag(jact.num_rows);

	thrust::equal_to<int> binary_pred;
	thrust::plus<oof_float> binary_op;


	memory_vector<oof_float> v(jact.values.size());
	square_functor<oof_float> sq_op;



	thrust::transform(jact.values.begin(),jact.values.end(),v.begin(),sq_op);
	//squared_plus_functor<oof_float> binary_op;

	thrust::reduce_by_key(jact.row_indices.begin(),jact.row_indices.end(),v.begin(),outkey.begin(),diag.begin(),binary_pred,binary_op);

	oof_float max=thrust::reduce(diag.begin(),diag.end(),0.0,thrust::maximum<oof_float>());
	return max;
}
memory_vector<oof_float> getGaussWeightVector(int row,oof_float mue,oof_float sigma){
	memory_vector<oof_float> x(row);
	memory_vector<oof_float> gaussweight(row);
	thrust::sequence(x.begin(),x.end());
	normal_distribution_functor<oof_float> func(mue,sigma);
	thrust::transform(x.begin(),x.end(),gaussweight.begin(),func);
	thrust::plus<oof_float> plus_func;
	oof_float init=0.0;
	oof_float sum=thrust::reduce(gaussweight.begin(),gaussweight.end(),init,plus_func);
	memory_vector<oof_float> temp(row);
	thrust::fill(temp.begin(),temp.end(),sum);
	thrust::divides<oof_float> divid_func;
	thrust::transform(gaussweight.begin(),gaussweight.end(),temp.begin(),gaussweight.begin(),divid_func);
	return gaussweight;
}



//it returns the global start
int Levmar::registerOptObj(int n){
	int returnInd=nrParaOpt;
	nrParaOpt+=n;
	return returnInd;
}

int Levmar::insertOptObj(void *obj,external_device_set_value_func_t *func,int n,int globalStart){
	for (int i=0;i<n;i++){
		ValueStruct_t v;
		v.m_func=func;
		v.ind=i;
		v.v1=obj;
		d_ValueVec[globalStart+i]=v;
	}
	return 1;
}
int Levmar::registerMeasObj(MeasurementStruct_t *m){
	m->globalStart=nrParaFunc;
	m->startJac=nrParaJac;
	nrParaFunc+=m->nr_row;
	nrParaJac+=m->nr_nonzero_jac;

	return 1;
}
int Levmar::insertMeasObj(MeasurementStruct_t *m){

	m->weight_value=thrust::raw_pointer_cast(&d_weightVec[m->globalStart]);
	m->func_result=thrust::raw_pointer_cast(&d_funcVec[m->globalStart]);
	m->jac_value=thrust::raw_pointer_cast(&d_jac_value[m->startJac]);
	m->jac_row_ind=thrust::raw_pointer_cast(&d_jac_row_ind[m->startJac]);
	m->jac_col_ind=thrust::raw_pointer_cast(&d_jac_col_ind[m->startJac]);
	h_MeasVec.push_back(*m);
	for (int i=0;i<m->nr_row;i++){
		h_inv_covVec[m->globalStart+i]=1.0;
	}

	return 1;
}
int Levmar::insertMeasObj(MeasurementStruct_t *m,std::vector<oof_float> &inv_cov){

	m->weight_value=thrust::raw_pointer_cast(&d_weightVec[m->globalStart]);
	m->func_result=thrust::raw_pointer_cast(&d_funcVec[m->globalStart]);
	m->h_func_result=thrust::raw_pointer_cast(&h_funcVec[m->globalStart]);
	m->jac_value=thrust::raw_pointer_cast(&d_jac_value[m->startJac]);
	m->jac_row_ind=thrust::raw_pointer_cast(&d_jac_row_ind[m->startJac]);
	m->jac_col_ind=thrust::raw_pointer_cast(&d_jac_col_ind[m->startJac]);
	h_MeasVec.push_back(*m);
	for (int i=0;i<m->nr_row;i++){
		if (i<inv_cov.size())
			h_inv_covVec[m->globalStart+i]=inv_cov[i];
		else
			h_inv_covVec[m->globalStart+i]=1.0;
	}

	return 1;
}
int Levmar::initMemory(){

	if (verbose>0)	std::cout<<"nrParameters:"<<nrParaOpt<<"\n";
	if (verbose>0)	std::cout<<"nrParaJac"<<nrParaJac<<"\n";
	if (verbose>0)	std::cout<<"nrParaFunc"<<nrParaFunc<<"\n";
	d_weightVec.resize(nrParaFunc);
	d_funcVec.resize(nrParaFunc);
	h_funcVec.resize(nrParaFunc);
	h_inv_covVec.resize(nrParaFunc);
	d_jac_row_ind.resize(nrParaJac);
	d_jac_col_ind.resize(nrParaJac);
	d_jac_value.resize(nrParaJac);
	//h_ValueVec.resize(nrParaOpt);
	d_ValueVec.resize(nrParaOpt);
	return 1;
}
int Levmar::init()
{
#ifdef useHostMemory
	int N=h_MeasVec.size();
	checkCudaErrors(cudaMallocHost((void**)&h_ptr, N*sizeof(MeasurementStruct_t)));
	checkCudaErrors(cudaHostGetDevicePointer((void **)&d_ptr, (void *)h_ptr,0));
	d_MeasVec_begin=thrust::device_pointer_cast(d_ptr);
	d_MeasVec_end=	d_MeasVec_begin+N;
	thrust::copy(h_MeasVec.begin(),h_MeasVec.end(),d_MeasVec_begin);
#else
	d_MeasVec.resize(h_MeasVec.size());
	thrust::copy(h_MeasVec.begin(),h_MeasVec.end(),d_MeasVec.begin());
#endif
	return 1;
}


void Levmar::printMemUsage(){
	size_t avail;
	size_t total;
	cudaMemGetInfo( &avail, &total );
	size_t used = total - avail;
	std::cout << "Device memory used: " << used << std::endl;
}



template <class LinearOperator1,
		  class LinearOperator2,
          class Vector,
          class Monitor,
          class Preconditioner>
void cgls(LinearOperator1& At,
		LinearOperator2& A,
        Vector& x,
        Vector& b,
        Vector& invCov,
        Monitor& monitor,
        Preconditioner& M, oof_float mu)
{
    CUSP_PROFILE_SCOPED();

    typedef typename LinearOperator1::value_type   ValueType;
    typedef typename LinearOperator1::memory_space MemorySpace;
    typedef typename LinearOperator2::value_type   ValueType2;
    typedef typename LinearOperator2::memory_space MemorySpace2;

   // assert(A.num_rows == A.num_cols);        // sanity check

    const size_t N = A.num_cols;
    const size_t N2 = A.num_rows;

    // allocate workspace
    cusp::array1d<ValueType,MemorySpace> y(N);
    cusp::array1d<ValueType,MemorySpace> y1(N2);

    cusp::array1d<ValueType,MemorySpace> z(N);
    cusp::array1d<ValueType,MemorySpace> r(N);
    cusp::array1d<ValueType,MemorySpace> p(N);

    // y <- Ax
    cusp::multiply(A, x, y1);
    blas::xmy(invCov,y1,y1);
    cusp::multiply(At, y1, y);
    blas::axpy(x,y,mu);
    // r <- b - A*x
    blas::axpby(b, y, r, ValueType(1), ValueType(-1));

    // z <- M*r
    cusp::multiply(M, r, z);

    // p <- z
    blas::copy(z, p);

    // rz = <r^H, z>
    ValueType rz = blas::dotc(r, z);

    while (!monitor.finished(r))
    {
        // y <- Ap
        //cusp::multiply(A, p, y);

        cusp::multiply(A, p, y1);
        blas::xmy(invCov,y1,y1);
        cusp::multiply(At, y1, y);
        blas::axpy(p,y,mu);

        // alpha <- <r,z>/<y,p>
        ValueType alpha =  rz / blas::dotc(y, p);
        //std::cout<<"alpha="<<alpha<<"\n";
        // x <- x + alpha * p
        blas::axpy(p, x, alpha);

        // r <- r - alpha * y
        blas::axpy(y, r, -alpha);

        // z <- M*r
        cusp::multiply(M, r, z);

        ValueType rz_old = rz;

        // rz = <r^H, z>
        rz = blas::dotc(r, z);
        //std::cout<<"rz="<<rz<<"\n";
        //std::cout<<"rz_old="<<rz_old<<"\n";

        // beta <- <r_{i+1},r_{i+1}>/<r,r>
        ValueType beta = rz / rz_old;
        if (beta==0.0)
        	break;
        //std::cout<<"beta="<<beta<<"\n";
        // p <- r + beta*p
        blas::axpby(z, p, p, ValueType(1), beta);

        ++monitor;
    }
}

int Levmar::run(){


	timespec ts;
	timespec te;
	clock_gettime(CLOCK_REALTIME, &ts);

	bool found=false;
	double Fx=0.0;
	double Fxnew=0.0;

	k=0;
	residuals.clear();

	//get row and col of jacobian matrix
	row=d_funcVec.size();
	col=d_ValueVec.size();

	//variables for the current state, update state, and function value
	cusp::array1d<oof_float, cusp_memory> x(col);
	cusp::array1d<oof_float, cusp_memory> xnew(col);
	cusp::array1d<oof_float, cusp_memory> f(d_funcVec.size());
	cusp::array1d<oof_float, cusp_memory> w_robust(row);

	cusp::array1d<oof_float, cusp_memory> invCov(row);
	cusp::array1d<oof_float, cusp_memory> invCovRobust(row);
	thrust::copy(h_inv_covVec.begin(),h_inv_covVec.end(),invCov.begin());

	//variables for matrix jac
	cusp::coo_matrix<int, oof_float, cusp_memory> jac(row,col,d_jac_row_ind.size());
	//variables for matrix W for weighting
	//cusp::coo_matrix<int, oof_float, cusp_memory> W(row,row,row);

	//further matrixes needed
	cusp::coo_matrix<int, oof_float, cusp_memory> jact;
	cusp::coo_matrix<int, oof_float, cusp_memory> jactW;
	cusp::coo_matrix<int, oof_float, cusp_memory> jtj;
    cusp::coo_matrix<int, oof_float, cusp_memory> jtj_mueI;
    cusp::coo_matrix<int, oof_float, cusp_memory> mueI(col, col,col);



	cusp::array1d<oof_float, cusp_memory> g(col);
	cusp::array1d<oof_float, cusp_memory> h(col);

    //calculate cost;
	#ifdef useHostMemory
		thrust::for_each(d_MeasVec_begin,d_MeasVec_end,wrapper_functor_with_jac());
		thrust::for_each(d_MeasVec_begin,d_MeasVec_end,robust_functor());
	#else
		thrust::for_each(d_MeasVec.begin(),d_MeasVec.end(),wrapper_functor_with_jac());
		thrust::for_each(d_MeasVec.begin(),d_MeasVec.end(),robust_functor());
	#endif
    thrust::copy(d_weightVec.begin(),d_weightVec.end(),w_robust.begin());
    //copy values from ValueStruct vector to x
	thrust::for_each(d_ValueVec.begin(), d_ValueVec.end(),
			wrapper_functor_set_get(false));


    thrust::transform(d_ValueVec.begin(),d_ValueVec.end(),x.begin(),get_value_functor());

	//copy values from funcVec to f
	thrust::copy(d_funcVec.begin(),d_funcVec.end(),f.begin());

	//create jacobi matrix from vectors
	createMatrix(d_jac_row_ind,d_jac_col_ind,d_jac_value,jac);

	//update matrix struct to have the real size
	jac.resize(row,col,jac.row_indices.size());
	double r=cusp::blas::dot(f,f);

	residuals.push_back(r);

	cusp::transpose(jac, jact);
	//cusp::multiply(jact,W,jactW);

	blas::xmy(invCov,w_robust,invCovRobust);
	blas::xmy(invCovRobust,w_robust,invCovRobust);
	blas::xmy(invCovRobust,f,f);
	//compute g
    cusp::multiply(jact, f, g);

    //-g
    cusp::blas::scal(g,-1.0);

    //compute jact*W*jac
    //cusp::multiply(jact,jac, jtj);

    //get initial mue
    {
        //cusp::array1d<oof_float,cusp_memory> d_diagA;
		//cusp::detail::extract_diagonal(jtj, d_diagA);
		//mue= tau*thrust::reduce(d_diagA.begin(), d_diagA.end(), (oof_float) 0, thrust::maximum<oof_float>());
    	//TODO:: weighted diagonal max nehmen
    	mue= tau*getDiagonalMax(jact);
    	mue=1.0;
	}
    if (verbose>2)	std::cout<<"initial mue:"<<mue<<"\n";
    //check if all values in g (infinity norm) are less than eps1
    max_g = thrust::transform_reduce(g.begin(), g.end(), abs_functor<oof_float>(), 0.0, thrust::maximum<oof_float>());

    //solution found
    if (max_g<eps1)
    	found =true;

    //compute residual same as dot product

    Fx=thrust::transform_reduce(
    		thrust::make_zip_iterator(thrust::make_tuple(f.begin(), invCovRobust.begin())),
    		thrust::make_zip_iterator(thrust::make_tuple(f.end(), invCovRobust.end())),
    		square_weight_functor(), 0.0, thrust::plus<oof_float>()) ;




    //Fx=  thrust::transform_reduce(f.begin(), f.end(), square<oof_float>(), 0.0, thrust::plus<oof_float>()) ;


    if (verbose>0) 	std::cout<<"iteration:"<<k<<" residual:"<<Fx<<"\n";
	clock_gettime(CLOCK_REALTIME, &te);
	if (verbose>3)	std::cout<<"Start Running Loop:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";

    while (!found && (k<kmax)){
    	k++;
    	clock_gettime(CLOCK_REALTIME, &te);
    	if (verbose>3) std::cout<<"Loop:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";

				//create diagonal matrix with mue at the diagonal
				//thrust::sequence(mueI.row_indices.begin(), mueI.row_indices.end());
				//thrust::sequence(mueI.column_indices.begin(), mueI.column_indices.end());
				//thrust::fill(mueI.values.begin(),mueI.values.end(),mue);

    			//(jactWjac+mueI)
    			//cusp::add(jtj,mueI, jtj_mueI);

        if (verbose>3) 	std::cout<<"solve: jact * W * jac    * h  = - g\n";
        thrust::fill(h.begin(),h.end(),0.0);
    //	thrust::copy(g.begin(),g.end(),h.begin());
        cusp::default_monitor<oof_float> monitor(h,cg_it,cg_thresh_rel,cg_thresh_abs);
       // cusp::precond::smoothed_aggregation<int,oof_float,	cusp::device_memory> M(jtj_mueI);
        //cusp::precond::bridson_ainv<oof_float,	cusp::device_memory> M(jtj_mueI);
        //cusp::precond::diagonal<oof_float, cusp_memory> M(jtj_mueI);

        cusp::coo_matrix<int, oof_float, cusp_memory> M;
        getDiagonalPrec(jact,invCovRobust,M);
        // SOLVE jact * W * jac    * h  = - g
        clock_gettime(CLOCK_REALTIME, &te);
    	if (verbose>3) std::cout<<"Solve Normal EQ start:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";
  	    //cusp::krylov::cg(jtj_mueI, h, g,monitor,M);


  	    cgls(jact,jac, h, g,invCovRobust,monitor,M,mue);
  	    clock_gettime(CLOCK_REALTIME, &te);
  	    if (verbose>3) std::cout<<"Solve Normal EQ end:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";
					//cusp::krylov::cg(jtj_mueI, h, g);
					//if (!monitor.converged()){
					//  	    	mue=mue*10.0;
					//  	    	continue;
					//  	    }

    	oof_float hnorm=cusp::blas::nrm2(h);
    	oof_float xnorm=cusp::blas::nrm2(x);

    	if (verbose>4) std::cout<<"hnorm:"<<hnorm<<"\n";
    	//check if change is to small, is so we are finished
    	if (hnorm<eps2*(xnorm+eps2)){
    		found=true;
    		if (verbose>4) std::cout<<"solution found "<<hnorm<<" < " <<eps2*(xnorm+eps2)<<"\n";
    	}else{
    		//xnew = x+h
    		cusp::blas::axpby(x,h,xnew,1.0,1.0);

    		//compute gain
    		cusp::array1d<oof_float, cusp_memory> temp(d_ValueVec.size());
    		//mue*h-g, g is alreade negated
    		cusp::blas::axpby(h,g,temp,mue,1.0);
    		oof_float gainLower=cusp::blas::dot(h,temp);

    		//schreibe werte aus x nach d_value
    		thrust::copy(xnew.begin(),xnew.end(),d_ValueVec.begin());

    		if (verbose>4) std::cout<<"copy xnew to original structs\n";
    		//schreibe werte in original structs
    		thrust::for_each(d_ValueVec.begin(),d_ValueVec.end(),wrapper_functor_set_get(true));
    		if (verbose>4) std::cout<<"evaluate function with xnew\n";
    		//evaluate function with xnew
    		clock_gettime(CLOCK_REALTIME, &te);
    		if (verbose>3)
    			std::cout<<"EVAL f and JAC start:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";
			#ifdef useHostMemory
				thrust::for_each(d_MeasVec_begin,d_MeasVec_end,wrapper_functor_with_jac());
				thrust::for_each(d_MeasVec_begin,d_MeasVec_end,robust_functor());
			#else
				thrust::for_each(d_MeasVec.begin(),d_MeasVec.end(),wrapper_functor_with_jac());
				thrust::for_each(d_MeasVec.begin(),d_MeasVec.end(),robust_functor());
			#endif

        	thrust::copy(d_weightVec.begin(),d_weightVec.end(),w_robust.begin());

        	clock_gettime(CLOCK_REALTIME, &te);
        	if (verbose>3) std::cout<<"EVAL F and JAC end:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";
         	//compute Fxnew from evaluation of F with xnew
    		thrust::copy(d_funcVec.begin(),d_funcVec.end(),f.begin());

    		//new residual
    		blas::xmy(invCov,w_robust,invCovRobust);
    		blas::xmy(invCovRobust,w_robust,invCovRobust);
    		blas::xmy(invCovRobust,f,f);
    		Fxnew=thrust::transform_reduce(
    		    		thrust::make_zip_iterator(thrust::make_tuple(f.begin(), invCovRobust.begin())),
    		    		thrust::make_zip_iterator(thrust::make_tuple(f.end(), invCovRobust.end())),
    		    		square_weight_functor(), 0.0, thrust::plus<oof_float>()) ;

//         	Fxnew=  thrust::transform_reduce(f.begin(), f.end(), square<oof_float>(), 0.0, thrust::plus<oof_float>()) ;
         	if (verbose>2) std::cout<<"proposed residual:"<<Fxnew<<"\n";
         	//gain
         	if (gainLower<0)
         		gainLower*=-1.0;
         	oof_float gain=(Fx-Fxnew)/gainLower;

         	if (verbose>4) std::cout<<"gain: "<<gain<<"\n";
         	if (verbose>4) std::cout<<"gainLower: "<<gainLower<<"\n";

         	//if positive gain update x
         	if (gain>0) {
         		//x=xnew
         		thrust::copy(xnew.begin(),xnew.end(),x.begin());
         		//fill jac matrix
         		createMatrix(d_jac_row_ind,d_jac_col_ind,d_jac_value,jac);
         		//copy residuals to f
         		thrust::copy(d_funcVec.begin(),d_funcVec.end(),f.begin());

        		cusp::transpose(jac, jact);
        		 //printMemUsage();

        	//	cusp::multiply(jact,W,jactW);
            	clock_gettime(CLOCK_REALTIME, &te);
            	if (verbose>3) std::cout<<"JacT * Jac start:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";

        		//cusp::multiply(jact,jac, jtj);


        		clock_gettime(CLOCK_REALTIME, &te);
        		if (verbose>3) std::cout<<"JacT * Jac end:"<<(double)(te.tv_sec - ts.tv_sec)<<"sec + "<<(double)(te.tv_nsec - ts.tv_nsec)/1000000 <<"ms\n";
        		blas::xmy(invCov,w_robust,invCovRobust);
        		blas::xmy(invCovRobust,w_robust,invCovRobust);
        		blas::xmy(invCovRobust,f,f);
        	    cusp::multiply(jact, f, g);

        	    cusp::blas::scal(g,-1.0);

        	    //infinity norm of g
        	    max_g = thrust::transform_reduce(g.begin(), g.end(), abs_functor<oof_float>(),(oof_float) 0.0, thrust::maximum<oof_float>());
        	    if (verbose>4) std::cout<<"max g"<<max_g<<"\n";

        	    if (max_g<eps1)
        	    	found =true;
        	    mue*=0.1;
        	    //mue=mue*max(0.33333333334,1.0-pow(2.0*gain-1.0,3));
        	    mue=std::max((oof_float)1e-10,mue);
        	    if (verbose>4) 	std::cout<<"mue:"<<mue<<"\n";

        	    v=2.0;
        	    Fx=thrust::transform_reduce(
        	    		thrust::make_zip_iterator(thrust::make_tuple(f.begin(), invCovRobust.begin())),
        	    		thrust::make_zip_iterator(thrust::make_tuple(f.end(), invCovRobust.end())),
        	    		square_weight_functor(), 0.0, thrust::plus<oof_float>()) ;
             	//Fx=  thrust::transform_reduce(f.begin(), f.end(), square<oof_float>(), 0.0, thrust::plus<oof_float>()) ;

			    if (verbose>1) 	std::cout<<"iteration:"<<k<<" residual:"<<Fx<<"\n";

        	    residuals.push_back(Fx);
		   // if (Fx<0.002)
			//break;
         	}

        	else{
        		residuals.push_back(Fx);
        	    mue=mue*10.0;
        	   // mue=mue*v;
        	    if (verbose>4) 	std::cout<<"mue:"<<mue<<"\n";
        		v=2*v;
        	}

    	}
    }

    if (verbose>0) 	std::cout<<"iteration:"<<k<<" residual:"<<Fx<<"\n";
	//copy solution

	//write value from x to d_value
	thrust::copy(x.begin(),x.end(),d_ValueVec.begin());

	//write values to original structs
	thrust::for_each(d_ValueVec.begin(),d_ValueVec.end(),wrapper_functor_set_get(true));
#ifdef useHostMemory
	thrust::for_each(d_MeasVec_begin,d_MeasVec_end,wrapper_functor_with_jac());
#else
	thrust::for_each(d_MeasVec.begin(),d_MeasVec.end(),wrapper_functor_with_jac());
#endif

	thrust::copy(d_funcVec.begin(),d_funcVec.end(),h_funcVec.begin());

	return 1;
}

Levmar::~Levmar() {
}


