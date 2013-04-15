'''

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


This file contains all methods to generate source code from your model designed with sympy.
'''
import sympy
from sympy import printing
import numpy
import sys
import numpy
import shutil
import glob
import os

NR_OPT_PARA=10

def countNonZero(jac):
    count=0
    for i in range(jac.shape[0]):
        for j in range(jac.shape[1]):
            if jac[i,j]!=0:
                count+=1
    return count

def genOptCodeStruct(name,variables,addVariables,addCode=''):
    """
    Method to generate the struct for your optimization objects

    :param name: string containing the name of the optimization objects
    :param variables: string containing the variable names of the optimization object, seperated by ",", (will by optimized), type will by oof_float
    :param addvariables: string containing the variable names of the optimization object seperated by "," (will not optimized), type will by oof_float
    :param addCode: String containing additional code that will be placed in the struct.
    
    Returns:    
        c++ code which will be placed in struct.h
            
    Examples:
        genOptCodeStruct('pt2','x,y','','int id;')
    """
    var=None
    if len(variables)>0:
        var=sympy.symbols(variables)
    addvar=None
    if len(addVariables)>0:
        addvar=sympy.symbols(addVariables)
    outstr='struct %s_t:OptStruct_t{\n'%name
    outstr+='\t%s_t(){\n'%name
    if var is not None:
        outstr+='\t\tnrPara=%d;\n'%len(var)
        for x in var:
            outstr+='\t\t%s=0.0;\n'%x
    else:
        outstr+='\t\tnrPara=0;\n'
    if addvar is not None:
        for x in addvar:
            outstr+='\t\t%s=0.0;\n'%x
    outstr+='\t}\n'
    if var is not None:
        for x in var:
            outstr+='\toof_float %s;\n'%x
    if addvar is not None:
        for x in addvar:
            outstr+='\toof_float %s;\n'%x

    outstr+='\n'
    outstr+='\toof_float get(int ind){\n'
    i=0
    if var is not None:
        for x in var:            
            outstr+='\t\t if (ind==%d) return %s;\n'%(i,x)
            i+=1
    if addvar is not None:
        for x in addvar:
            outstr+='\t\t if (ind==%d) return %s;\n'%(i,x)
            i+=1  
    outstr+='\t\t return 0.0;\n'
    outstr+='\t}\n'
    
    outstr+='\n'
    outstr+='\tvoid set(int ind,oof_float _inputValue){\n'
    i=0
    if var is not None:
        for x in var:            
            outstr+='\t\t if (ind==%d) %s=_inputValue;\n'%(i,x)
            i+=1
    if addvar is not None:
        for x in addvar:
            outstr+='\t\t if (ind==%d) %s=_inputValue;\n'%(i,x)
            i+=1  
    outstr+='\t}\n'
    
    outstr+=addCode
    outstr+='};'
    return outstr

def genOptCodeSet(name,variables):
    """
    Method to generate C++ code to set or get a variable by a number

    :param name: string containing the name of the optimization objects
    :param variables: string containing the variable names of the optimization object, seperated by ",", (will by optimized), type will by oof_float
   
    Returns:

	c++ code which will be placed in struct_func_src.h
    
    Examples:
        genOptCodeSet('pt2','x,y')


    """
    var=None
    if len(variables)>0:
        var=sympy.symbols(variables)

    outstr='__host__ __device__\n'
    outstr+='void set_value_%s_func(void * %s,int ind, oof_float &value,bool set){\n'%(name,name)
    if var is not None:    
        outstr+='\t%s_t * v=(%s_t*)%s;\n'%(name,name,name)
        outstr+='\tif (set){\n'
        for i,x in zip(range(len(var)),var):
            outstr+='\t\tif (ind == %d)\n'%i
            outstr+='\t\t\tv->%s=value;\n'%str(x)
        outstr+='\t}else{\n'
        for i,x in zip(range(len(var)),var):
            outstr+='\t\tif (ind == %d)\n'%i
            outstr+='\t\t\tvalue=v->%s;\n'%str(x)
        outstr+='\t}\n'
    outstr+='}\n'
    
    return outstr

def genOptCodeFuncPtr(name):
    """
    Method generates C++ code. The source code registers a function pointer on the GPU.

    :param name: string containing the name of the optimization objects
   
    Returns:

    c++ code which will be placed in struct_func_src.h
    
    Examples:
        genOptCodeFuncPtr('pt2')

    """
    
    outstr='__global__ void set_value_%s_ptr('%name
    outstr+='external_device_set_value_func_t* dev_func_ptr){\n'
    outstr+='\t*dev_func_ptr = set_value_%s_func;\n'%name
    outstr+='}\n';
    outstr+='\n';
    outstr+='void set_value_%s_ptr_h('%name
    outstr+='external_device_set_value_func_t* dev_func_ptr){\n'
    outstr+='\t*dev_func_ptr = set_value_%s_func;\n'%name
    outstr+='}\n';
    return outstr
    
def genOptCodeFuncPtr_h(name):    
    """
    Method generates C++ code. Header of the function generated in genOptCodeFuncPtr

    :param name: string containing the name of the optimization objects
   
    Returns:

    c++ code which will be placed in struct_func_src.h
    
    Examples:
        genOptCodeFuncPtr_h('pt2')

    """
    outstr='__global__ void set_value_%s_ptr(external_device_set_value_func_t* dev_func_ptr);'%name
    outstr+='\n'
    outstr+='void set_value_%s_ptr_h(external_device_set_value_func_t* dev_func_ptr);'%name
    return outstr 

def genOptCodeStructFile(names,variables,addVariables,addCode=None,filename='structs',openof_root=None):
    """
    Method generates C++ code and writes it to structs.h. To generate the code the genOptCodeStruct function is used.

    :param names: List of strings containing the names of the optimization objects
    :param variables: List of strings containing the variable names of the optimization object, seperated by ",", (will by optimized), type will by oof_float
    :param addVariables: List of strings containing the variable names of the optimization object seperated by "," (will not optimized), type will by oof_float
    :param addCode: List of strings containing additional code that will be placed in the structs.h.
    :param filename: The filename of the structs can be adjusted. This can be necessary for example if you want to use to OpenOF libraries in the same project.


    Returns:

    c++ code which will be placed in structs.h
    
    Examples:
    
        names=['pt2','pt3','cam','cam_in']

        variables=['x,y','X,Y,Z','CX,CY,CZ,q1,q2,q3,q4','fx,fy,u0,v0,k1,k2,k3']
    
        addVariables=['','','',['s']]
    
        addCode=['int id;','int id;int r,g,b;','','']
    
        genOptCodeStructFile(names,variables,addVariables,addCode)
    """
    
    
    outstr='#ifndef %s_H_\n'%filename.upper()
    outstr+='#define %s_H_\n'%filename.upper()
    outstr+='#include "levmar_structs_opt.h"\n\n'
    if addCode is None:
        addCode=numpy.zeros((len(names)),dtype=str)
    for n,v,a,c in zip(names,variables,addVariables,addCode):
        outstr+='\n'
        outstr+=genOptCodeStruct(n, v,a,c)
        outstr+='\n'
    outstr+='#endif'
    f=open(openof_root+'/src/Model/'+filename+'.h','w')
    f.write(outstr)
    f.close
    f=open(openof_root+'/src/Model/'+'structs_in.h','w')
    f.write('#include "%s.h"'%filename)
    f.close
    return outstr

def genOptCodeStructFuncFile(names,variables,openof_root):
    """
    Method generates C++ code and writes it to struct_func_src.h. To generate the code, the genOptCodeFuncPtr, genOptCodeSet, genOptCodeFuncPtr_h function is used.

    :param names: List of strings containing the names of the optimization objects
    :param variables: List of strings containing the variable names of the optimization object, seperated by ",", (will by optimized), type will by oof_float


    Returns:

    c++ code which will be placed in struct_func_src.h
    
    Examples:
    
        names=['pt2','pt3','cam','cam_in']

        variables=['x,y','X,Y,Z','CX,CY,CZ,q1,q2,q3,q4','fx,fy,u0,v0,k1,k2,k3']
    
        genOptCodeStructFuncFile(names,variables)
    """
    
    
    outstr='#ifndef STRUCTS_FUNC_SRC_H_\n'
    outstr+='#define STRUCTS_FUNC_SRC_H_\n'
    outstr+='\n'
    
    for n,v in zip(names,variables):
        outstr+='\n'
        outstr+=genOptCodeSet(n,v)
        outstr+='\n'
        
        
        
    for n,v in zip(names,variables):
        outstr+='\n'
        outstr+=genOptCodeFuncPtr(n)
        outstr+='\n'
        
    outstr+='#endif'
    f=open(openof_root+'/src/Model/'+'structs_func_src.h','w')
    f.write(outstr)
    f.close
    outstr_src=outstr
    
    outstr='#ifndef STRUCTS_FUNC_H_H_\n'
    outstr+='#define STRUCTS_FUNC_H_H_\n'
    
    for n,v in zip(names,variables):
        outstr+='\n'
        outstr+=genOptCodeFuncPtr_h(n)
        outstr+='\n'
        
    outstr+='#endif'
    f=open(openof_root+'/src/Model/'+'structs_func_h.h','w')
    f.write(outstr)
    f.close
    
    

    
    return outstr,outstr_src
     

def genOptCodeMeas(name,optNames,jac):
    
    """
    Method generates C++ code, for the optimization measurement struct.

    :param name: Name of the measurement
    :param optNames: Names of the optimization objects, must be identical with the structs
    :param jac: Jacobian Matrix of the measurement
    


    Returns:

    c++ code which will be placed in struct_meas.h
    
    Examples:
  
        genOptCodeMeas('ProjectMetric',['pt2','pt3','cam','cam_in'],jac.shape[0],jac.shape[1]))
        
        with jac being the computed jacobian of the underlying function
    """
    
    
    #if type(optNames)==str:
    #optNames=tuple(optNames)
    nrFunc=jac.shape[0]
    nrPara=jac.shape[1]
    nrNonZeroJac=countNonZero(jac)
    outStr='struct Measurement_%s_t:MeasurementStruct_t{\n'%name
    #outStr+='\tMeasurement_%s_t(external_device_func_t* func,external_device_func_jac_t* func_jac,'%name
    outStr+='\tMeasurement_%s_t(external_device_func_t* func,external_device_func_jac_t* func_jac,void **_v)'%name
#    for i in range(len(optNames)):
#        outStr+='%s_t *_v%d' % (optNames[i],i+1)
#        if i!=len(optNames)-1:
#            outStr+=','
#    outStr+='){\n'
    outStr+='{\n'
    outStr+='\t\tm_func=func;\n\t\tm_func_jac=func_jac;\n'
    #for i in range(len(optNames)):
    #    outStr+='\t\tv%d=(void *)_v%d;\n'%(i+1,i+1)
    
#    outStr+='\t\tvoid **v_temp=(void**)malloc(%d * sizeof(void*));\n'%len(optNames)
#    for i in range(len(optNames)):
#        outStr+='\t\tv_temp[%d]=(void *)_v%d;\n'%(i,i+1);
#    outStr+='#ifdef USE_GPU\n'
#    outStr+='\t\tcudaMalloc((void**)&v, %d*sizeof(void*));\n'%len(optNames)    
#    outStr+='\t\tcudaMemcpy(v,v_temp,%d*sizeof(void*),cudaMemcpyHostToDevice);\n'%len(optNames)
#    outStr+='\t\t delete[] v_temp;\n'
#    outStr+='#else\n'
#    outStr+='\t\tv=v_temp;\n';
 #   outStr+='#endif\n'
    outStr+='\t\tv=_v;\n';
        
    #outStr+='\t\tnr_opt_obj=%d;\n'%len(optNames)
    outStr+='\t\tnr_row=%d;\n'%nrFunc
    outStr+='\t\tnr_para=%d;\n'%nrPara
    outStr+='\t\tnr_nonzero_jac=%d;\n'%nrNonZeroJac
    outStr+='\t}\n';
    outStr+='};\n';
    
    return outStr
def genOptCodeMeasFile(str,openof_root): 
    """
    Method writes structs_meas.h.

    :param str: Lists of strings
    
    """
    
    
        
    outStr='#ifndef STRUCTS_MEAS_H_\n'
    outStr+='#define STRUCTS_MEAS_H_\n'
    for s in str:
        outStr+='\n'
        outStr+=s
        outStr+='\n'
    outStr+='#endif'
    f=open(openof_root+'/src/Model/'+'structs_meas.h','w')
    f.write(outStr)
    f.close
def genOptCodeMeasIterator(measFuncName,names,measFuncInputJac,measFuncInputObj):
    """
    Iterator that runs over all MeasFuncName. It can be used to create a list for genOptCodeMeasFile
    :param measFuncName: List with names of all measurements
    :param names: List with names of all objects
    :param measFuncInputJac: List of Jacobian matrices, corresponding to the measFuncName list
    :param measFuncInputObj: List of Lists. Each inner list correspond to one measurement function and contains the indices of objects that are part of this measurement. The indice belongs the names list.
    """
    
    for i,n,jac in zip(range(len(measFuncName)),measFuncName,measFuncInputJac):
        yield genOptCodeMeas(n,[names[j] for j in measFuncInputObj[i]],jac)

def genOptCodeFuncMeasIterator(measFuncName,names,variables,addVariables,measFuncInputF,measFuncInputJac,measFuncInputObj,measFuncInputObjOpt,addFunc=[],debug=False):
    """
    Iterator that can be used to create a list for genOptCodeFuncMeasFile.
    :param measFuncName: List with names of all measurements
    :param names: List with names of all objects
    :param variables: List of optimization variables belonging to names
    :param addVariables: List of constant variables belonging to names
    :param measFuncInputF: List of functions, corresponding to measFuncName
    :param measFuncInputJac: List of Jacobian matrices, corresponding to the measFuncName list
    :param measFuncInputObj: List of Lists. Each inner list correspond to one measurement function and contains the indices of objects that are part of this measurement. The indice belongs the names list.
    :param measFuncInputObjOpt: List of Lists. Each inner list correspond to one measurement function and contains the indices of objects that are part of the optimization. The indice does not belong to the name list but to the measFuncInputObj list. If it is empty all objects are part of the optimization.
    :param addFunc: If subfunctions are needed these can be added here
    :
    
    """
    for str in addFunc:
        yield str
    
    for i,n,f,jac,opt in zip(range(len(measFuncName)),measFuncName,measFuncInputF,measFuncInputJac,measFuncInputObjOpt):
        yield genOptCodeFuncMeas(n,f,[names[j] for j in measFuncInputObj[i]],\
                                 [variables[j] for j in measFuncInputObj[i]],\
                                 [addVariables[j] for j in measFuncInputObj[i]],debug=debug)
        yield genOptCodeJacMeas(n,jac,[names[j] for j in measFuncInputObj[i]],\
                          [variables[j] for j in measFuncInputObj[i]],\
                          [addVariables[j] for j in measFuncInputObj[i]],\
                          measFuncInputObjOpt[i])
    for n in measFuncName:
        yield genOptCodeFuncMeasPtr(n)
        yield genOptCodeJacMeasPtr(n)

        
        

    
def genOptCodeFuncMeas(name,f,namesOpt,variablesOpt,addVariables,debug=False):
    """
    Method generates C++ code, for the actual measurement.

    :param name: Name of the function
    :param f: function as sympy matrix
    :param namesOpt: Names of the optimization objects, must be identical with the structs
    :param variablesOpt: variables of each optimization object.
    :param addVariables: additional variables of each optimization object need for the function but will not be optimized.


    Returns:

    c++ code which will be placed in meas_func_src.h
    
    Examples:
  
        f,jac=getFandJacProj(mode='metric')
        names=('pt2','pt3','cam','cam_in')
        variables=('x,y','X,Y,Z','CX,CY,CZ,q1,q2,q3,q4','fx,fy,u0,v0,k1,k2,k3')
        addVariables=('','','',['s'])

        genOptCodeFuncMeas('ProjectMetric',f,names,variables,addVariables)
        
    """
    
#    namesOpt=tuple(namesOpt)
 #   variablesOpt=tuple(variablesOpt)
    outStr='__host__ __device__\n'
    outStr+='void %s('%name
    
    #for i in range(1,NR_OPT_PARA+1):
    outStr+='void ** v,'
    #void * _v1, void * _v2,void * _v3, void * _v4,void * _v5
    outStr+='oof_float* func_result){\n'
    for i,no,vs,avs in zip(range(len(namesOpt)),namesOpt,variablesOpt,addVariables):
        #print vs
        
        outStr+='\n'
        outStr+='\t%s_t *v%d=(%s_t *) v[%d];\n'%(no,i+1,no,i)
        outStr+='\n'
        if vs is not None and len(vs)>0:
            var=sympy.symbols(vs)
            for v in var:
                try:
                    free_sym=[]
                    for ff in f:
                        for fff in list(ff.free_symbols):
                            free_sym.append(fff)
                    temp_v=sympy.symbols('v%d_%s'%(i+1,v))
                    if free_sym.index(temp_v)>-1:
                        outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
                except:
                    pass
            
        if avs is not None and len(avs)>0:
            var=sympy.symbols(avs)
            for v in var:
                try:
                    free_sym=[]
                    for ff in f:
                        for fff in list(ff.free_symbols):
                            free_sym.append(fff)
                    temp_v=sympy.symbols('v%d_%s'%(i+1,v))
                    if free_sym.index(temp_v)>-1:                
                        outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
                except:
                    print v, "not used in function", name
                    pass
            
                 
    #print sympy.simplify.cse_main.cse(f)
    
    #ex,eq= sympy.cse(f)
    ex,eq= sympy.cse(f)
    if sympy.__version__>='0.7.2':
        eq=eq[0]
    
    for e in ex:
        outStr+="\toof_float "+str(e[0])+'='+printing.ccode(e[1])+';\n'
        
    for i in range(len(eq)):
        #print "jac_value["+str(i)+"]="+printing.ccode(eq[i])+";"
        
        outStr+="\tfunc_result["+str(i)+"]="+printing.ccode(eq[i])+";\n"
    if debug:
        outStr+='\tprintf("%s:'%name
        for i in range(len(eq)-1):
            outStr+='%f, '
        outStr+='%f\\n"'
        for i in range(len(eq)):
            outStr+=',func_result['+str(i)+']'
        outStr+=');\n'
        
    outStr+='}\n'
    outStr=outStr.replace('pow', 'powf')
    return outStr
    
def genOptCodeJacMeas(name,jac,namesOpt,variablesOpt,addVariables,jacOrder=None):
    """
    Method generates C++ code, for the jacobian matrix.

    :param name: Name of the function
    :param jac: jacobian matrix as sympy matrix
    :param namesOpt: Names of the optimization objects, must be identical with the structs
    :param variablesOpt: variables of each optimization object.
    :param addVariables: additional variables of each optimization object need for the function but will not be optimized.
    :param jacOrder: order of the optimization objets within the jacobian matrix. If None the order will be identical to the order in namesOpt and all object will be optimized. Which object should be optimized can be select (see Example).


    Returns:

    c++ code which will be placed in meas_func_src.h
    
    Example:
  
        f,jac=getFandJacProj(mode='metric')
        names=('pt2','pt3','cam','cam_in')
        variables=('x,y','X,Y,Z','CX,CY,CZ,q1,q2,q3,q4','fx,fy,u0,v0,k1,k2,k3')
        addVariables=('','','',['s'])

        genOptCodeJacMeas('ProjectMetric',jac,names,variables,addVariables,[2,1])
        
    """
    
    
    
    if jacOrder is None or len(jacOrder)==0:
        jacOrder=range(len(namesOpt))
    outStr='__host__ __device__\n'
    outStr+='void %s_jac('%name
    #for i in range(1,NR_OPT_PARA+1):
    outStr+='void ** v,'
     #void * _v2,void * _v3, void * _v4,void * _v5,
    outStr+='int rowStart,oof_float * jac_value, int* jac_row_ind,int* jac_col_ind){\n'
    for i,no,vs,avs in zip(range(len(namesOpt)),namesOpt,variablesOpt,addVariables):
        #print vs
        
        outStr+='\n'
        outStr+='\t%s_t *v%d=(%s_t *) v[%d];\n'%(no,i+1,no,i)
        outStr+='\n'
#        if vs is not None  and len(vs)>0:
#            var=sympy.symbols(vs)
#            for v in var:
#                outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
#            
#        if avs is not None  and len(avs)>0:
#            var=sympy.symbols(avs)
#            for v in var:
#                outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
                
        if vs is not None and len(vs)>0:
            var=sympy.symbols(vs)
            for v in var:
                try:
                    free_sym=[]
                    for ff in jac:
                        for fff in list(ff.free_symbols):
                            free_sym.append(fff)
                    temp_v=sympy.symbols('v%d_%s'%(i+1,v))
                    if free_sym.index(temp_v)>-1:
                        outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
                except:
                    pass
            
        if avs is not None and len(avs)>0:
            var=sympy.symbols(avs)
            for v in var:
                try:
                    free_sym=[]
                    for ff in jac:
                        for fff in list(ff.free_symbols):
                            free_sym.append(fff)
                    temp_v=sympy.symbols('v%d_%s'%(i+1,v))
                    if free_sym.index(temp_v)>-1:                
                        outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
                except:
                    print v, "not used in jacobian",name
                    pass                
                
                
    ex,eq= sympy.cse(jac)
    if sympy.__version__>='0.7.2':
        eq=eq[0]
    for e in ex:
        outStr+="\toof_float "+str(e[0])+'='+printing.ccode(e[1])+';\n'
    

    j=0
    for i in range(len(eq)):
        #print "jac_value["+str(i)+"]="+printing.ccode(eq[i])+";"
        if eq[i]==0:
            continue
        outStr+="\tjac_value["+str(j)+"]="+printing.ccode(eq[i])+";\n"
        rowInd=i/jac.shape[1]
        colInd=i%jac.shape[1]
        count=0
        for ja in jacOrder:
            var=sympy.symbols(variablesOpt[ja])
            if colInd<len(var)+count:
                outStr+='\tjac_row_ind['+str(j)+']=rowStart+%d;\n'%rowInd
                outStr+='\tjac_col_ind['+str(j)+']=v%d->globalStart+%d;\n'%(ja+1,colInd-count)
                break
            else:
                count+=len(var) 
            
        j+=1
        
    outStr+='}\n'
    outStr=outStr.replace('pow', 'powf')
    return outStr    

def genOptCodeJacMeasOld(name,jac,namesOpt,variablesOpt,addVariables,jacOrder=None):
    """
    Method generates C++ code, for the jacobian matrix.

    :param name: Name of the function
    :param jac: jacobian matrix as sympy matrix
    :param namesOpt: Names of the optimization objects, must be identical with the structs
    :param variablesOpt: variables of each optimization object.
    :param addVariables: additional variables of each optimization object need for the function but will not be optimized.
    :param jacOrder: order of the optimization objets within the jacobian matrix. If None the order will be identical to the order in namesOpt and all object will be optimized. Which object should be optimized can be select (see Example).


    Returns:

    c++ code which will be placed in meas_func_src.h
    
    Example:
  
        f,jac=getFandJacProj(mode='metric')
        names=('pt2','pt3','cam','cam_in')
        variables=('x,y','X,Y,Z','CX,CY,CZ,q1,q2,q3,q4','fx,fy,u0,v0,k1,k2,k3')
        addVariables=('','','',['s'])

        genOptCodeJacMeas('ProjectMetric',jac,names,variables,addVariables,[2,1])
        
    """
    
    
    if jacOrder is None:
        jacOrder=range(len(namesOpt))
    outStr='__host__ __device__\n'
    outStr+='void %s_jac('%name
    for i in range(1,NR_OPT_PARA+1):
        outStr+='void * _v%d,'%i
     #void * _v2,void * _v3, void * _v4,void * _v5,
    outStr+='int rowStart,oof_float * jac_value, int* jac_row_ind,int* jac_col_ind){\n'
    for i,no,vs,avs in zip(range(len(namesOpt)),namesOpt,variablesOpt,addVariables):
        #print vs
        var=sympy.symbols(vs)
        outStr+='\n'
        outStr+='\t%s_t *v%d=(%s_t *) _v%d;\n'%(no,i+1,no,i+1)
        outStr+='\n'
        for v in var:
            outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
            
        if avs is not None  and len(avs)>0:
            var=sympy.symbols(avs)
            for v in var:
                outStr+='\toof_float v%d_%s=v%d->%s;\n'%(i+1,v,i+1,v)
    ex,eq= sympy.cse(jac)
    for e in ex:
        outStr+="\toof_float "+str(e[0])+'='+printing.ccode(e[1])+';\n'
        
    for i in range(len(eq)):
        #print "jac_value["+str(i)+"]="+printing.ccode(eq[i])+";"
        outStr+="\tjac_value["+str(i)+"]="+printing.ccode(eq[i])+";\n"
    outStr+='\tint nrPara=0;\n'
    for i in jacOrder:
        outStr+='\tnrPara+=v%d->nrPara;\n'%(i+1)
    
    outStr+='\tfor(int i=0;i<%d;i++)\n'%jac.shape[0]
    outStr+='\t\tfor (int j=0;j<nrPara;j++)\n'
    outStr+='\t\t\tjac_row_ind[j+i*nrPara]=rowStart+i;\n'
    
    outStr+='\tint offset=0;\n'
    
    for i in jacOrder:
        
        outStr+='\tfor(int i=0;i<v%d->nrPara;i++)\n'%(i+1)
        outStr+='\t\tfor(int j=0;j<%d;j++)\n'%jac.shape[0]
        
        #outStr+='\t\tjac_col_ind[i+offset]=v%d->globalStart+i;\n'%(i+1)
        outStr+='\t\t\tjac_col_ind[i+offset+j*nrPara]=v%d->globalStart+i;\n'%(i+1)
        #outStr+='\t}\n'
        outStr+='\toffset+=v%d->nrPara;\n'%(i+1)
    
    
        
    outStr+='}\n'
    outStr=outStr.replace('pow', 'powf')
    return outStr    


def genOptCodeFuncMeasPtr(name):
    """
        Method generates the function pointer on the GPU for the method created in genOptCodeFuncMeas.
        
        :param name: Name of the measurement function.
        
            Returns:

    c++ code which will be placed in meas_func_src.h
    
    Example:

        genOptCodeFuncMeasPtr('ProjectMetric')
    
    """
    
    outStr='__global__ void %s_function_ptr(external_device_func_t* dev_func_ptr){\n'%name
    outStr+='\t*dev_func_ptr = %s;\n'%name
    outStr+='}\n'
    outStr+='\n'
    outStr+=' void %s_function_ptr_h(external_device_func_t* dev_func_ptr){\n'%name
    outStr+='\t*dev_func_ptr = %s;\n'%name
    outStr+='}\n'
    return outStr
    
def genOptCodeJacMeasPtr(name):
    """
        Method generates the function pointer on the GPU for the method created in genOptCodeJacMeas.
        
        :param name: Name of the measurement function.
        
            Returns:

    c++ code which will be placed in meas_func_src.h
    
    Example:

        genOptCodeJacMeasPtr('ProjectMetric')
    
    """
    
    outStr='__global__ void %s_jac_ptr(external_device_func_jac_t* dev_func_ptr){\n'%name
    outStr+='\t*dev_func_ptr = %s_jac;\n'%name
    outStr+='}\n'
    outStr+='\n'
    outStr+='void %s_jac_ptr_h(external_device_func_jac_t* dev_func_ptr){\n'%name
    outStr+='\t*dev_func_ptr = %s_jac;\n'%name
    outStr+='}\n'
    return outStr

def genOptCodeFuncMeasPtrH(name):
    """
        Method generates the function header for the method created in genOptCodeFuncMeasPtr.
        
        :param name: Name of the measurement function.
        
            Returns:

    c++ code which will be placed in meas_func_h.h
    
    Example:

        genOptCodeFuncMeasPtrH('ProjectMetric')
    
    """
    
    
    outStr='__global__ void %s_function_ptr(external_device_func_t* dev_func_ptr);\n'%name
    outStr+='\n'
    outStr+='void %s_function_ptr_h(external_device_func_t* dev_func_ptr);\n'%name
    return outStr
    
def genOptCodeJacMeasPtrH(name):
    """
        Method generates the function header for the method created in genOptCodeJacMeasPtr.
        
        :param name: Name of the measurement function.
        
            Returns:

    c++ code which will be placed in meas_func_h.h
    
    Example:

        genOptCodeJacMeasPtrH('ProjectMetric')
    
    """
    
    outStr='__global__ void %s_jac_ptr(external_device_func_jac_t* dev_func_ptr);\n'%name
    outStr+='\n'
    outStr+='void %s_jac_ptr_h(external_device_func_jac_t* dev_func_ptr);\n'%name
    return outStr

def genOptCodeFuncMeasFile(str,openof_root):     
    """
    Method writes meas_func_src.h.

    :param str: Lists of strings
    
    """
    outStr='#ifndef MEAS_FUNC_SRC_H_\n'
    outStr+='#define MEAS_FUNC_SRC_H_\n'
    outStr+='#include "structs_in.h"\n'
    for s in str:
        outStr+='\n'
        outStr+=s
        outStr+='\n'
    outStr+='#endif'
    f=open(openof_root+'/src/Model/'+'meas_func_src.h','w')
    f.write(outStr)
    f.close
def genOptCodeFuncMeasFileH(str,openof_root):    
    """
    Method writes meas_func_h.h.

    :param str: Lists of strings.
    
    """
 
    outStr='#ifndef MEAS_FUNC_H_H_\n'
    outStr+='#define MEAS_FUNC_H_H_\n'
    for s in str:
        outStr+='\n'
        outStr+=s
        outStr+='\n'
    outStr+='#endif'
    f=open(openof_root+'/src/Model/'+'meas_func_h.h','w')
    f.write(outStr)
    f.close
def genOptCodeFuncMeasHIterator(measFuncName):
    """
    Iterator that can be used for genOptCodeFuncMeasFileH.
    :param measFuncName: List with names of measurements.
    """
    for n in measFuncName:
        yield genOptCodeFuncMeasPtrH(n)
        yield genOptCodeJacMeasPtrH(n)
        
def genOptCodeFillClass(names,namesMeasF,namesMeasJac,optNamesInd,namePerMeasInd,advance=False,name='Minimizer'):
    """
        Method generates the c++ code that is necessary to fill the optimization with values.
        
        :param names: List of strings containing the names of the optimization objects.
        :param namesMeasF: List of strigns containing the names of the measurement functions.
        :param namesMeasJac: List of strigns containing the names of the function computing the jacobie matrix.
        :param optNamesInd: List of indizes of the optimization objects which will be optimized.
        :param namePerMeasInd: List of List of the index of the optimization objects which will go in a measurement. The outer list has the size of the number of measurement functions. 
        :param advance: If True, the minimize method has two extra attributes (kmax: maximum number of iterations of LM, residuals: vector that contains the residual of each iteration).
        :param name: default=minimize, name of the method.
        
        Returns:

    c++ code which will be placed in e.g. optimize.cu
    
    Example:
    
    """

    outStr='%s::%s(){\n'%(name,name)
    outStr+='\tverbose=1;\n'
    outStr+='\tkmax=100;\n'
    outStr+='\teps1=1e-10;\n'
    outStr+='\teps2=1e-10;\n'

    outStr+='\tcg_it=150;\n'
    outStr+='\tcg_thresh_rel=1e-12;\n'
    outStr+='\tcg_thresh_abs=1e-12;\n'
    outStr+='}\n'

    outStr+='void %s::run(){\n'%name

    outStr+='\tLevmar myLevmar;\n'
    if advance:
        outStr+='\tmyLevmar.kmax=kmax;\n'
    outStr+='\tmyLevmar.verbose=verbose;\n'
    outStr+='\tmyLevmar.kmax=kmax;\n'
    outStr+='\tmyLevmar.eps1=eps1;\n'
    outStr+='\tmyLevmar.eps2=eps2;\n'
    outStr+='\tmyLevmar.cg_it=cg_it;\n'
    outStr+='\tmyLevmar.cg_thresh_rel=cg_thresh_rel;\n'
    outStr+='\tmyLevmar.cg_thresh_abs=cg_thresh_abs;\n'
#    outStr+='\tmyLevmar.cg_thresh2=cg_thresh2;\n\n'

    outStr+='\tcudaError_t status;\n'
    outStr+='\n'
    outStr+='#ifdef USE_GPU\n'
    for n in namesMeasF:
            outStr+='\texternal_device_func_t * func%s;\n'%n
            outStr+='\tstatus = cudaMalloc((void**)&func%s, sizeof(external_device_func_t));\n'%n
            outStr+='\tassert(status == cudaSuccess);\n'
            outStr+='\t%s_function_ptr<<<1, 1>>>(func%s);\n'%(n,n)
            outStr+='\n'
    
    for n in namesMeasJac:
            outStr+='\texternal_device_func_jac_t * jac%s;\n'%n
            outStr+='\tstatus = cudaMalloc((void**)&jac%s, sizeof(external_device_func_jac_t));\n'%n
            outStr+='\tassert(status == cudaSuccess);\n'
            outStr+='\t%s_jac_ptr<<<1, 1>>>(jac%s);\n'%(n,n)
            outStr+='\n'
    
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\texternal_device_set_value_func_t * set_value_%s;\n'%n
        outStr+='\tstatus = cudaMalloc((void**)&set_value_%s, sizeof(external_device_set_value_func_t));\n'%n
        outStr+='\tassert(status == cudaSuccess);\n'
        outStr+='\tset_value_%s_ptr<<<1, 1>>>(set_value_%s);\n'%(n,n)
        outStr+='\n' 
    outStr+='#else\n'
    for n in namesMeasF:
            outStr+='\texternal_device_func_t * func%s;\n'%n
            outStr+='\tfunc%s = (external_device_func_t *) malloc ( sizeof(external_device_func_t) );\n'%n
            outStr+='\tassert(func%s != NULL);\n'%n
            #outStr+='\t*func%s=%s;\n'%(n,n)
            outStr+='\t%s_function_ptr_h(func%s);\n'%(n,n)
            outStr+='\n'
    
    for n in namesMeasJac:
            outStr+='\texternal_device_func_jac_t * jac%s;\n'%n
            outStr+='\tjac%s = (external_device_func_jac_t *) malloc ( sizeof(external_device_func_jac_t) );\n'%n
            outStr+='\tassert(jac%s != NULL);\n'%n
            #outStr+='\t*jac%s=%s_jac;\n'%(n,n)
            outStr+='\t%s_jac_ptr_h(jac%s);\n'%(n,n)
            outStr+='\n'
    
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\texternal_device_set_value_func_t * set_value_%s = NULL;\n'%n
        outStr+='\tset_value_%s = (external_device_set_value_func_t *) malloc ( sizeof(external_device_set_value_func_t));\n'%n
        outStr+='\tassert(NULL != set_value_%s);\n'%n
        #outStr+='\t*set_value_%s = set_value_%s_func;\n'%(n,n)
        outStr+='\tset_value_%s_ptr_h(set_value_%s);\n'%(n,n)
        outStr+='\n'           
    outStr+='#endif\n'            
      
        
    outStr+='\n'
    for n in names:
        outStr+='\tthrust::host_vector<%s_t> h_%s=%s;\n'%(n,n,n)
    outStr+='\n'        
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tfor(int i=0;i<h_%s.size();i++)\n'%n
        outStr+='\t\th_%s[i].globalStart=myLevmar.registerOptObj(h_%s[i].nrPara);\n'%(n,n)
    outStr+='\n'        
    for n in names:
        outStr+='\tmemory_vector<%s_t> d_%s=h_%s;\n'%(n,n,n)    
    outStr+='\n'
    outStr+='\tthrust::host_vector<MeasurementStruct_t> h_MeasVec(measurementCombinations.size());\n'
    outStr+='\n'
    outStr+='\tint nr_pointer=0;\n'
    outStr+='\tfor(int i=0;i<measurementCombinations.size();i++){\n'
    outStr+='\t\tnr_pointer+=measurementCombinations[i].v.size();\n'
    outStr+='\t}\n'
    outStr+='\tint current_pointer_index=0;\n'
    
    outStr+='\tvoid **v_temp=(void**)malloc(nr_pointer * sizeof(void*));\n'
    outStr+='\tvoid **v;\n'
    outStr+='#ifdef USE_GPU\n'
    outStr+='\tcudaMalloc((void**)&v, nr_pointer*sizeof(void*));\n'
    outStr+='#else\n'
    outStr+='\tv=v_temp;\n'
    outStr+='#endif\n'
    outStr+='\tfor(int i=0;i<measurementCombinations.size();i++){\n'
    outStr+='\t\tint type=measurementCombinations[i].type;\n'
    for nf,nj,ind,i in zip(namesMeasF,namesMeasJac,namePerMeasInd,range(len(namesMeasF))):
        outStr+='\t\tif (type==%d){\n'%i
        outStr+='\t\t\tMeasurement_%s_t m(func%s,jac%s,&v[current_pointer_index]);\n'%(nf,nf,nj)
        for nameInd,j in zip(ind,range(len(ind))):
            outStr+='\t\t\t\tv_temp[current_pointer_index+%d]=thrust::raw_pointer_cast(&d_%s[measurementCombinations[i].v[%d]]);\n'%(j,names[nameInd],j)
        
        outStr+='\t\t\tm.robust=measurementCombinations[i].robust;\n'
        outStr+='\t\t\tm.robust_para_a=measurementCombinations[i].robust_para_a;\n'
        outStr+='\t\t\tm.robust_para_b=measurementCombinations[i].robust_para_b;\n'
        outStr+='\t\t\tm.robust_para_c=measurementCombinations[i].robust_para_c;\n'
        outStr+='\t\t\tmyLevmar.registerMeasObj(&m);\n'
        outStr+='\t\t\th_MeasVec[i]=m;\n'
        outStr+='\t\t}\n'
    outStr+='\t\tcurrent_pointer_index+=measurementCombinations[i].v.size();\n'
    outStr+='\t}\n'  
    outStr+='#ifdef USE_GPU\n'
    outStr+='cudaMemcpy(v,v_temp,nr_pointer*sizeof(void*),cudaMemcpyHostToDevice);\n'
    outStr+='#endif\n'        
    outStr+='\tmyLevmar.initMemory();\n\n'

    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tfor(int i=0;i<h_%s.size();i++)\n'%n
        outStr+='\t\tmyLevmar.insertOptObj(thrust::raw_pointer_cast(&d_%s[i]),set_value_%s,h_%s[i].nrPara,h_%s[i].globalStart);\n'%(n,n,n,n)          
    outStr+='\n'
    outStr+='\tfor (int i=0;i<h_MeasVec.size();i++){\n'
    outStr+='\t\tmyLevmar.insertMeasObj(&h_MeasVec[i],measurementCombinations[i].inv_cov);\n'
    outStr+='\t}\n'                                                    
    outStr+='\n'
    outStr+='\tmyLevmar.init();\n'
    outStr+='\tmyLevmar.run();\n'
    outStr+='\n'
    outStr+='\tfor (int i=0;i<h_MeasVec.size();i++)'
    outStr+='\t\tfor (int j=0;j<h_MeasVec[i].nr_row;j++)'
    outStr+='\t\t\tmeasurementCombinations[i].res.push_back(h_MeasVec[i].h_func_result[j]);'

    outStr+='\n'
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tthrust::copy(d_%s.begin(),d_%s.end(),%s.begin());\n'%(n,n,n)
        #outStr+='\tstd::copy(h_%s.begin(),h_%s.end(),%s.begin());\n'%(n,n,n)
        outStr+='\n'
    outStr+='\n'
    #if advance:
    outStr+='\tresiduals.resize(myLevmar.residuals.size());\n'
    outStr+='\tstd::copy(myLevmar.residuals.begin(),myLevmar.residuals.end(),residuals.begin());\n'
    
#    outStr+='\tif (calcVar){\n'
#    outStr+='\t\tmyLevmar.copyVarToStruct();\n'
#    for ind in optNamesInd:
#        n=names[ind]
#        outStr+='\t\tvar_%s.resize(%s.size());'%(n,n)
#        outStr+='\t\tthrust::copy(d_%s.begin(),d_%s.end(),var_%s.begin());\n'%(n,n,n)
        #outStr+='\tstd::copy(h_%s.begin(),h_%s.end(),%s.begin());\n'%(n,n,n)
#        outStr+='\n'
#    outStr+='\t}\n'
    
    outStr+='#ifdef USE_GPU\n'
    outStr+='\tstatus = cudaFree(v);\n'
    outStr+='\tdelete [] v_temp;\n'
    for n in namesMeasF:
        outStr+='\tstatus = cudaFree(func%s);\n'%n
    for n in namesMeasJac:
        outStr+='\tstatus = cudaFree(jac%s);\n'%n
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tstatus = cudaFree(set_value_%s);\n'%n
    outStr+='#else\n'
    for n in namesMeasF:
        outStr+='\tfree(func%s);\n'%n
    for n in namesMeasJac:
        outStr+='\tfree(jac%s);\n'%n
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tfree(set_value_%s);\n'%n
    outStr+='#endif\n'    
        
    
    outStr+='}\n'
    outStr+='\n'
    outStr+='%s::~%s(){\n'%(name,name)
    
    outStr+='}\n'
    return outStr


    
def genOptCodeFillClassH(names,namesMeasF,namesMeasJac,optNamesInd,namePerMeasInd,advance=False,name='Minimizer'):
    """
        Method generates the c++ code that is necessary to fill the optimization with values.
        
        :param names: List of strings containing the names of the optimization objects.
        :param namesMeasF: List of strigns containing the names of the measurement functions.
        :param namesMeasJac: List of strigns containing the names of the function computing the jacobie matrix.
        :param optNamesInd: List of indizes of the optimization objects which will be optimized.
        :param namePerMeasInd: List of List of the index of the optimization objects which will go in a measurement. The outer list has the size of the number of measurement functions. 
        :param advance: If True, the minimize method has two extra attributes (kmax: maximum number of iterations of LM, residuals: vector that contains the residual of each iteration).
        :param name: default=minimize, name of the method.
        
        Returns:

    c++ code which will be placed in e.g. optimizate.h
    
    Example:
    
    """
    
    outStr='class %s{\n\n'%name
    outStr+='public:\n'
    outStr+='\tint kmax;\n'
    outStr+='\tint verbose;\n'
    outStr+='\toof_float eps1;\n'
    outStr+='\toof_float eps2;\n'
    outStr+='\tint cg_it;\n'
    outStr+='\toof_float cg_thresh_rel;\n'
    outStr+='\toof_float cg_thresh_abs;\n'
#    outStr+='\tbool calcVar;\n'
    
    for n in names:
        outStr+='\tstd::vector<%s_t> %s;\n'%(n,n)
#        outStr+='\tstd::vector<%s_t> var_%s;\n'%(n,n)
    
    outStr+='\tstd::vector<MeasurementCombinations_t> measurementCombinations;\n'
    outStr+='\tstd::vector<double> residuals;\n'    
    outStr+='\n\n'
    outStr+='\t %s();\n'%name
    outStr+='\t void run();\n'
    outStr+='\t ~%s();\n'%name
    
    outStr+='};\n'
    return outStr

def genOptCodeFill(names,namesMeasF,namesMeasJac,optNamesInd,namePerMeasInd,advance=False,name='minimize'):
    """
        Method generates the c++ code that is necessary to fill the optimization with values.
        
        :param names: List of strings containing the names of the optimization objects.
        :param namesMeasF: List of strigns containing the names of the measurement functions.
        :param namesMeasJac: List of strigns containing the names of the function computing the jacobie matrix.
        :param optNamesInd: List of indizes of the optimization objects which will be optimized.
        :param namePerMeasInd: List of List of the index of the optimization objects which will go in a measurement. The outer list has the size of the number of measurement functions. 
        :param advance: If True, the minimize method has two extra attributes (kmax: maximum number of iterations of LM, residuals: vector that contains the residual of each iteration).
        :param name: default=minimize, name of the method.
        
        Returns:

    c++ code which will be placed in optimization.cu
    
    Example:

        genOptCodeOptimizeFile(genOptCodeFill(names,['ProjectMetric','ProjectFocal','Quat'],['ProjectMetric','ProjectFocal','Quat'],[3,2,1],[[0,1,2,3],[0,1,2,3],[2]]))
    
    """

    outStr='void %s('%name

    for n in names:
        outStr+='std::vector<%s_t> &%s,'%(n,n)
    if advance:    
        outStr+='std::vector<MeasurementCombinations_t> &measurementCombinations,int kmax,std::vector<double> &residuals)\n{\n'
    else:
        outStr+='std::vector<MeasurementCombinations_t> &measurementCombinations)\n{\n'
    outStr+='\tLevmar myLevmar;\n'
    if advance:
        outStr+='\tmyLevmar.kmax=kmax;\n'
    outStr+='\tmyLevmar.verbose=1;\n'
    outStr+='\tcudaError_t status;\n'
    outStr+='\n'
    outStr+='#ifdef USE_GPU\n'
    for n in namesMeasF:
            outStr+='\texternal_device_func_t * func%s;\n'%n
            outStr+='\tstatus = cudaMalloc((void**)&func%s, sizeof(external_device_func_t));\n'%n
            outStr+='\tassert(status == cudaSuccess);\n'
            outStr+='\t%s_function_ptr<<<1, 1>>>(func%s);\n'%(n,n)
            outStr+='\n'
    
    for n in namesMeasJac:
            outStr+='\texternal_device_func_jac_t * jac%s;\n'%n
            outStr+='\tstatus = cudaMalloc((void**)&jac%s, sizeof(external_device_func_jac_t));\n'%n
            outStr+='\tassert(status == cudaSuccess);\n'
            outStr+='\t%s_jac_ptr<<<1, 1>>>(jac%s);\n'%(n,n)
            outStr+='\n'
    
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\texternal_device_set_value_func_t * set_value_%s;\n'%n
        outStr+='\tstatus = cudaMalloc((void**)&set_value_%s, sizeof(external_device_set_value_func_t));\n'%n
        outStr+='\tassert(status == cudaSuccess);\n'
        outStr+='\tset_value_%s_ptr<<<1, 1>>>(set_value_%s);\n'%(n,n)
        outStr+='\n' 
    outStr+='#else\n'
    for n in namesMeasF:
            outStr+='\texternal_device_func_t * func%s;\n'%n
            outStr+='\tfunc%s = (external_device_func_t *) malloc ( sizeof(external_device_func_t) );\n'%n
            outStr+='\tassert(func%s != NULL);\n'%n
            #outStr+='\t*func%s=%s;\n'%(n,n)
            outStr+='\t%s_function_ptr_h(func%s);\n'%(n,n)
            outStr+='\n'
    
    for n in namesMeasJac:
            outStr+='\texternal_device_func_jac_t * jac%s;\n'%n
            outStr+='\tjac%s = (external_device_func_jac_t *) malloc ( sizeof(external_device_func_jac_t) );\n'%n
            outStr+='\tassert(jac%s != NULL);\n'%n
            #outStr+='\t*jac%s=%s_jac;\n'%(n,n)
            outStr+='\t%s_jac_ptr_h(jac%s);\n'%(n,n)
            outStr+='\n'
    
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\texternal_device_set_value_func_t * set_value_%s = NULL;\n'%n
        outStr+='\tset_value_%s = (external_device_set_value_func_t *) malloc ( sizeof(external_device_set_value_func_t));\n'%n
        outStr+='\tassert(NULL != set_value_%s);\n'%n
        #outStr+='\t*set_value_%s = set_value_%s_func;\n'%(n,n)
        outStr+='\tset_value_%s_ptr_h(set_value_%s);\n'%(n,n)
        outStr+='\n'           
    outStr+='#endif\n'            
      
        
    outStr+='\n'
    for n in names:
        outStr+='\tthrust::host_vector<%s_t> h_%s=%s;\n'%(n,n,n)
    outStr+='\n'        
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tfor(int i=0;i<h_%s.size();i++)\n'%n
        outStr+='\t\th_%s[i].globalStart=myLevmar.registerOptObj(h_%s[i].nrPara);\n'%(n,n)
    outStr+='\n'        
    for n in names:
        outStr+='\tmemory_vector<%s_t> d_%s=h_%s;\n'%(n,n,n)    
    outStr+='\n'
    outStr+='\tthrust::host_vector<MeasurementStruct_t> h_MeasVec(measurementCombinations.size());\n'
    outStr+='\n'
    outStr+='\tint nr_pointer=0;\n'
    outStr+='\tfor(int i=0;i<measurementCombinations.size();i++){\n'
    outStr+='\t\tnr_pointer+=measurementCombinations[i].v.size();\n'
    outStr+='\t}\n'
    outStr+='\tint current_pointer_index=0;\n'
    outStr+='\tvoid **v_temp=(void**)malloc(nr_pointer * sizeof(void*));\n'
    outStr+='\tvoid **v;\n'
    outStr+='#ifdef USE_GPU\n'
    outStr+='\tcudaMalloc((void**)&v, nr_pointer*sizeof(void*));\n'
    outStr+='#else\n'
    outStr+='\tv=v_temp;\n'
    outStr+='#endif\n'
    outStr+='\tfor(int i=0;i<measurementCombinations.size();i++){\n'
    outStr+='\t\tint type=measurementCombinations[i].type;\n'
    for nf,nj,ind,i in zip(namesMeasF,namesMeasJac,namePerMeasInd,range(len(namesMeasF))):
        outStr+='\t\tif (type==%d){\n'%i
        outStr+='\t\t\tMeasurement_%s_t m(func%s,jac%s,&v[current_pointer_index]);\n'%(nf,nf,nj)
        for nameInd,j in zip(ind,range(len(ind))):
            outStr+='\t\t\t\tv_temp[current_pointer_index+%d]=thrust::raw_pointer_cast(&d_%s[measurementCombinations[i].v[%d]]);\n'%(j,names[nameInd],j)
        
        outStr+='\t\t\tm.robust=measurementCombinations[i].robust;\n'
        outStr+='\t\t\tm.robust_para_a=measurementCombinations[i].robust_para_a;\n'
        outStr+='\t\t\tm.robust_para_b=measurementCombinations[i].robust_para_b;\n'
        outStr+='\t\t\tm.robust_para_c=measurementCombinations[i].robust_para_c;\n'
        outStr+='\t\t\tmyLevmar.registerMeasObj(&m);\n'
        outStr+='\t\t\th_MeasVec[i]=m;\n'
        outStr+='\t\t}\n'
    outStr+='\t\tcurrent_pointer_index+=measurementCombinations[i].v.size();\n'
    outStr+='\t}\n'  
    
            
    outStr+='#ifdef USE_GPU\n'
    outStr+='cudaMemcpy(v,v_temp,nr_pointer*sizeof(void*),cudaMemcpyHostToDevice);\n'
    outStr+='#endif\n'
    outStr+='\tmyLevmar.initMemory();\n\n'

    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tfor(int i=0;i<h_%s.size();i++)\n'%n
        outStr+='\t\tmyLevmar.insertOptObj(thrust::raw_pointer_cast(&d_%s[i]),set_value_%s,h_%s[i].nrPara,h_%s[i].globalStart);\n'%(n,n,n,n)          
    outStr+='\n'
    outStr+='\tfor (int i=0;i<h_MeasVec.size();i++){\n'
    outStr+='\t\tmyLevmar.insertMeasObj(&h_MeasVec[i],measurementCombinations[i].inv_cov);\n'
    outStr+='\t}\n'                                                    
    outStr+='\n'
    outStr+='\tmyLevmar.init();\n'
    outStr+='\tmyLevmar.run();\n'
    outStr+='\n'
    outStr+='\tfor (int i=0;i<h_MeasVec.size();i++)'
    outStr+='\t\tfor (int j=0;j<h_MeasVec[i].nr_row;j++)'
    outStr+='\t\t\tmeasurementCombinations[i].res.push_back(h_MeasVec[i].h_func_result[j]);'
    outStr+='\n'
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tthrust::copy(d_%s.begin(),d_%s.end(),%s.begin());\n'%(n,n,n)
        #outStr+='\tstd::copy(h_%s.begin(),h_%s.end(),%s.begin());\n'%(n,n,n)
        outStr+='\n'
    outStr+='\n'
    if advance:
        outStr+='\tresiduals.resize(myLevmar.residuals.size());\n'
        outStr+='\tstd::copy(myLevmar.residuals.begin(),myLevmar.residuals.end(),residuals.begin());\n'

    outStr+='\tdelete [] v_temp;\n'
    outStr+='#ifdef USE_GPU\n'
    outStr+='\tstatus = cudaFree(v);\n'
    
    for n in namesMeasF:
        outStr+='\tstatus = cudaFree(func%s);\n'%n
    for n in namesMeasJac:
        outStr+='\tstatus = cudaFree(jac%s);\n'%n
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tstatus = cudaFree(set_value_%s);\n'%n
    outStr+='#else\n'
    for n in namesMeasF:
        outStr+='\tfree(func%s);\n'%n
    for n in namesMeasJac:
        outStr+='\tfree(jac%s);\n'%n
    for ind in optNamesInd:
        n=names[ind]
        outStr+='\tfree(set_value_%s);\n'%n
    outStr+='#endif\n'    
        
    
    outStr+='}'
    return outStr
def genOptCodeOptimizeFile(src,name="optimize",openof_root=None):
    """
    Method writes optimizer.cu.

    :param str: string containing the optimization object
    
    """
 
    outStr='#include "%s.h"\n'%name
    outStr+='#include <thrust/host_vector.h>\n'
    outStr+='#ifdef USE_GPU\n'
    outStr+='\t#include <thrust/system/cuda/vector.h>\n'
    outStr+='#else\n'
    outStr+='\t#include <thrust/system/tbb/vector.h>\n'
    outStr+='#endif\n'
    outStr+='#include <iostream>\n'
    outStr+='#include <Levmar.h>\n'
    outStr+='#include <assert.h>\n'
    #outStr+='#include "structs_func_src.h"\n'
    #outStr+='#include "meas_func_src.h"\n'
    outStr+='#include "structs_meas.h"\n'
    outStr+='#include <vector>\n'
    outStr+='#include <cstdlib>\n'
    outStr+='\n'
    if type(src)==list:
        for s in src:
            outStr+='\n'    
            outStr+=s
            outStr+='\n'
    else:
        outStr+=src

    f=open(openof_root+'/src/Model/'+name+'.cu','w')
    f.write(outStr)
    f.close
        
def genOptCodeOptimizeFileH(namesMeasF,src,name="optimize",filenameStruct='structs',openof_root=None):
    """
    Method writes optimizer.h.

    :param str: string containing the optimization object
    
    """
    
    if type(src)==list:
        sp=src[0].splitlines()
        s= sp[0]+';'
        sclass=src[1]
    else:
        sp=src.splitlines()
        s= sp[0]+';'
        sclass=''
    outStr='#ifndef '+name.upper()+'_H_\n'
    outStr+='#define '+name.upper()+'_H_\n'
    outStr+='#include <vector>\n'
    outStr+='#include "%s.h"\n'%filenameStruct

    outStr+='\n'
    outStr+='enum MeasTyp {'
    for i,n in zip(range(len(namesMeasF)),namesMeasF):
        outStr+='e'+n
        if i<len(namesMeasF)-1:
            outStr+=','
    outStr+='};\n\n'
            
    
    outStr+=s
    outStr+='\n\n\n'
    outStr+=sclass
    outStr+='\n'
    outStr+='#endif'
    f=open(openof_root+'/src/Model/'+name+'.h','w')
    f.write(outStr)
    f.close
def genOptCodeExtractHeader(namesMeasF,src):
    """
    Method extracts the header of the minimize function and created an enumeration for the measurement functions.
    :param neamsMeasF: names of the measurement functions
    :param src: source of the method that fills the optimization with data.    
    """
    outStr='enum MeasTyp {'
    for i,n in zip(range(len(namesMeasF)),namesMeasF):
        outStr+='e'+n
        if i<len(namesMeasF):
            outStr+=','
    outStr+='};\n\n'

    if type(src)==list:
        for s in src:
            outStr+='\n'
            if s.find('class')>=0:
                outStr+=s
            else:
                sp=s.splitlines()
                outStr+= sp[0]+';'
            outStr+='\n'
    else:            
        sp=src.splitlines()
        outStr+= sp[0]+';'
    return outStr      
def genOptCodePython(names,nameLib,filenameOpt,fct,openof_root):
    """
    Method that generates the python interfaces for SWIG.
    :param names: Names of the objects.
    :param nameLib: name of the library that will be generated
    :param filenameOpt: name of the file generated in genOptCodeOptimizeFile
    :param fct: header of the function that starts the optimization
    """
    outStr='%%module %s\n'%nameLib
    outStr+='%{\n'
    outStr+=' /* Includes the header in the wrapper code */\n'
    outStr+='#define SWIG_FILE_WITH_INIT\n'
    outStr+='#include <vector>\n'
    outStr+='#include "levmar_structs_opt.h"\n'
    outStr+='#include "structs_in.h"\n'
    outStr+='#include "%s.h"\n'%filenameOpt
    outStr+='%}\n'
    
    outStr+='%include "std_vector.i"\n'
    outStr+='%include "levmar_structs_opt.h"\n'
    outStr+='%include "structs_in.h"\n'
    outStr+='\n'
    for n in names:
        outStr+='%%template(%sVector) std::vector<%s_t>;\n'%(n,n)
    outStr+='%template(MeasurementCombinationsVector) std::vector<MeasurementCombinations_t>;\n'
    outStr+='%template(doubleVector) std::vector<double>;\n'
    outStr+='%template(floatVector) std::vector<float>;\n'
    outStr+='%template(intVector) std::vector<int>;\n'
    outStr+='\n'
    outStr+=fct
    f=open(openof_root+'/src/Model/'+filenameOpt+'.i','w')
    f.write(outStr)
    f.close 
def genOptCodeCMakeLists(libName,libNamePy,filenameOpt,openof_root):
    """
    creates the CMakeLists.txt for library.
    :param libName: name of the C++ library.
    :param libNamePy: name of the Python library.
    :param filenameOpt: name of the file that contains the minimize method.
    """
    
    fin=open(openof_root+'/Generator/CMakeLists.tpl','r')
    fout=open(openof_root+'CMakeLists.txt','w')
    str=fin.read()
    strOut=str%(libName,libName,filenameOpt,libName,filenameOpt,filenameOpt,libNamePy,filenameOpt,filenameOpt,libNamePy,libName,libNamePy)
    fout.write(strOut)
    fin.close()
    fout.close()
    
def genAll(names,variables,addVariables,addCode,measFuncName,measFuncInputF,measFuncInputObj,measFuncInputObjOpt,measFuncObjOpt,addFunc=[],filename_structs='structs',min_function_name='minimize',min_class_name='Minimizer',filename_cpp='optimize',libname_python='OptimizerPy',libname_cpp='Optimizer',debug=False,clear_model_folder=False,openof_root='../'):
    measFuncInputJac=[]
    for i,f,o in zip(range(len(measFuncInputF)),measFuncInputF,measFuncInputObjOpt):
        var=''
        for vind in o:
            vaddon='v%d_'%(vind+1)
            print vind
            print variables[measFuncInputObj[i][vind]]
            if type(variables[measFuncInputObj[i][vind]])==list:
                varstr=vaddon+variables[measFuncInputObj[i][vind]][0]
            else:
                varstr=vaddon+variables[measFuncInputObj[i][vind]]
            varstr=varstr.replace(',',','+vaddon)
            if var=='':
                var+=varstr
            else:
                var+=','+varstr
        print var
	if len(var.split(','))==1:
		var_sympy=sympy.Matrix([sympy.symbols(var)])
		print var_sympy
	else:
	        var_sympy=sympy.Matrix(sympy.symbols(var))
        
        
        jac=f.jacobian(var_sympy)
        measFuncInputJac.append(jac)
    
    if clear_model_folder:
        for f in glob.glob(openof_root+'/src/Model/*'):
            if f.find('README')<0:
	            os.remove(f)
    genOptCodeStructFuncFile(names,variables,openof_root=openof_root)
    genOptCodeStructFile(names,variables,addVariables,addCode,filename=filename_structs,openof_root=openof_root)
    
    genOptCodeMeasFile((str for str in genOptCodeMeasIterator(measFuncName,names,measFuncInputJac,measFuncInputObj)),openof_root=openof_root)
    genOptCodeFuncMeasFile((str for str in genOptCodeFuncMeasIterator(measFuncName, names, variables, addVariables, measFuncInputF, measFuncInputJac, measFuncInputObj, measFuncInputObjOpt,addFunc,debug=debug)),openof_root=openof_root)
    genOptCodeFuncMeasFileH((str for str in genOptCodeFuncMeasHIterator(measFuncName)),openof_root=openof_root)                                                                     
    
    
    genOptCodeOptimizeFile([genOptCodeFill(names,measFuncName,measFuncName,measFuncObjOpt,measFuncInputObj,advance=True,name=min_function_name)\
                            ,genOptCodeFillClass(names,measFuncName,measFuncName,measFuncObjOpt,measFuncInputObj,name=min_class_name)],name=filename_cpp,openof_root=openof_root)
    genOptCodeOptimizeFileH(measFuncName,[genOptCodeFill(names,measFuncName,measFuncName,measFuncObjOpt,measFuncInputObj,advance=True,name=min_function_name)\
                                          ,genOptCodeFillClassH(names,measFuncName,measFuncName,measFuncObjOpt,measFuncInputObj,name=min_class_name)],name=filename_cpp,filenameStruct=filename_structs,openof_root=openof_root)
    header=genOptCodeExtractHeader(measFuncName,[genOptCodeFill(names,measFuncName,measFuncName,measFuncObjOpt,measFuncInputObj,advance=True,name=min_function_name)\
                                          ,genOptCodeFillClassH(names,measFuncName,measFuncName,measFuncObjOpt,measFuncInputObj,name=min_class_name)])
    genOptCodePython(names,libname_python,filename_cpp,header,openof_root=openof_root)
    genOptCodeCMakeLists(libname_cpp,libname_python,filename_cpp,openof_root=openof_root)
    
#    shutil.copy('CMakeLists.txt',openof_root)
#    os.remove('CMakeLists.txt')

        
#    for f in glob.glob('*.h'):
#        shutil.copy(f,openof_root+'/src/Model/')
#        os.remove(f)
#    for f in glob.glob('*.cu'):
#        shutil.copy(f,openof_root+'/src/Model/')
#        os.remove(f)
#    for f in glob.glob('*.i'):
#        shutil.copy(f,openof_root+'/src/Model/')
#        os.remove(f)   
    
 
    
