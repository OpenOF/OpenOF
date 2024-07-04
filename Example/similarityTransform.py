'''

 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *

 
'''

import sys
sys.path.append('../Generator')
import Quaternion
import sympy
from genOptStruct import genAll

def getFSimTrans():
    
    XS,YS,ZS=sympy.symbols('v1_X,v1_Y,v1_Z')
    XD,YD,ZD=sympy.symbols('v2_X,v2_Y,v2_Z')
    X,Y,Z,q1,q2,q3,q4,s=sympy.symbols('v3_CX,v3_CY,v3_CZ,v3_q1,v3_q2,v3_q3,v3_q4,v3_s')
    
    
    q=Quaternion.Quaternion()
    
    q.setValues(q1,q2,q3,q4)
    
    pt3=Quaternion.Quaternion()

    pt3.setValues(0,XS,YS,ZS)
    
    pt3d=s*q*pt3*q.conj()
    
    f=sympy.Matrix([pt3d[1]-XD+X,pt3d[2]-YD+Y,pt3d[3]-ZD+Z])
    
    return f

def getFQuat():
    
    X,Y,Z,q1,q2,q3,q4,s=sympy.symbols('v1_CX,v1_CY,v1_CZ,v1_q1,v1_q2,v1_q3,v1_q4,v1_s')
    q=Quaternion.Quaternion()
    

    f=(1-(q1**2+q2**2+q3**2+q4**2))**2
       
    fm=sympy.Matrix([f])
    
    return fm


if __name__ == '__main__':
    print(sympy.__version__)
    
    '''define the error/cost function of the similarity transformation'''
    f1=getFSimTrans()
    '''define the error/cost function of the quaternion, forcing the quaternion to be of unit length'''
    f2=getFQuat()
    
    '''If enabled a print statement is included in each function'''
    debug=False
    
    '''names of the different optimization objects'''
    names=['pt3','transform']
    '''variable parameters of the optimization objects''' 
    variables=['X,Y,Z','CX,CY,CZ,q1,q2,q3,q4,s']
    '''additional non variable parameters'''
    addVariables=('','')
    ''' additional c++ code which is included in the struct belonging to the optimization object'''
    addCode=('','')

    '''Give each function a unique name'''
    measFuncName=['SimilarityTransformation','Quat']
    '''Include the definition of the the cost funciton'''
    measFuncInputF=[f1,f2]
    '''define the input structs for each error/cost function e.g. [0,0,1] means [pt,pt,trans]'''
    measFuncInputObj=[[0,0,1],[1]]
    '''define which objects are estimated for each error/cost function, e.g. [2] means that for the SimilarityTransformation function only the values of trans are allowed to be changed.'''
    measFuncInputObjOpt=[[2],[0]]
    '''defines which input structs in the overall optimization are allowed to change.'''
    measFuncObjOpt=[1]

    '''generate the complete code'''
    genAll(names,\
           variables,\
           addVariables,\
           addCode,\
           measFuncName,\
           measFuncInputF,\
           measFuncInputObj,\
           measFuncInputObjOpt,\
           measFuncObjOpt,\
           filename_structs='structs_simtrans',\
           min_function_name='minimize_simtrans',\
           min_class_name='MinimizerSimTrans',\
           filename_cpp='optimizeSimTrans',\
           libname_python='OptimizerSimTransPy',\
           libname_cpp='OptimizerSimTrans',\
           clear_model_folder=True,\
           openof_root='../')

