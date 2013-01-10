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
 
'''


import sys
sys.path.append('../Generator')
import Quaternion
import sympy
from genOptStruct import genAll

def getFProj():
    
    q=Quaternion.Quaternion()
    q1=sympy.symbols('v3_q1')
    q2=sympy.symbols('v3_q2')
    q3=sympy.symbols('v3_q3')
    q4=sympy.symbols('v3_q4')
    
    q.setValues(q1,q2,q3,q4)
    
    pt3=Quaternion.Quaternion()
    
    X=sympy.symbols('v2_X')
    Y=sympy.symbols('v2_Y')
    Z=sympy.symbols('v2_Z')
    
    pt3.setValues(0.,X,Y,Z)
    
    CX,CY,CZ=sympy.symbols('v3_CX,v3_CY,v3_CZ')
    
    pt3[1]=pt3[1]-CX
    pt3[2]=pt3[2]-CY
    pt3[3]=pt3[3]-CZ
    pt3c2=q*pt3*q.conj()
    
    fx,fy,u0,v0,s,k1,k2,k3=sympy.symbols('v4_fx,v4_fy,v4_u0,v4_v0,v4_s,v4_k1,v4_k2,v4_k3')
    
    
    ptx=pt3c2[1]/pt3c2[3]
    pty=pt3c2[2]/pt3c2[3]

    ptx2=ptx*fx+u0
    pty2=pty*fy+v0

    pt=sympy.Matrix([ptx2,pty2])

    ptx,pty=sympy.symbols('v1_x, v1_y')
    pt[0]=((ptx-pt[0]))
    pt[1]=((pty-pt[1]))
    
    return pt

def getFInv():
    
    q=Quaternion.Quaternion()
    q1=sympy.symbols('v3_q1')
    q2=sympy.symbols('v3_q2')
    q3=sympy.symbols('v3_q3')
    q4=sympy.symbols('v3_q4')
    
    q.setValues(q1,q2,q3,q4)
    
    pt3=Quaternion.Quaternion()
    
    X=sympy.symbols('v2_X')
    Y=sympy.symbols('v2_Y')
    Z=sympy.symbols('v2_Z')
    d=sympy.symbols('v2_d')
    pt3.setValues(0,X/d,Y/d,Z/d)
    
    CX,CY,CZ=sympy.symbols('v3_CX,v3_CY,v3_CZ')
    pt3c=q.conj()*pt3*q
    pt3c[1]=pt3c[1]+CX
    pt3c[2]=pt3c[2]+CY
    pt3c[3]=pt3c[3]+CZ
    
    CX2,CY2,CZ2=sympy.symbols('v4_CX,v4_CY,v4_CZ')
    qq2=Quaternion.Quaternion()
    q12=sympy.symbols('v4_q1')
    q22=sympy.symbols('v4_q2')
    q32=sympy.symbols('v4_q3')
    q42=sympy.symbols('v4_q4')
    
    qq2.setValues(q12,q22,q32,q42)
    
    pt3c[1]=pt3c[1]-CX2
    pt3c[2]=pt3c[2]-CY2
    pt3c[3]=pt3c[3]-CZ2
    pt3c2=qq2*pt3c*qq2.conj()
    
    fx,fy,u0,v0,s,k1,k2,k3=sympy.symbols('v5_fx,v5_fy,v5_u0,v5_v0,v5_s,v5_k1,v5_k2,v5_k3')
    
    
    ptx=pt3c2[1]/pt3c2[3]
    pty=pt3c2[2]/pt3c2[3]
    
    ptx2=ptx*fx+u0
    pty2=pty*fy+v0

    pt=sympy.Matrix([ptx2,pty2])

    ptx,pty=sympy.symbols('v1_x, v1_y')
    pt[0]=((ptx-pt[0]))
    pt[1]=((pty-pt[1]))
    
    return pt

def getFQuat():
    q=Quaternion.Quaternion()
    q1=sympy.symbols('v1_q1')
    q2=sympy.symbols('v1_q2')
    q3=sympy.symbols('v1_q3')
    q4=sympy.symbols('v1_q4')
    

    CX,CY,CZ=sympy.symbols('v1_CX v1_CY v1_CZ')
    f=(1-(q1**2+q2**2+q3**2+q4**2))**2
       
    fm=sympy.Matrix([f])
    
    return fm


if __name__ == '__main__':
    print sympy.__version__
    if sympy.__version__=='0.7.2':
        print "Please use sympy version 0.7.1"
        print "sympy.cse returns an error."
        print "a bugfix in sympy is discuessed here:"
        print "http://www.mail-archive.com/sympy@googlegroups.com/msg15986.html"
    
        
    f1=getFInv()
    f2=getFProj()
    f3=getFQuat()
       
    names=['pt2','invDepth','cam','cam_in','pt3']
    variables=['x,y',['d'],'CX,CY,CZ,q1,q2,q3,q4','fx,fy,u0,v0,k1,k2,k3','X,Y,Z']
    addVariables=('','X,Y,Z','',['s'],'')
    addCode=('int id;','int id;\nint r,g,b;','','','int id;\nint r,g,b;')

    measFuncName=['ProjectMetricInv','ProjectMetric','Quat']
    measFuncInputF=[f1,f2,f3]
    
    measFuncInputObj=[[0,1,2,2,3],[0,4,2,3],[2]]
    measFuncInputObjOpt=[[1,2,3],[1,2],[0]]
    measFuncObjOpt=[1,2,4]
    
    genAll(names,\
       variables,\
       addVariables,\
       addCode,\
       measFuncName,\
       measFuncInputF,\
       measFuncInputObj,\
       measFuncInputObjOpt,\
       measFuncObjOpt,\
       filename_structs='structs',\
       min_function_name='minimize',\
       min_class_name='Minimizer',\
       filename_cpp='optimize',\
       libname_python='OptimizerGPUPy',\
       libname_cpp='OptimizerGPU',\
       clear_model_folder=True)
    
    

