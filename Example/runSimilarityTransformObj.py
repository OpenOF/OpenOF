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
sys.path.append('../')
import numpy

'''import the optimization library generated with OpenOF'''
import OptimizerSimTransPy as opt


def similarityTransformation(pt1,pt2):
    '''code example which shows how to use the optimization in an object oriented sense'''
    
    min=opt.MinimizerSimTrans()
    '''verbose=0 : no print statements; verbose=1 : initial and final values; verbose=2 : residual in each iteration'''
    min.verbose=1
    '''termination values of the levenberg marquardt algorithm'''

    min.eps1=1e-10
    min.eps2=1e-10
    '''maximum number of iterations'''
    min.kmax=100
    min.cg_it=150
    min.cg_thresh_abs=1e-12
    
    '''vector of the 3D points struct, to fill the values'''
    pt1V=min.pt3
    '''vector which contains all measurement combination'''
    meas=min.measurementCombinations
    '''create a vector which contains transformation structs with one element'''
    min.transform=opt.transformVector(1)
    transV=min.transform
    transV[0].CX=0.0
    transV[0].CY=0.0
    transV[0].CZ=0.0
    transV[0].q1=1.0
    transV[0].q2=0.0
    transV[0].q3=0.0
    transV[0].q4=0.0
    transV[0].s=1.0
    
    '''can be enabled to have a better initial estimate'''
    #delta=numpy.mean(pt1,axis=0)-numpy.mean(pt2,axis=0)
    #print delta
    #transV[0].CX=-delta[0]
    #transV[0].CY=-delta[1]
    #transV[0].CZ=-delta[2]
    
    '''insert the points in the struct and generate a measurement struct for each point pair'''
    for i,p1,p2 in zip(list(range(pt1.shape[0])),pt1,pt2):
        print(p1,p2)
        p1a=opt.pt3_t()
        p1a.X=p1[0]
        p1a.Y=p1[1]
        p1a.Z=p1[2]        
        pt1V.push_back(p1a)

        p2a=opt.pt3_t()
        p2a.X=p2[0]
        p2a.Y=p2[1]
        p2a.Z=p2[2]
        pt1V.push_back(p2a)
        
        m=opt.MeasurementCombinations_t()
        m.type=opt.eSimilarityTransformation
        #m.robust=opt.OF_SORT_WEIGHT
        #m.robust_para_a=1
        m.v.push_back(2*i)
        m.v.push_back(2*i+1)
        m.v.push_back(0)
        meas.push_back(m)
                
    '''inclue a measurement struct for the Quaternion'''
    m=opt.MeasurementCombinations_t()
    m.type=opt.eQuat
    '''the value corresponds to the index in the corresponding vector'''
    m.v.push_back(0)
    '''a inverse covariance can be set for each error/cost function'''
    #m.inv_cov.push_back(10);
    meas.push_back(m)
    min.verbose=0
    min.kmax=100
    min.run()
    
    res=min.residuals
    print("Residuals per function:")
    
    '''each measurement combination contains its final residual'''
    for m in meas:
        i=0
        for r in m.res:
            i+=1 
            print(i,":", r)
            
    
    print("Residual",res[res.size()-1]) 
    return numpy.array([transV[0].CX,transV[0].CY,transV[0].CZ,transV[0].q1,transV[0].q2,transV[0].q3,transV[0].q4,transV[0].s]),res[res.size()-1]

def quat2rot(arg1=None,arg2=None,arg3=None,arg4=None):
    '''
        quat2rot transforms a Quaternion to a rotation matrix
        as arguments its supports either a numpy array (arg1)
        or w,x,y,z = (arg1,arg2,arg3,arg4)
    '''
    if arg2 is None:
        assert arg1.shape==(4,)
        qw=arg1[0]
        qx=arg1[1]
        qy=arg1[2]
        qz=arg1[3]
    else:
        assert arg1 is not None
        assert arg2 is not None
        assert arg3 is not None
        assert arg4 is not None
        qw=arg1
        qx=arg2
        qy=arg3
        qz=arg4
    
    R=numpy.zeros((3,3))
    R[0,0]=1 - 2*qy*qy - 2*qz*qz
    R[0,1]=2*qx*qy - 2*qz*qw
    R[0,2]=2*qx*qz + 2*qy*qw
    R[1,0]=2*qx*qy + 2*qz*qw
    R[1,1]=1 - 2*qx*qx - 2*qz*qz    
    R[1,2]=2*qy*qz - 2*qx*qw
    R[2,0]=2*qx*qz - 2*qy*qw
    R[2,1]=2*qy*qz + 2*qx*qw
    R[2,2]=1 - 2*qx*qx - 2*qy*qy
    
    
    return R

def transform(trans,pt3):
    t=trans[0:3]
    qa=trans[3:7]
    R=quat2rot(qa)
    s=trans[7]
    
    pt3out=s*numpy.dot(R,pt3.T).T
    
    pt3out[:,0]+=t[0]
    pt3out[:,1]+=t[1]
    pt3out[:,2]+=t[2]
        
    return pt3out

if __name__ == '__main__':

    pt1=numpy.random.random((50,3))

    pt2=pt1.copy()
    '''apply noise to the original points'''
    pt1+=0.001*numpy.random.randn(50,3)
    R=quat2rot(0.70710678,  0. ,  0.70710678,  0.      )
    '''rotate points'''
    pt2=numpy.dot(R,pt2.T).T
    '''scale points'''
    pt2[:,0:3]*=10.0
    '''translate points'''
    pt2[:,0]+=10.0

    '''estimate the transformation which was applied'''
    trans,res=similarityTransformation(pt1,pt2)
    
    '''transform points pt1''' 
    pt1new=transform(trans,pt1)

    '''print estimated transformation'''
    print("translation:", trans[0:3])
    print("rotation quaternion:", trans[3:7])
    print("scale:", trans[7])
