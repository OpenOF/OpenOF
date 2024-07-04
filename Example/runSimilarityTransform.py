'''

 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *


'''

import sys
sys.path.append('../')
import numpy

#from IO.undist import undist

import OptimizerSimTransPy as opt
def similarityTransformation(pt1,pt2):
    
    '''see runSimilarityTransformObj.py for detailed documentation'''
    pt1V=opt.pt3Vector()
    
    meas=opt.MeasurementCombinationsVector(0)
    transV=opt.transformVector(1)
    transV[0].CX=0.0
    transV[0].CY=0.0
    transV[0].CZ=0.0
    transV[0].q1=1.0
    transV[0].q2=0.0
    transV[0].q3=0.0
    transV[0].q4=0.0
    transV[0].s=1.0
    delta=numpy.mean(pt1,axis=0)-numpy.mean(pt2,axis=0)
    #print delta
   # transV[0].CX=-delta[0]
   # transV[0].CY=-delta[1]
   # transV[0].CZ=-delta[2]
    
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
        m.v.push_back(2*i)
        m.v.push_back(2*i+1)
        m.v.push_back(0)
        meas.push_back(m)
                
    res=opt.doubleVector()
    m=opt.MeasurementCombinations_t()
    m.type=opt.eQuat
    m.v.push_back(0)
    meas.push_back(m)
    opt.minimize_simtrans(pt1V,transV,meas,1000,res)
    
        
    print("Residual",res[res.size()-1]) 
    return numpy.array([transV[0].CX,transV[0].CY,transV[0].CZ,transV[0].q1,transV[0].q2,transV[0].q3,transV[0].q4,transV[0].s]),res[res.size()-1]

def quat2rot(arg1=None,arg2=None,arg3=None,arg4=None):
    '''
        quat2dcm transforms a Quaternion to a rotation matrix
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
    C=trans[0:3]
    qa=trans[3:7]
    R=quat2rot(qa)
    s=trans[7]
    
    pt3out=s*numpy.dot(R,pt3.T).T
    
    pt3out[:,0]+=C[0]
    pt3out[:,1]+=C[1]
    pt3out[:,2]+=C[2]
        
    return pt3out

if __name__ == '__main__':


    pt1=numpy.random.random((50,3))

    pt2=pt1.copy()
    '''apply noise to the original points'''
    pt1+=0.0*numpy.random.randn(50,3)
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
