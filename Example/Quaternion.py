'''

 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *

 
'''


import sympy
import sympy.core.numbers
import numpy
import math

class Quaternion(object):
    '''
    classdocs
    '''


    def __init__(self):
        '''
        Constructor
        '''
        w,x,y,z=sympy.symbols('w x y z')
        self.q=sympy.Matrix([w,x,y,z])
    def __getitem__(self,offset):
        return self.q[offset]
    def __setitem__(self,offset,value):
        self.q[offset]=value
    def __add__(self,right):
        q=Quaternion()
        q.setValues(self[0]+right[0], self[1]+right[1], self[2]+right[2],self[3]+right[3])
        return q
    def __sub__(self,right):
        q=Quaternion()
        q.setValues(self[0]-right[0], self[1]-right[1], self[2]-right[2],self[3]-right[3])
        return q
    def __mul__(self,right):
        if isinstance(right, Quaternion):
            return self.mul(right)
        else:
            qout=Quaternion()
            qout.setValues(self[0]*right, self[1]*right, self[2]*right, self[3]*right)
            return qout
    __rmul__ = __mul__        
            
    def __str__(self):
        return self.q.__str__()
    
    def setValues(self,q,x=None,y=None,z=None):
        if x is None:
            self.q=q
        else:
            self.q[0]=q
            self.q[1]=x
            self.q[2]=y
            self.q[3]=z
        
    def toRotationSympy(self):
        
        w=self.q[0]
        x=self.q[1]
        y=self.q[2]
        z=self.q[3]
        rot1temp=sympy.Matrix([[w,z,-y,x], [-z, w, x, y],[y, -x, w, z],[ -x, -y, -z, w]])
        rot2temp=sympy.Matrix([[w, z, -y, -x], [-z, w, x, -y],[y, -x, w, -z],[x, y, z, w]])
        rot=rot1temp*rot2temp
        rot=rot[0:3,0:3].transpose()
        return rot

    def toRotationSympy2(self):
        qw=self.q[0]
        qx=self.q[1]
        qy=self.q[2]
        qz=self.q[3]
        #rot1temp=sympy.Matrix([[w,z,-y,x], [-z, w, x, y],[y, -x, w, z],[ -x, -y, -z, w]])
        R=sympy.Matrix()
        R=R.zeros((3,3))
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
    def toRotationNumpy(self):
        w=self.q[0]
        x=self.q[1]
        y=self.q[2]
        z=self.q[3]
        
        #assert isinstance(w,sympy.core.numbers)
        #assert isinstance(x,sympy.core.numbers)
        #assert isinstance(y,sympy.core.numbers)
        #assert isinstance(z,sympy.core.numbers)
        
        rot1temp=sympy.Matrix([[w,z,-y,x], [-z, w, x, y],[y, -x, w, z],[ -x, -y, -z, w]])
        
        rot2temp=sympy.Matrix([[w, z, -y, -x], [-z, w, x, -y],[y, -x, w, -z],[x, y, z, w]])
        rot=rot1temp*rot2temp
       # R=sympy.Matrix.zeros((3,3))
       # R[0,0]=1 - 2*qy*qy - 2*qz*qz
       # R[0,1]=2*qx*qy - 2*qz*qw
       # R[0,2]=2*qx*qz + 2*qy*qw
       # R[1,0]=2*qx*qy + 2*qz*qw
       # R[1,1]=1 - 2*qx*qx - 2*qz*qz    
       # R[1,2]=2*qy*qz - 2*qx*qw
       # R[2,0]=2*qx*qz - 2*qy*qw
       # R[2,1]=2*qy*qz + 2*qx*qw
       # R[2,2]=1 - 2*qx*qx - 2*qy*qy
    
        
        rot=rot[0:3,0:3].transpose()
        return numpy.array(rot)
    def conj(self):
        qconj=Quaternion()
        qconj.setValues(self[0], -self[1], -self[2], -self[3])
        return qconj
    def mul(self,q2):
        a=self[0];
        b=self[1];
        c=self[2];
        d=self[3];
    
        e=q2[0];
        f=q2[1];
        g=q2[2];
        h=q2[3];
        qoutw=a*e - b*f - c*g- d*h
        qoutx=(b*e + a*f + c*h - d*g)
        qouty=(a*g - b*h + c*e + d*f)
        qoutz=(a*h + b*g - c*f + d*e)
        qout=Quaternion()
        qout.setValues(qoutw,qoutx,qouty,qoutz)
        return qout
    
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
def rot2quat(m):
    tr = m[0,0] + m[1,1] + m[2,2]

    if (tr > 0): 
        S = math.sqrt(tr+1.0) * 2 
        qw = 0.25 * S
        qx = (m[2,1] - m[1,2]) / S
        qy = (m[0,2] - m[2,0]) / S 
        qz = (m[1,0] - m[0,1]) / S 
    elif ((m[0,0] > m[1,1]) and (m[0,0] > m[2,2])): 
        S = math.sqrt(1.0 + m[0,0] - m[1,1] - m[2,2]) * 2 
        qw = (m[2,1] - m[1,2]) / S
        qx = 0.25 * S
        qy = (m[0,1] + m[1,0]) / S 
        qz = (m[0,2] + m[2,0]) / S
    elif m[1,1] > m[2,2]: 
        S = math.sqrt(1.0 + m[1,1] - m[0,0] - m[2,2]) * 2
        qw = (m[0,2] - m[2,0]) / S
        qx = (m[0,1] + m[1,0]) / S 
        qy = 0.25 * S
        qz = (m[1,2 ]+ m[2,1]) / S
    else: 
        S = math.sqrt(1.0 + m[2,2] - m[0,0] - m[1,1]) * 2
        qw = (m[1,0] - m[0,1]) / S;
        qx = (m[0,2] + m[2,0]) / S;
        qy = (m[1,2] + m[2,1]) / S;
        qz = 0.25 * S;
        
    return numpy.array([qw,qx,qy,qz])
def quatMul(q1,q2):
    a=q1[0];
    b=q1[1];
    c=q1[2];
    d=q1[3];

    e=q2[0];
    f=q2[1];
    g=q2[2];
    h=q2[3];
    qoutw=a*e - b*f - c*g- d*h
    qoutx=(b*e + a*f + c*h - d*g)
    qouty=(a*g - b*h + c*e + d*f)
    qoutz=(a*h + b*g - c*f + d*e)
    return numpy.array([qoutw,qoutx,qouty,qoutz])
    #return qoutw, qoutx, qouty, qoutz

def angle2quat(ax,ay,az,angle):
    
    qw = numpy.cos(angle/2)
    qx = ax * numpy.sin(angle/2)
    qy = ay * numpy.sin(angle/2)
    qz = az * numpy.sin(angle/2)
    return numpy.array([qw,qx,qy,qz])
    #return qw, qx, qy, qz

def conjugate(q):
    qout=q.copy()
    qout[1:4]=-qout[1:4]
    return qout
