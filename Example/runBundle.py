import sys
sys.path.append('../')
import numpy

import OptimizerGPUPy as opt
def savePly(filename,pt3):
	num_points=pt3.size()
	f=open(filename,'w+')
	f.write('ply\n')
	f.write('format ascii 1.0\n')
	f.write('element vertex %d\n'%num_points)
	f.write('property float x\n')
	f.write('property float y\n')
	f.write('property float z\n')
	f.write('property uchar red\n')
	f.write('property uchar green\n')
	f.write('property uchar blue\n')
	f.write('end_header\n')

	for p in pt3:
		f.write('%f %f %f %d %d %d\n'%(p.X,p.Y,p.Z,p.r,p.g,p.b))

	f.close()
	


def runBundle(filename):
    '''create a minimizer object'''
    min=opt.Minimizer()   

    f=open(filename)
    lines=f.readlines()
    

    unique_calib=lines[0]
    unique_calib=lines[0].find('FixedK')>0

    if unique_calib:
        print('Bundle with a unique calibration')
	print('the distortion parameters are currently not modelled')
        calib=opt.cam_in_t()
        calib.fx=float(lines[0].split(' ')[2])
        calib.fy=float(lines[0].split(' ')[4])
        '''the princple point is already substracted from the measurements in the nvm file '''
	calib.u0=0.0
        calib.v0=0.0
	calib.k1=0.0
#        calib.u0=float(lines[0].split(' ')[3])
#        calib.v0=float(lines[0].split(' ')[5])
        
        min.cam_in.push_back(calib)

    '''read number of cameras'''
    num_cams=int(lines[2])
    print('number of cams %d'%num_cams)
    
    '''read all cameras'''
    for l in lines[3:3+num_cams]:
        cam=opt.cam_t()
        ls=l.split(' ')

	'''fill the camera object with data'''
        cam.q1=float(ls[1])
        cam.q2=float(ls[2])
        cam.q3=float(ls[3])
        cam.q4=float(ls[4])
        
        cam.CX=float(ls[5])
        cam.CY=float(ls[6])
        cam.CZ=float(ls[7])
        
	'''insert the camera into the minimizer'''
        min.cam.push_back(cam)
	if not unique_calib:
		calib=opt.cam_in_t()
		calib.fx=float(ls[0].split('\t')[1])
		calib.fy=float(ls[0].split('\t')[1])
		calib.u0=0.0
		calib.v0=0.0
		calib.k1=float(ls[8])
		min.cam_in.push_back(calib)
	

	'''create a measurement combination which forces the quaternion of the camera to be of unit length'''
        m=opt.MeasurementCombinations_t()
	m.type=opt.eQuat
	m.v1=min.cam.size()-1
	'''set a high weight'''
	m.inv_cov.push_back(1000)
        min.measurementCombinations.push_back(m)

    
    num_points=int(lines[4+num_cams])
    print('number of 3D points %d'%num_points)

    '''read all 3D points'''
    for l in lines[5+num_cams:5+num_cams+num_points]:
        ls=l.split(' ')
        pt=opt.pt3_t()
        pt.X=float(ls[0])
        pt.Y=float(ls[1])
        pt.Z=float(ls[2])
        pt.r=int(ls[3])
        pt.g=int(ls[4])
        pt.b=int(ls[5])
        
        
        min.pt3.push_back(pt)

        num_obs=int(ls[6])
	'''read the measurements for each 3D point'''
        for i in range(num_obs):
            pt2=opt.pt2_t()
            cam_ind=int(ls[7+i*4])
            feature_ind=float(ls[8+i*4])
            pt2.x=float(ls[9+i*4])
            pt2.y=float(ls[10+i*4])
            
            min.pt2.push_back(pt2)
            
            '''create a measurement combination for each observation'''
            m=opt.MeasurementCombinations_t()
            m.type=opt.eProjectMetric
            m.v1=min.pt2.size()-1
            m.v2=min.pt3.size()-1
            m.v3=cam_ind
            '''calibration index is zero since only one calibration is used'''
            if unique_calib:
		m.v4=0
            else:		
            	m.v4=cam_ind
            
            min.measurementCombinations.push_back(m)
            
            
    f.close()
        
    min.verbose=2
    min.eps1=1e-6
    min.eps2=1e-6
    min.kmax=25
    min.run()
    
    '''compute pixel error'''
    pixel_err=0.0
    for i in range(num_cams,min.measurementCombinations.size()):
	err_x=min.measurementCombinations[i].res[0]
	err_y=min.measurementCombinations[i].res[1]
	pixel_err+=numpy.sqrt(err_x**2+err_y**2)

    pixel_err/=min.measurementCombinations.size()-num_cams
    print('Mean Pixel Error:%f'%pixel_err)
    print('save 3D points > points.ply')
    savePly('points.ply',min.pt3)
    print('finished')


if __name__=='__main__':
    print('read castle.nvm')
    runBundle('./data/castle.nvm')
    
    
    
