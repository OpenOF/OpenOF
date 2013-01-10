OpenOF: Framework for Sparse Non-linear Least Squares Optimization on a GPU
===================================================
Version 0.1

Hardware requirements:

- Nvidia GPU with compute capability >= 2.0

Software requirements:

The framework was tested using Ubuntu 12.04, 64 bit and libraries with the following version number.
- CUDA 5
- Python 2.7
- sympy 0.7.1  It is recommended to use this version as the sympy.cse function in version 0.7.2 does not work for some expression (http://www.mail-archive.com/sympy@googlegroups.com/msg15986.html).
- Thrust 1.6
- CUSP 0.3
- numpy 1.6.1

If you use OpenOF for you research please cite:

Cornelius Wefelscheid and Olaf Hellwich
OpenOF: Framework for sparse non-linear least squares optimization on a gpu.
VISAPP 2013

License:
   OpenOF - Open Optimization Framework
   Copyright (C) 2012 C. Wefelscheid


   OpenOF is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   OpenOF is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with OpenOF.  If not, see <http://www.gnu.org/licenses/>.

===================================================
Getting started:

OpenOF comes with two examples, similarity transformation and bundle adjustement. To get started we recommend locking in the simpler model of the similarity transformation. To get the similarity transformation example running you need to do the following steps.
- Move to OpenOF folder

$ cd Example
$ python similarityTransform.py
$ cd ..
$ cmake .
$ make

A C++ library as well as a Python library was created.

An example how to use the Python library can be found under ./Example
To run the Example proceed as follows:

$ cd Example

$ python runSimilarityTranformObj.py
or
$ python runSimilarityTranform.py


=================================================
Trouble Shooting:

In case you have version 0.7.2 of sympy you can downgrade with the following command:
>> sudo pip install sympy==0.7.1

Cuda 5 does not come with the Thrust 1.6, please exchange the thrust version usually under /usr/local/cuda/include/thrust



