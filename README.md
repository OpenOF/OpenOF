OpenOF: Framework for Sparse Non-linear Least Squares Optimization on a GPU
===================================================
Version 0.2

License update to Apache License to allow commercial use.

Hardware requirements:
------------------------------------------

- Nvidia GPU with compute capability >= 2.0

Software requirements:
------------------------------------------

The framework was tested using Ubuntu 18.04, 64 bit and libraries with the following version number.
- CUDA 9.1
- Python 3.5
- sympy 1.3
- Thrust 1.6
- CUSP 0.5.1
- numpy 1.14.2
- swig 3.0.8

Cite:
------------------------------------------
If you use OpenOF for your research please cite:

Cornelius Wefelscheid and Olaf Hellwich
OpenOF: Framework for sparse non-linear least squares optimization on a gpu.
VISAPP 2013 [[pdf]](http://www.cv.tu-berlin.de/fileadmin/fg140/OpenOF.pdf)

License:
------------------------------------------
   OpenOF - Open Optimization Framework
   Copyright (C) 2012 C. Wefelscheid


   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

      http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

Getting started:
------------------------------------------

OpenOF comes with two examples, similarity transformation and bundle adjustement. To get started we recommend looking into the simpler model of the similarity transformation. To get the similarity transformation example running you need to execute the following commands.

----------------------
> $ git clone https://github.com/OpenOF/OpenOF.git
>
> $ cd OpenOF
>
> $ git submodule init
>
> $ git submodule update
>
> $ cd Example
> 
> $ python3 similarityTransform.py
> 
> $ cd ..
>
> $ cmake .
>
> $ make

---------------------

A C++ library as well as a Python library was created.

An example how to use the Python library can be found under ./Example.
To run the Example proceed as follows:

------------
> $ cd Example
> 
> $ python3 runSimilarityTransformObj.py

------------

or

------------
> $ python3 runSimilarityTransform.py

------------

A more complex example is under Example/bundle.py (design) and Example/runBundle.py (runtime).

Docker:
------------------------------------------
Use the provided docker from inside vscode.
Install the remote development extension within vscode.

Inside the docker you can continue from getting startet step:
"cd Example"

Trouble Shooting:
------------------------------------------

- In case you have version 0.7.2 of sympy you can downgrade with the following command

------------
>
> $ sudo pip install sympy==0.7.1

------------
or upgrade to version 0.7.3.
- Cuda 5 does not come with Thrust 1.6, please exchange the thrust version usually under /usr/local/cuda/include/thrust.
- Cuda 5 is not working with gcc version 4.7 or higher. Try using an older gcc version.
- Compiler error: /bin/sh: 1: Syntax error: "(" unexpected.
      
      Make sure the path of the library doesn't contain any brackets.





