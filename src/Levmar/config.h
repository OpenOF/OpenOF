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



#ifndef CONFIG_H_
#define CONFIG_H_


//Define if either double or fload should be used
typedef double oof_float;

/*If enabled the host memory is uses for the measurement structs, for everything else the memory which is assigned below is used */
//#define useHostMemory

//for small least squares adjustment the cpu can be used as well (EXPERIMENTAL) (comment line below)
#define USE_GPU

#endif /* FUNCTORS_H_ */
