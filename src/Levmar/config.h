/*
 *   OpenOF - Open Optimization Framework
 *   Copyright (C) 2012 C. Wefelscheid
 *
 *   This file is part of OpenOF.
 *

 */



#ifndef CONFIG_H_
#define CONFIG_H_


//Define if either double or float should be used
typedef double oof_float;

/*If enabled the host memory is uses for the measurement structs, for everything else the memory which is assigned below is used */
//#define useHostMemory

//for small least squares adjustment the cpu can be used as well (EXPERIMENTAL) (comment line below)
#define USE_GPU

#endif /* FUNCTORS_H_ */
