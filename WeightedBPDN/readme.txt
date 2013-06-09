This package contains MATLAB code for solving iterative reweighted 
and adaptive reweighted L1 problems using the homotopy methods 
described in our paper:

"Fast and accurate algorithms for re-weighted L1 norm minimization," by 
M. Salman Asif and Justin Romberg
School of ECE, Georgia Tech.

For usage details, consult demo_wtBPDN.m and demo_adpWBPDN.m

This package also includes the following algorithms by
other authors (to allow running the comparative tests):

The SpaRSA algorithm, which can be downloaded from
http://www.lx.it.pt/~mtf/SpaRSA/

The YALL1 algorithm, which can be downloaded from
http://yall1.blogs.rice.edu/

The SPGL1 algorithm, which can be downloaded from
http://www.cs.ubc.ca/~mpf/spgl1/

(I have included copies of these solvers in this package as well)


This code is in development stage; any comments or bug reports
are very welcome.

Contact: sasif@gatech.edu         

Code website: http://users.ece.gatech.edu/~sasif/homotopy

%------------------------------------------------------------

To reproduce results in the paper, use the following scripts

job_wtBPDN_WAVE.m (Blocks and HeaviSine)

job_adpWBPDN.m (grayscale images)

** May need to compile some mex files for matrix-vector product
 noiselet and wavelet transforms ** 
    See compile.m in the main directory. 

The package also includes a modified version of SpaRSA (renamed as SpaRSA_adpW)
that adaptively changes weights at every continuation step. 
adp_wt = 1 sets the flag for adaptive reweighting

See SpaRSA_adpW.m for the details. 


%------------------------------------------------------------

Copyright (2012) M. Salman Asif 

This file is part of L1 homotopy toolbox.

The code distributed under the terms of the GNU General Public
License 3.0.

Permission to use, copy, modify, and distribute this software
for any purpose without fee is hereby granted, provided that
this entire notice is included in all copies of any software
which is or includes a copy or modification of this software
and in all copies of the supporting documentation for such
software. This software is being provided "as is", without any
express or implied warranty. In particular, the authors do not
make any representation or warranty of any kind concerning the
merchantability of this software or its fitness for any
particular purpose."

%------------------------------------------------------------
