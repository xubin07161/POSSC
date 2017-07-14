# POSSC
Polarimetric SAR Image Filtering based on Patch Ordering and Simultaneous Sparse Coding

%  POSSC is an algorithm for POLSAR image filtering. 
%  This algorithm reproduces the results from the article:
%  [1] B. Xu et al. 'Polarimetric SAR Image Filtering based on Patch
%      Ordering and Simultaneous Sparse Coding'
%  Please refer to this paper for a more detailed description of the algorithm.
%
%  BASIC USAGE EXAMPLES:
%
%     1) Using the default parameters
% 
%      img_filtered = POLSAR_POSSC(img,'ENL',ENL)
% 
%  INPUT ARGUMENTS (OPTIONAL):
%
%     1) img : The input POLSAR image. Each pixel should be the vector form
%     of covariance matrix.
%
%     2) ENL : The equivalent number of looks. The ENL can be obtained by 
%              by supervised or unsupervised estimation. For a homogeneous 
%              region, the ENL can be calculated by 
%              ENL=Tr(<C>)^2/(<Tr(CC)>-Tr(<C><C>)).
%
%  OUTPUTS:
%     1) img_filtered  : The filtered image. Each pixel is the vector form
%     of covariance matrix.                                              

-------------------------------------------------------------------
 Copyright
-------------------------------------------------------------------

Copyright (c) 2014 Bin Xu
All rights reserved.
This work should only be used for nonprofit purposes.

-------------------------------------------------------------------
 Requirements
-------------------------------------------------------------------

All the functions and scripts were tested on MATLAB 2011b,
the operation is not guaranteed with older version of MATLAB.

The version requires Windows 64 bit operating system.

-------------------------------------------------------------------
 Execution times
-------------------------------------------------------------------
Intel(R) Core(TM) i5-2300 2.80Ghz 64bit; 8Gb Memory

dim.    |      time      |  
400x600 |  around 60 sec |  

-------------------------------------------------------------------
 Feedback
-------------------------------------------------------------------

If you have any comment, suggestion, or question, please do
contact Bin Xu at xubin07161@gmail.com
