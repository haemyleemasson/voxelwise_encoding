
# Voxel-wise encoding model 
With this matlab code you will be able to run encoding model analysis with the banded ridge regression method proposed by ([Nunez-Elizalde et al., 2019](https://www.sciencedirect.com/science/article/abs/pii/S1053811919302988?via%3Dihub)). 
Instead of the ridge function (a built-in matlab function), bridge function is used to have more than 1 ridge penalty. In this function, three methods are available, but they give (almost) identical results. 
 1. Standard Tikhonov Regression method 
 2. Standard Form Solution 
 3. Singular Value Decomposition (the fastest solution!)
 
bridge function provide beta estimates and corresponding lambda values. 
 # Dependencies
[MATLAB](https://www.mathworks.com/products/matlab.html), [COSMOMVPA toolbox](http://www.cosmomvpa.org/), [NIFTI toolbox](https://www.nitrc.org/projects/nifti/)

