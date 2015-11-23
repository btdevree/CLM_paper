# include <math.h>
# include <matrix.h>
# include <mex.h>

/*
 * Gauss_STORM_image_MEX Creates a STORM image from the list of (x,y) data; 
 *  mirrors the MATLAB function Gauss_STORM_image.cpp. 
 *
 * Creates a STORM image as a matrix of floating-point doubles. Assumes a
 *  Cartesian coordinate system. 
 *
 * Note: Type checking is done by the calling MATLAB function!!!
 *
 * Inputs
 *  resolution: floating-point double, number of nanometers per final image
 *      pixel
 *  xy_data: n by 2 array of floating-point doubles, center of each
 *      gaussian pdf
 *  covar_inv:  2 by 2 array of floating-point doubles, the inverse of the
 *      2D gaussian covariance matrix
 *  covar_inv: floating-point double, the determinant of the 2D gaussian
 *      covariance matrix
 *  calc_cutoff_pixels_x\y: 32 bit integer, extent in number of pixels on 
 *      either side of the center pixel to calculate the pdf distribution
 *      out to
 *  total_number_pixels_x\y: 32 bit integer, number of pixels in the x and
 *      y dimension of the final image.
 *  xmesh\ymesh: array of floating-point doubles with the coordinate of
 *      the center of each pixel, generated with meshgrid.
*/

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    plhs[0] = mxCreateDoubleMatrix (10, 10, mxREAL);
}