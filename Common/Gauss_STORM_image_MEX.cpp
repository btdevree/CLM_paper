// Include needed libraries
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

// Define the input and output names for convenience
# define RESOLUTION_IN  prhs[0]
# define XY_DATA_IN     prhs[1]
# define COVAR_INV_IN   prhs[2]
# define COVAR_DET_IN   prhs[3]
# define CUTOFF_PX_X_IN prhs[4]
# define CUTOFF_PX_Y_IN prhs[5]
# define XMESH_IN       prhs[6]
# define YMESH_IN       prhs[7]
# define IMAGE_OUT      plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    // Define variables
    int m, n, cutoff_x, cutoff_y;
    mxArray* image_ptr;
    double* image, xy_data, covar_inv, xmesh, ymesh;
    double value, resolution, covar_det;
    
    //Assign array variables from pointers
    xy_data = mxGetPr(XY_DATA_IN);
    covar_inv = mxGetPr(COVAR_INV_IN);
    xmesh = mxGetPr(XMESH_IN);
    ymesh = mxGetPr(YMESH_IN);
        
    // Get the dimensions of the image 
    m = mxGetM(XMESH_IN); 
    n = mxGetN(XMESH_IN);
        
    /*
    xmesh = mxGetPr(XMESH_IN); // Get pointer to xmesh
    ymesh = mxGetPr(YMESH_IN); // Get pointer to ymesh
    */
    mexPrintf("Hello world\n");

    image_ptr = mxCreateDoubleMatrix (m, n, mxREAL);
    image = mxGetPr(image_ptr);
    image[0] = 13;
    value = ymesh[0];
    
    mexPrintf("%g\n", value);
    
    IMAGE_OUT = image_ptr;
    
}




