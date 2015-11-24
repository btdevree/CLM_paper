// Include needed libraries
# include <math.h>
# include <matrix.h>
# include <mex.h>
# include <algorithm> // min_element, max_element

# define _USE_MATH_DEFINES // to get pi from math.h

/*
 * Gauss_STORM_image_MEX Creates a STORM image from the list of (x,y) data; 
 *  mirrors the MATLAB function Gauss_STORM_image.cpp. 
 *
 * Creates a STORM image as a matrix of floating-point doubles. Assumes a
 *  Cartesian coordinate system, but the image does not need to contain the 
 *  origin.
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
# define XY_DATA_IN     prhs[0]
# define RESOLUTION_IN  prhs[1]
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
    int cutoff_x, cutoff_y, x_ctrpx, y_ctrpx, x_minpx, y_minpx, x_maxpx, y_maxpx;
    mxArray *image_ptr;
    double *image, *xy_data, *covar_inv, *xmesh, *ymesh;
    double resolution, covar_det;
    mwSize image_m, image_n, xy_data_m; // Equlivent to int or size_t, depending on system and -largeArrayDims compiler flag
    double x, y, half_res, left_edge, right_edge, bottom_edge, top_edge;
    double delta_x, delta_y, gauss_scale_factor, value;
    
    
    //Assign array variables from pointers
    xy_data = mxGetPr(XY_DATA_IN);
    covar_inv = mxGetPr(COVAR_INV_IN);
    xmesh = mxGetPr(XMESH_IN);
    ymesh = mxGetPr(YMESH_IN);
        
    // Get the dimensions of the image and data
    image_m = mxGetM(XMESH_IN); 
    image_n = mxGetN(XMESH_IN);
    xy_data_m = mxGetM(XY_DATA_IN);
    
    // Get the cutoff distances and resolution
    cutoff_x = static_cast<int>(mxGetScalar(CUTOFF_PX_X_IN));
    cutoff_y = static_cast<int>(mxGetScalar(CUTOFF_PX_Y_IN));
    resolution = mxGetScalar(RESOLUTION_IN);
    covar_det = mxGetScalar(COVAR_DET_IN);
    
    // Make the output image
    image_ptr = mxCreateDoubleMatrix(image_m, image_n, mxREAL);
    image = mxGetPr(image_ptr);
    
    // Calculate Cartesian boundries of the image
    half_res = resolution * 0.5;
    left_edge = *(std::min_element(xmesh, xmesh + image_m*image_n) - half_res);
    right_edge = *(std::max_element(xmesh, xmesh + image_m*image_n) + half_res);
    bottom_edge = *(std::min_element(ymesh, ymesh + image_m*image_n) - half_res);
    top_edge = *(std::max_element(ymesh, ymesh + image_m*image_n) + half_res);
    
    // Precalculate the scaling factor for the 2D Gausian
    gauss_scale_factor = (resolution * resolution) / (2 * M_PI * sqrt(covar_det));
    
    // Print stuff
    mexPrintf("image_m = %u, image_n = %u, xy_data_m = %u\n", image_m, image_n, xy_data_m);
    mexPrintf("x min = %.2f, x max = %.2f, y min = %.2f, y max = %.2f\n", left_edge, right_edge, bottom_edge, top_edge);
    
    // Loop through each datapoint
    for(mwSize xy_data_index = 0; xy_data_index < xy_data_m; xy_data_index++) 
    {    
        // Copy the x and y value
        x = xy_data[xy_data_index];
        y = xy_data[xy_data_index + xy_data_m];
        
        mexPrintf("index = %u, x = %.2f, y = %.2f\n", xy_data_index, x, y);

        // Calculate the center pixel for the xy_datapoint, uses 0-based 2D indexing from top left corner
        x_ctrpx = static_cast<int>floor((x - left_edge) / resolution);
        y_ctrpx = static_cast<int>(image_m - ceil((y - bottom_edge) / resolution));
        
        // Check if the center pixel is out of bounds
        if(x_ctrpx < 0 || x_ctrpx >= image_n || y_ctrpx < 0 || y_ctrpx >= image_m)
        {
            continue // Skip datapoint if it lies outside of image bounds
        }
        
        // Calculate the pixel indices for the box of values we are going to compute
        x_minpx = x_ctrpx - x_cutoff;
        x_minpx = (x_minpx < 0) ? 0 : x_minpx; // Ternary operator works like: x = (condition) ? (value_if_true) : (value_if_false);
        x_maxpx = x_ctrpx + x_cutoff;
        x_maxpx = (x_maxpx >= image_n) ? image_n - 1 : x_maxpx;
        y_minpx = y_ctrpx - y_cutoff;
        y_minpx = (y_minpx < 0) ? 0 : y_minpx; 
        y_maxpx = y_ctrpx + y_cutoff;
        y_maxpx = (y_maxpx >= image_m) ? image_m - 1 : y_maxpx;
        
        // Loop through each pixel in the box
        for(int x_index = x_minpx; x_index <= x_maxpx; x_index++);
        {
            for(int y_index = y_minpx; y_index <= y_maxpx; y_index++);
            {
                 // Calculate value of 2D Gaussian at the pixel
                delta_x = xmesh[x_index * image_m + y_index] - x;
                delta_y = ymesh[x_index * image_m + y_index] - y;
                exponent = -0.5 * (delta_x * delta_x * covar_inv[0] 
                        + delta_x * delta_y * (covar_inv[1] + [covar_inv[2]) + delta_y * delta_y * covar_inv[3]);
                value = gauss_scale_factor * exp(exponent);
                
                // Add value to the image matrix
                image[x_index * image_m + y_index] = image[x_index * image_m + y_index] + value;
            }
        }
    }
   
    
    // Testing junk
    mexPrintf("Hello world\n");

    
    image[4] = 13;
    value = cutoff_x;
        
    IMAGE_OUT = image_ptr;
    
}




