/* Gauss_STORM_image_MEX Creates a STORM image from the list of (x,y) data; 
 *  mirrors the MATLAB function Gauss_STORM_image.cpp. 
 *
 * Creates a STORM image as a matrix of floating-point doubles. Assumes a
 *  Cartesian coordinate system, but the image does not need to contain the 
 *  origin.
 *
 * Note: Type checking is done by the calling MATLAB function!!! The C++ 
 *  function can handle xy data points outside of the image boundries, but
 *  most other errors are not checked. Don't give it bad values.
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
 *  x_vector\y_vector: Column vector arrays of floating-point doubles with 
 *      the coordinate of the center of each pixel in the x or y direction.
 *      Values are in increasing order (i.e., starting from the bottom-left 
 *      corner of the Cartesian plane of the image). 
 * Outputs
 *  image: MATLAB array of floating-point double values. 
*/

// Include needed libraries
# include <math.h>
# include <matrix.h>
# include <mex.h>
# include <algorithm> // min_element, max_element

// Get the math constants from the math.h header
# define _USE_MATH_DEFINES // to get pi

// Define the input and output names for convenience
# define XY_DATA_IN     prhs[0]
# define RESOLUTION_IN  prhs[1]
# define COVAR_INV_IN   prhs[2]
# define COVAR_DET_IN   prhs[3]
# define CUTOFF_PX_X_IN prhs[4]
# define CUTOFF_PX_Y_IN prhs[5]
# define X_VECTOR_IN    prhs[6]
# define Y_VECTOR_IN    prhs[7]
# define IMAGE_OUT      plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    // Define variables
    int cutoff_x, cutoff_y, x_ctrpx, y_ctrpx, x_minpx, y_minpx, x_maxpx, y_maxpx;
    mxArray *image_ptr;
    double *image, *xy_data, *covar_inv, *x_vector, *y_vector;
    double resolution, covar_det;
    mwSize image_m, image_n, xy_data_m; // Equlivent to int or size_t, depending on system and -largeArrayDims compiler flag
    double x, y, half_res, left_edge, right_edge, bottom_edge, top_edge;
    double delta_x, delta_y, gauss_scale_factor, exponent, value;
    mwSize xy_data_index; // We use the index values outside of the loop definition, so we want to delcare them along with all other variables (-fpermissive issue)
    int x_index = -10, y_index = -10;
    
    //Assign array variables from pointers
    xy_data = mxGetPr(XY_DATA_IN);
    covar_inv = mxGetPr(COVAR_INV_IN);
    x_vector = mxGetPr(X_VECTOR_IN);
    y_vector = mxGetPr(Y_VECTOR_IN);
        
    // Get the dimensions of the image and data
    image_n = mxGetM(X_VECTOR_IN); //Both x and y vectors given as column vectors
    image_m = mxGetM(Y_VECTOR_IN);
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
    left_edge = *(std::min_element(x_vector, x_vector + image_n)) - half_res;
    right_edge = *(std::max_element(x_vector, x_vector + image_n)) + half_res;
    bottom_edge = *(std::min_element(y_vector, y_vector + image_m)) - half_res;
    top_edge = *(std::max_element(y_vector, y_vector + image_m)) + half_res;
    
    // Precalculate the scaling factor for the 2D Gausian
    gauss_scale_factor = pow(resolution, 2) / (2 * M_PI * sqrt(covar_det));
        
    // Loop through each datapoint
    for(xy_data_index = 0; xy_data_index < xy_data_m; xy_data_index++) 
    {    
        // Copy the x and y value
        x = xy_data[xy_data_index];
        y = xy_data[xy_data_index + xy_data_m];
        
        // Calculate the center pixel for the xy_datapoint, uses 0-based 2D indexing from top left corner
        x_ctrpx = static_cast<int>(floor((x - left_edge) / resolution));
        y_ctrpx = static_cast<int>(image_m) - static_cast<int>(ceil((y - bottom_edge) / resolution));
        
        // Check if the center pixel is out of bounds
        if(x_ctrpx < 0 || x_ctrpx > image_n || y_ctrpx < 0 || y_ctrpx > image_m)
        {
            continue; // Skip datapoint if it lies outside of image bounds
        }
        
       // Calculate the pixel indices for the box of values we are going to compute
        x_minpx = x_ctrpx - cutoff_x;
        x_minpx = (x_minpx < 0) ? 0 : x_minpx; // Ternary operator works like: x = (condition) ? (value_if_true) : (value_if_false);
        x_maxpx = x_ctrpx + cutoff_x;
        x_maxpx = (x_maxpx >= image_n) ? image_n - 1 : x_maxpx;
        y_minpx = y_ctrpx - cutoff_y;
        y_minpx = (y_minpx < 0) ? 0 : y_minpx; 
        y_maxpx = y_ctrpx + cutoff_y;
        y_maxpx = (y_maxpx >= image_m) ? image_m - 1 : y_maxpx;
        
        // Loop through each pixel in the box 
        for(x_index = x_minpx ; x_index <= x_maxpx; x_index++)
        {
            for(y_index = y_minpx; y_index <= y_maxpx; y_index++)
            {
                // Calculate value of 2D Gaussian at the pixel
                delta_x = x_vector[x_index] - x;
                delta_y = top_edge - y_vector[y_index] - y; // y_vector runs opposite direction of the pixel indices
                exponent = -0.5 * (delta_x * delta_x * covar_inv[0] + delta_x * delta_y * (covar_inv[1] + covar_inv[2]) 
                                                                                    + delta_y * delta_y * covar_inv[3]);
                value = gauss_scale_factor * exp(exponent);
                
                // Add value to the image matrix
                image[x_index * image_m + y_index] = image[x_index * image_m + y_index] + value;
            }
        } // Close pixel loops    
    } //Close datapoint loop
    
    // Set the output pointer to the image
    IMAGE_OUT = image_ptr;
}
