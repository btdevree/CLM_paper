/* sample_inverse_cdf_MEX Convert uniform [0-1] random numbers into samples 
 * from a cdf function. Mirrors the MATLAB function Gauss_STORM_image.cpp
 *
 *  Calculates the indices as the highest index value where the random 
 *      number is greater than or equal to the cdf value at that index.  
 *
 * Note: Type checking is done by the calling MATLAB function!!! The C++ 
 *  function can't handle numbers outside of the expected range. Don't give 
 *  it bad values.
 *
 * Inputs:
 *  sample_rands: column vector of floating-point doubles, array of random 
 *      numbers between 0 and 1.
 *  cdf: column_vector of floating-point doubles, monotonically
 *      non-decreasing values between 0 and 1.
 * Output:
 *  indices: mwSize MATLAB array (usually 32 bit integer), the index values
 *      of the cdf entry that corrospond to the values of the random 
 *      samples.
*/

// Include needed libraries
# include <math.h>  // round
# include <matrix.h>
# include <mex.h>

// Define the input and output names for convenience
# define SAMPLE_RANDS_IN prhs[0]
# define CDF_IN          prhs[1]
# define INDICES_OUT     plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    // Define variables
    mxArray *indices_ptr;
    double *sample_rands, *cdf;
    double rand_value, test_value;
    mwSize *indices; // mwSize is equlivent to int or size_t, depending on system and -largeArrayDims compiler flag
    mwSize rands_m, cdf_m, rand_index, lower_lim, upper_lim, test_index; 
            
    //Assign array variables from pointers
    sample_rands = mxGetPr(SAMPLE_RANDS_IN);
    cdf = mxGetPr(CDF_IN);
    
    // Get the length of the sample_rands and cdf vectors
    rands_m = mxGetM(SAMPLE_RANDS_IN); 
    cdf_m = mxGetM(CDF_IN);
    
    // Make the output array
    indices_ptr = mxCreateNumericMatrix(rands_m, 1, mxINT32_CLASS, mxREAL);
    indices = static_cast<mwSize*>(mxGetData(indices_ptr));
    
    // Loop through each sample_rand value
    for(rand_index = 0; rand_index < rands_m; rand_index++) 
    {    
        // Copy the sample_rands value
        rand_value = sample_rands[rand_index];
           
        // Reset the upper and lower index limits
        upper_lim = cdf_m - 1;
        lower_lim = 0; 
        
        // Search for the corrosponding cdf value by repeated comparison at halfway between the min and max limits
        while(upper_lim - lower_lim > 1)
        {
            // Average the upper and lower limits and get the corrosponding cdf value
            test_index = static_cast<mwSize>(round((upper_lim + lower_lim)/2));
            test_value = cdf[test_index];

            // If the random value is smaller than the test value, move the upper limit to the test index
            if(rand_value < test_value)
            {
                upper_lim = test_index;
            }
            // If the random value is equal to or larger than the test value, move the lower limit to the test index
            else
            {
                lower_lim = test_index;
            }
        } //end while
        
        // Write the lower limit into the output matrix
        indices[rand_index] = lower_lim + 1; // Convert indices to MATLAB (first index = 1)
    } //end for
        
    // Set the output pointer to the image
    INDICES_OUT = indices_ptr;
}
