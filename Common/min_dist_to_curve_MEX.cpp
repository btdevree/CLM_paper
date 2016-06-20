/* MIN_DIST_TO_CURVE Calculates the minimum distance between each of the x
 * and y coordinates to any point on the sampled curves.
 *
 * Inputs:
 *   coords: floating-point double n by 2 matrix of x and y coordinates to
 *       measure the distances from.
 *   curve: floating-point double m by 2 matrix of x and y coordinates that
 *       define a curve to measure the coords against.
 * Output:
 *   distances: floating-point double n by 2 matrix of Eucledian distances 
 *       from the points given in coords to any point contained in curve.\
*/

// Include needed libraries
# include <matrix.h>
# include <mex.h>
# include <limits>
# include <math.h>

// Define the input and output names for convenience
# define COORDS_IN      prhs[0]
# define CURVE_IN       prhs[1]
# define DISTANCES_OUT  plhs[0]

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    // Define variables
    mwSize coord_length, curve_length; // Equlivent to int or size_t, depending on system and -largeArrayDims compiler flag
    mxArray *distances_ptr;
    double *coords, *curve, *distances;
    double curve_x, curve_y, coord_x, coord_y, new_distance;
    double inf = std::numeric_limits<double>::infinity();
    mwSize curve_index, coord_index; // We use the index values outside of the loop definition, so we want to delcare them along with all other variables (-fpermissive issue)
    
    //Assign array variables from pointers
    coords = mxGetPr(COORDS_IN);
    curve = mxGetPr(CURVE_IN);
    
    // Get lengths of matrices
    coord_length = mxGetM(COORDS_IN);
    curve_length = mxGetM(CURVE_IN);
    
    // Create output matrix
    distances_ptr = mxCreateDoubleMatrix(coord_length, 1, mxREAL);
    distances = mxGetPr(distances_ptr);
    
    // Fill output with "infinite" values so the first iteration will fill the matrix with all the measured distances 
    for(coord_index = 0; coord_index < coord_length; coord_index++)
    {
        distances[coord_index] = inf;
    }
    
    // Calculate the distances from all coord points to the curve points
    for(curve_index = 0; curve_index < curve_length; curve_index++) 
    {    
        // Get the current curve point
        curve_x = curve[curve_index];
        curve_y = curve[curve_length + curve_index];
        
        // Get the distance between the curve point and all the coord points
        for(coord_index = 0; coord_index < coord_length; coord_index++) 
        {  
            // Get the current coord point
            coord_x = coords[coord_index];
            coord_y = coords[coord_length + coord_index];

            // Find distance from the current curve point to the current coord point
            new_distance = sqrt(pow(curve_x - coord_x, 2) + pow(curve_y - coord_y, 2));
            
            // If the new distance is smaller than the old one, replace the old one
            if(new_distance < distances[coord_index])
            {
                distances[coord_index] = new_distance;
            } 
        } // Close coord loop    
    } //Close curve loop

    // Set the output pointer
    DISTANCES_OUT = distances_ptr;
}
