function [expected_radav_xcor] = alaising_test_analytical_solution(distance_vector, params, passed_vars)
%ALAISING_TEST_ANALYTICAL_SOLUTION Numerical calculation of the expected 
%   radially averaged crosscorrelation computed with analytical functions.

% Get required parameters
ch2_seperation_length = params.ch2_crosscor_params(1);
analytical_calc_cutoff = passed_vars.analytical_calc_cutoff; % sigma
analytical_arc_length_sampling = passed_vars.analytical_arc_length_sampling;
pdf_sigma = params.STORM_precision;

% Calculate analytical solution (work in nanaometers)
% Calc average spot density
ch1_total_density = params.number_events_ch1 / (pi * passed_vars.cell_radius^2); % spots/nm^2
ch2_total_density = params.number_events_ch2 / (pi * (passed_vars.cell_radius + params.ch2_crosscor_params(1))^2); % spots/nm^2
ch1_rand_density = (params.number_events_ch1 - 1) / (pi * passed_vars.cell_radius^2); % spots/nm^2 Density of all other spots in the channel besides the one modeled in the psf
ch2_rand_density = (params.number_events_ch2 - 1) / (pi * (passed_vars.cell_radius + params.ch2_crosscor_params(1))^2); % spots/nm^2 Density of all other spots in the channel besides the one modeled in the psf
covar_matrix = [pdf_sigma.^2, 0; 0, pdf_sigma.^2];
covar_inv = inv(covar_matrix);
covar_det = det(covar_matrix);

% Calc mesh x and y values
x_vector = [-flipud(distance_vector);distance_vector(2:end)];
y_vector = [flipud(distance_vector);-distance_vector(2:end)];
[x_mesh, y_mesh] = meshgrid(x_vector, y_vector);

% Calc psf of ch2
X = [x_mesh(:).';y_mesh(:).'];
X_shift_ch2 = X - repmat([ch2_seperation_length; 0], 1, length(X)); 
ch2_pdf = (2 .* pi .* sqrt(covar_det))^-1 .* exp(-0.5 .* sum(X_shift_ch2.' * covar_inv .* X_shift_ch2.', 2));

% Radially average to get the expected c(r) curve
% Loop through the desired indices
expected_radav_xcor = zeros(size(distance_vector));
parfor rad_ind = 1:length(distance_vector) 
    radius = distance_vector(rad_ind);
    
     % Special case where radius is 0
    if radius == 0
        x_circ = 0;
        y_circ = 0;
        num_uneval_points = 0;
    else

        % Spread out sampling grid evenly along circumference 
        % from 0 to pi, only calculate one half of the circle because the 
        % other half is symmetric. Don't calculate when centers of the
        % pdf are farther away than the analytical_calc_cultoff
        total_circ_length = 2 .* pi .* radius; % circumfrence in nm
        half_num_points = floor((total_circ_length./2)./analytical_arc_length_sampling) + 1;
        if half_num_points < 5; half_num_points = 5; end; % Ensure at least a few points at small radii when pixels are big
        total_num_points = (2 .* half_num_points) - 2;
        actual_interval = total_circ_length ./ total_num_points;

        % Decide which points to average.
        % Compute law of cosines expression a^2 + b^2 - c^2 / 2ab
        cosine_law_product = (radius.^2 + ch2_seperation_length.^2 -...
            (analytical_calc_cutoff .* pdf_sigma).^2) ./ (2 .* radius .* ch2_seperation_length);
        
        % Must re-define these outside of the if statements of in order for the parfor to run
        half_circ_lengths = [];
        repeat_last_circ_point = false;
        num_uneval_points = 0;
       
        % If the product is <-1, the two pdf centers are always within
        % the cutoff length and we need to calculate the entire half circle
        if cosine_law_product <= -1
            half_circ_lengths = linspace(0, total_circ_length ./ 2, half_num_points);
            repeat_last_circ_point = false;
            num_uneval_points = 0;

        % If the product is between -1 and 1, there exists an angle where the two centers
        % are past the cutoff and we need to go from 0 radians to the closest point to that angle
        elseif -1 < cosine_law_product && cosine_law_product < 1
            max_angle = acos(cosine_law_product);
            max_circumference = max_angle .* radius;
            num_half_circ_points = floor(max_circumference ./ actual_interval);
            half_circ_lengths = linspace(0, num_half_circ_points .* actual_interval, num_half_circ_points);
            repeat_last_circ_point = true;
            num_uneval_points = total_num_points - num_half_circ_points .* 2 - 1;

        % If the product is >1, there is no angle where the two centers
        % are less than the cutoff and we don't need to evaluate anything
        elseif 1 < cosine_law_product
            half_circ_lengths = [];
            num_uneval_points = total_num_points;
            repeat_last_circ_point = false; % Doesn't matter for an empty matrix
        end

        % Calculate radians and convert to x,y coordinates
        radians = half_circ_lengths./radius;
        [x_circ, y_circ] = pol2cart(radians, radius);
    end
   
    % Calculate psf of ch1 at the specified radial positions
    xcor_circ = zeros(length(x_circ), 1);
    for circ_ind = 1:length(x_circ)
        xshift = x_circ(circ_ind);
        yshift = y_circ(circ_ind);
        shift = [xshift; yshift];
        X_shift_ch1 = X - repmat(shift, 1, length(X)); 
        ch1_pdf = (2 .* pi .* sqrt(covar_det))^-1 .* exp(-0.5 .* sum(X_shift_ch1.' * covar_inv .* X_shift_ch1.', 2));

        % Calculate the xcor at the specified radial positions
        xcor = sum(ch1_pdf .* ch2_pdf) ./ sum(ch1_pdf .* ch2_total_density) +...
            sum(ch1_rand_density .* ch2_pdf) ./ (ch1_rand_density .* ch2_total_density .* length(ch1_pdf)) +...
            sum(ch1_pdf .* ch2_rand_density) ./ (ch1_total_density .* ch2_rand_density .* length(ch1_pdf)) +...
            (ch1_rand_density .* ch2_rand_density) ./ (ch1_total_density .* ch2_total_density);
        
        % Add to results vector
        xcor_circ(circ_ind) = xcor;
    end

    % Duplicate the 2nd half of the curve to save computation - not relevent for zero radius
    if radius == 0
        otherhalf_xcor = [];
    elseif repeat_last_circ_point
        otherhalf_xcor = xcor_circ(2:end); % Don't repeat 0 radians
    elseif ~repeat_last_circ_point
        otherhalf_xcor = xcor_circ(2:end-1); % Don't repeat 0 and pi radians
    end

    % Add both halves together
    xcor_circ = [xcor_circ; otherhalf_xcor];

    % Average the values
    eval_sum = sum(xcor_circ, 1); % Gives zero if xcor_circ is empty
    uneval_sum = num_uneval_points .* ((ch1_rand_density .* ch2_rand_density) ./ (ch1_total_density .* ch2_total_density));
    expected_radav_xcor(rad_ind) = (eval_sum + uneval_sum) ./ (length(xcor_circ) + num_uneval_points);
end
end

