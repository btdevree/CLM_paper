% Script for obtaining Direct vs. Fourier correlation algorithm execution time comparisons

% Start parallel pool
parpool(8);

% Create number of points vector
number_points_binning = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 1e5; 3e5; 1e6; 3e6; 1e7; 3e7];
number_points_pdf_parallel = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 1e5; 3e5; 1e6; 3e6; 1e7; 3e7];
number_points_pdf = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 1e5; 3e5; 1e6; 3e6; 1e7];
number_points_direct = [1e2; 3e2; 1e3; 3e3; 1e4; 2.1e4];
number_points_memory = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 4.3e4];


% Run the algorithms; ignore memory usage, it's not really useful
[binning_Fourier_times, ~] = calc_Fourier_timing(number_points_vector, 'binning', true);
[pdf_parallel_Fourier_times, ~] = calc_Fourier_timing(number_points_vector, 'Gaussian_pdf_parallel', true);
[pdf_Fourier_times, ~] = calc_Fourier_timing(number_points_vector, 'Gaussian_pdf', true);
[pdist2_direct_times, ~] = calc_direct_timing(number_points_direct, 'pdist2', true);
[custom_memory_direct_times, ~] = calc_direct_timing(number_points_memory, 'custom_memory', true);


