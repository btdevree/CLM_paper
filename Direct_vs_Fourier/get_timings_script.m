% Script for obtaining Direct vs. Fourier correlation algorithm execution time comparisons

% Create number of points vector
number_points_vector = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 1e5; 3e5; 1e6; 3e6; 1e7];
number_points_direct = [1e2; 3e2; 1e3; 3e3; 1e4; 2.1e4];
number_points_memory = [1e2; 3e2; 1e3; 3e3; 1e4; 3e4; 4.2e4];

% Run the algorithms
[binning_Fourier_times, binning_Fourier_memory] = calc_Fourier_timing(number_points_vector, 'binning', true);
[pdf_Fourier_times, pdf_Fourier_memory] = calc_Fourier_timing(number_points_vector, 'Gaussian_pdf', true);
[pdist2_direct_times, pdist2_direct_memory] = calc_direct_timing(number_points_direct, 'pdist2', true);
[custom_memory_direct_times, custom_memory_direct_memory] = calc_direct_timing(number_points_memory, 'custom_memory', true);


