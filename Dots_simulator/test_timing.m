% test timing

% Load a parameter structure from current dirctory with deisred settings
load('dots_params.mat') % loads as 'params'

% Create a pdf map for the image
[testmap, ~] = dots_pdf_map(25, 100, 5, params);
cdf = cumsum(testmap(:));

% Create rand list
testrands = rands(1e7, 1);

tic
[~] = sample_inverse_cdf(testrands, cdf);
mattime = toc;

tic
[~] = sample_inverse_cdf_MEX(testrands, cdf);
mextime = toc;

disp(['MATLAB native time = ', num2str(mattime), ' sec, C++ MEX time = ', num2str(mextime), ' sec.']);
