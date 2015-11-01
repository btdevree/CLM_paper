function [ parameter_structure ] = test_movie_parameters_dv
%TEST_MOVIE_PARAMETERS Creates a parameter structure template with some 
%   default options. For use with a dualview configuration.

% Documentation given inline with commands

% Initialize structure
params = struct();

% ------ Basic image and cell dimensions -----------

% Define Cartesian coordinate bounds in nanometers
% given as [min x, min y, max x, max y]
choice1 = [0, 0, 51200, 25600]; % Dualview
choice2 = [0, 0, 51200, 51200]; % Full frame
choice3 = [0, 0, 25600, 25600]; % Quarter frame
params.bounds = choice3;

% Define cell center in nanometers
choice1 = 'centered'; % Centered in the coordinate bounds at the beginning of the movie.
choice2 = [20000, 10000]; % Direct specification as [x, y]
params.cell_center = choice1;

% Define cell radius in nanometers
choice1 = '80%'; % x percent of the maximum radius that fits on the image at the beginning of the movie
choice2 = 10000; % Direct specification
params.cell_radius = choice1;

% ------ Cellular event distribution -----------------

% Define number of events in cell per movie
params.number_events_ch1 = 1000;
params.number_events_ch2 = 0;

% Define channel 1 autocorrelation model
choice1 = 'random'; % random distribution of events
params.ch1_autocor = choice1;
% NOTE: Would like to have better options for ch1 distribution
% (double/triple detections, "holes", etc.)

% Define channel 2 assignment distribution model
choice1 = 'random'; % ch2 events assigned to a ch1 event at random
choice2 = 'evenly_distributed'; % ch2 event assigned as evenly as possible to ch1 events
params.ch2_distribution = choice1;

% Define channel 2 assignment distribution model parameters
choice1 = [.8]; % fraction of ch2 points that are assigned to ch1 events
choice2 = [3]; % maximum number of ch2 points that are assigned to each ch1. Note: fractional values allowed (eg, 2.5 means 1/2 of ch1 events have 2 ch2 events assigned, and 1/2 have 3 events assigned).
params.ch2_distribution_params = choice1;

% Define channel 2 crosscorrelation model
choice1 = 'random'; % no relation to channel 1, random distribution
choice2 = 'exact'; % all channel 2 events are exactly a certain distance from a channel 1 event
choice3 = 'Gaussian'; % Gaussian distribution of channel 2 events from the channel 1 event centers
choice4 = 'critical'; % channel 2 events obey a 2D critical phase distribution
params.ch2_crosscor = choice2;

% Define channel 2 crosscorrelation model parameters
%choice1 = null % no parameters are needed, argument is ignored
choice2 = [0]; % distance between ch1 and ch2 events, in nanometers
choice3 = [0, 30]; % mean distance from ch1 event (mu) and standard deviation (sigma) given as [mu, sigma] in nanometers
%choice4 = ???? % TODO parameters for the critical model
params.ch2_crosscor_params = choice2;

% -------- Temporal event distribution ------------------

% Define event length distribution model
choice1 = 'random'; % Event can start and stop at any time, length obeys exponential distribution
choice2 = 'discrete'; % Events always start and stop at integer multiples of the frame times, length obeys discrete Poisson distribution
choice3 = 'one_frame'; % Events are always completed in a single frame
params.event_length_model = choice1;

% Define event length model parameters
choice1 = [1.3]; % define average lifetime (tau) of events in frames/event
choice2 = [0.8]; % define expected length of events in frames/event
%choice3 = null % no parameters are needed, argument is ignored
params.event_length_model_params = choice1;

% ------ Event parameters --------------

% Define average photon budgets; only constant flux models for now
params.ch1_photon_budget_params = [1500]; % Number of total photons detected for an average length channel 1 event
params.ch2_photon_budget_params = [1000]; % Number of total photons detected for an average length channel 2 event

% Define optical psf model
choice1 = 'Gaussian'; % Approximates the psf as a Gaussian distribution, calculated from the difference of the Gaussian cdf
choice2 = 'Airy'; % Calculates the psf using an ideal Airy disk from an isotropic emitter, integrated over idealized square ccd pixels
params.psf_model = choice1;

% Define optical psf model parameters for channel 1
choice1 = [69, 5]; % Gaussian model, given as [standard deviation(sigma) in nm, calculation limit in multiples of sigma]
choice2 = [200, 6]; % Airy disk model, given as [radius of the first dark ring in nm, calculation limit in multiples of the first ring radius]
params.psf_model_params_ch1 = choice1;

% Define optical psf model parameters for channel 2
choice1 = [79, 5]; % Gaussian model, given as [standard deviation(sigma) in nm, calculation limit in multiples of sigma]
choice2 = [230, 6]; % Airy disk model, given as [radius of the first dark ring in nm, calculation limit in multiples of the first ring radius]
params.psf_model_params_ch2 = choice1;

% Define resampling uncertainty model
choice1 = 'Gaussian_displacement'; % displacement of an exact psf by a zero-mean Gaussian distribution
choice2 = 'shot_noise'; % corruption of psf by noise in Poission distribution
choice3 = 'none'; % no resampling uncertainty
params.resampling_model = choice2;

% Define resampling uncertainty model parameters
choice1 = [20]; % standard deviation of gaussian displacement, in nm
%choice2 = null % no parameters are needed, argument is ignored
%choice3 = null % no parameters are needed, argument is ignored
params.resampling_model_params = choice1;

% ------- Background events ---------------------

% Define the number of random points to add to each channel
params.number_background_events_ch1 = 1000;
params.number_background_events_ch2 = 0;

% ------ Background -----------------------------

% Define overall background noise model
choice1 = 'Gaussian'; 
choice2 = 'Poisson';
params.background_noise_model = choice2;

% Define overall background noise parameters for channel 1
choice1 = [0, 0]; % Gaussian model, given as [mean(mu), standard deviation(sigma)] in detected photons/pixel
choice2 = [.5]; % Poisson model, given as [expected value(lambda)] in detected photons/pixel
params.background_noise_params_ch1 = choice2;

% Define overall background noise parameters for channel 2
choice1 = [0, 0]; % Gaussian model, given as [mean(mu), standard deviation(sigma)] in detected photons/pixel
choice2 = [1]; % Poisson model, given as [expected value(lambda)] in detected photons/pixel
params.background_noise_params_ch2 = choice2;

% Define cell background noise model
choice1 = 'Gaussian';
choice2 = 'Poisson'; 
params.cell_noise_model = choice2;

% Define cell background noise parameters for channel 1
choice1 = [0, 0]; % Gaussian model, given as [mean(mu), standard deviation(sigma)] in detected photons/pixel
choice2 = [5]; % Poisson model, given as [expected value(lambda)] in detected photons/pixel
params.cell_noise_params_ch1 = choice2;

% Define cell background noise parameters for channel 2
choice1 = [0, 0]; % Gaussian model, given as [mean(mu), standard deviation(sigma)] in detected photons/pixel
choice2 = [2]; % Poisson model, given as [expected value(lambda)] in detected photons/pixel
params.cell_noise_params_ch2 = choice2;

% ----- Camera and movie parameters ------------

% Define baseline clamp value
choice1 = [400]; % iXon EM Camera value
params.baseline_clamp = choice1;

% Define EM gain model
choice1 = 'Gaussian'; % Gaussian approximation of EM gain register (large signals)
choice2 = 'Gamma'; % Gamma function model of EM gain register (large and small signals)
choice3 = 'none'; % No EM gain
params.EM_gain_model = choice2;

% Define EM gain settings
params.EM_gain = 250;

% Define A/D conversion noise
choice1 = 50; %iXon EM at 10 MHz, 4.9x preamp gain, given in RMS e- 
choice2 = 38; %iXon EM at 5 MHz,  4.9x preamp gain, given in RMS e-
choice3 = 30; %iXon EM at 3 MHz,  4.9x preamp gain, given in RMS e-
choice4 = 6.5; %iXon EM at 1 MHz as conventional CCD, 2.4x preamp gain, given in RMS e- 
params.A_to_D_noise = choice2;

% Define A/D converstion ratio
choice1 = 11.7; %iXon EM at 10 MHz, 4.9x preamp gain, given in RMS e-/ADU
choice2 = 9.8; %iXon EM at 5 MHz,  4.9x preamp gain, given in RMS e-/ADU
choice3 = 9.9; %iXon EM at 3 MHz,  4.9x preamp gain, given in RMS e-/ADU
choice4 = 0.7; %iXon EM at 1 MHz as conventional CCD, 2.4x preamp gain, given in RMS e-/ADU
params.e_per_ADU = choice2;

% Define movie length
params.number_frames = 500;

% Define movie pixel size in nanometers
params.movie_pixel_size = 100;

% ----- STORM parameters-----------------------

% Define STORM image pixel size 
params.STORM_pixel_size = 7; % in nanometers

% Define STORM image localization precision
params.STORM_precision = 25; % in nanometers

% Define colors for the STORM display image
params.ch1_color = [0, 1, 0]; % green
params.ch2_color = [1, 0, 1]; % purple

% Flag to save movies to file
params.save_STORM_images_flag = false;

% ----- Cleanup ------------------------------

% Rename params to something more descriptive for the function return
parameter_structure = params;



end

