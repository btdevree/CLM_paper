% Batch run script for overnight run

clear all
parpool(7)

FD_curves_script
clear all
cd ..
cd Direct_vs_Fourier

get_timings_script

delete(gcp)