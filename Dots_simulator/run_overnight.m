% Overnight script

parpool(5)

tic

run_dots_correlation_no_over
magic(100)
sectime = toc();
mintime = sectime/60;
disp(['Time elapsed = ', num2str(mintime), ' minutes']);

run_dots_correlation_low_over
sectime = toc();
mintime = sectime/60;
disp(['Time elapsed = ', num2str(mintime), ' minutes']);

run_dots_correlation_high_over
sectime = toc();
mintime = sectime/60;
disp(['Time elapsed = ', num2str(mintime), ' minutes']);

delete(gcp)