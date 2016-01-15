% Overnight script

parpool(6)

run_dots_correlation_no_over
run_dots_correlation_low_over
run_cots_correlation_high_over

delete(gcp)