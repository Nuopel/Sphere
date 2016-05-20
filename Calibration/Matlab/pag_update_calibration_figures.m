function h_sensitivity_bar = pag_update_calibration_figures(sensitivity,mic_no,h_sensitivity)
% This updates figure h_sensitivity which shows the sensitivities.

sfigure(h_sensitivity);
h_sensitivity_bar = bar(mic_no,sensitivity);