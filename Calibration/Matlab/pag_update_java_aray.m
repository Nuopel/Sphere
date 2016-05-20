function dblArray = pag_update_java_aray(M,sensitivity,mic_no,last_calibration,routing,h_sensitivity_table,dblArray)
% This updates the java array which shows the sensitivities.

for ii = 1:M
    dblArray(ii,1) = java.lang.Integer(mic_no(ii)); % Micro number
    dblArray(ii,2) = java.lang.Double(sensitivity(ii)); % Sensitivity
    dblArray(ii,3) = java.lang.String(last_calibration(ii,:)); % Calibration date
    dblArray(ii,4) = java.lang.Integer(routing(ii)); % Routing
end