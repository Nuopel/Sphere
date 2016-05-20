function dblArray = pag_update_array(M,sensitivity,mic_no,last_calibration,routing,dblArray)
% This updates the array which shows the sensitivities.

dblArray(:,1) = mat2cell(mic_no,ones(1,M),1); % Micro number
dblArray(:,2) = mat2cell(sensitivity,ones(1,M),1); % Sensitivity
dblArray(:,3) = mat2cell(last_calibration,ones(1,M),20); % Calibration date
dblArray(:,4) = mat2cell(routing,ones(1,M),1); % Routing