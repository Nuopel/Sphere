function [ pos_sca ] = closest( select,vect )
% Give back the closest position of the value considered SELECT from the
% vector VECT

[~, pos_sca] = min(abs(vect-select)); %index of closest value


end
