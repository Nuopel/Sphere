function [ Antenna ] = AntennArray(pas_m,nbrx_sca,nbry_sca)


Antenna.x=-(pas_m*(nbrx_sca-1))/2:pas_m:(pas_m*(nbrx_sca-1))/2; % mic position x
if nbry_sca<0
    Antenna.y=-(pas_m*(nbry_sca+1))/2:-pas_m:(pas_m*(nbry_sca+1))/2; % mic position y

else
Antenna.y=-(pas_m*(nbry_sca-1))/2:pas_m:(pas_m*(nbry_sca-1))/2; % mic position y
end
[Antenna.Y_mat,Antenna.X_mat]=meshgrid(Antenna.y,Antenna.x);% create 2 matrix of coordinate associated to the grid
Antenna.coord_vect = [reshape(Antenna.X_mat,1,numel(Antenna.X_mat));reshape(Antenna.Y_mat,1,numel(Antenna.Y_mat))]; % transform the grid in a single row of pair coordinate of the mic



end

