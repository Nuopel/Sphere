function [ Antenna ] = AntennArray(pas_m,nbrx_sca,nbry_sca)


Antenna.x=-(pas_m*(nbrx_sca-1))/2:pas_m:(pas_m*(nbrx_sca-1))/2; % mic position x
if nbry_sca<0
    Antenna.y=-(pas_m*(nbry_sca+1))/2:-pas_m:(pas_m*(nbry_sca+1))/2; % mic position y

else
Antenna.y=-(pas_m*(nbry_sca-1))/2:pas_m:(pas_m*(nbry_sca-1))/2; % mic position y
end
[Antenna.Y_mat,Antenna.X_mat]=meshgrid(Antenna.y,Antenna.x);% create 2 matrix of coordinate associated to the grid
Antenna.coord_vect = [reshape(Antenna.X_mat,1,numel(Antenna.X_mat));reshape(Antenna.Y_mat,1,numel(Antenna.Y_mat))]; % transform the grid in a single row of pair coordinate of the mic

Antenna.Rmicro = sqrt(Antenna.coord_vect(1,:).^2+Antenna.coord_vect(2,:).^2);
[ Antenna.theta, Antenna.phi,Antenna.R] = cart2sph( Antenna.coord_vect(1,:)', Antenna.coord_vect(2,:)',zeros(length(Antenna.coord_vect(2,:)),1) ) ;

[Antenna.Rmicro_sort, Antenna.index ]= sort(Antenna.R);

for ii=1:nbrx_sca*nbry_sca;
    Antenna.Rmicro_sc(ii) = Antenna.Rmicro_sort(ii);
    Antenna.index_sc(ii) = Antenna.index(ii);
    if mod(ii,2) == 0
        Antenna.Rmicro_sc = circshift(Antenna.Rmicro_sc',[1,0]);Antenna.Rmicro_sc = Antenna.Rmicro_sc';
        Antenna.index_sc = circshift(Antenna.index_sc',[1 0]);Antenna.index_sc = Antenna.index_sc';    
        
    end
end

end

