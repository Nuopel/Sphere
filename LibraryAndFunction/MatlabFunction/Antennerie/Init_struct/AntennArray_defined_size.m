function [ Antenna ] = AntennArray_defined_size(nbrx_sca,Size,offset)
%% [ Antenna ] = AntennArray_defined_size(nbrx_sca,Size,offset)
% Create a planar array structure of size Size, with a number of point per
% row and line nbrx_sca.
% It can be added an offset containing offset.x and offset.y
nbry_sca=nbrx_sca;
pas_m=Size/nbry_sca;

%% Option
if nargin==2
    offset.x=0;
    offset.y=0;
end

Antenna.N_mic=nbrx_sca*nbry_sca;
Antenna.x=-(pas_m*(nbrx_sca-1))/2:pas_m:(pas_m*(nbrx_sca-1))/2;Antenna.x=Antenna.x+offset.x; % mic position x

%case of negative points
if nbry_sca<0
    Antenna.y=-(pas_m*(nbry_sca+1))/2:-pas_m:(pas_m*(nbry_sca+1))/2;Antenna.y=Antenna.y+offset.y; % mic position y

else
Antenna.y=-(pas_m*(nbry_sca-1))/2:pas_m:(pas_m*(nbry_sca-1))/2-offset.y; % mic position y
end
[Antenna.Y_mat,Antenna.X_mat]=meshgrid(Antenna.y,Antenna.x);% create 2 matrix of coordinate associated to the grid
Antenna.coord_vect = [reshape(Antenna.X_mat,1,numel(Antenna.X_mat));reshape(Antenna.Y_mat,1,numel(Antenna.Y_mat))]; % transform the grid in a single row of pair coordinate of the mic

%% antenna in spherical coordinate
[ Antenna.theta, Antenna.phi,Antenna.R] = cart2sph( Antenna.coord_vect(1,:)', Antenna.coord_vect(2,:)',zeros(length(Antenna.coord_vect(2,:)),1) ) ;

%% Sort data with respect to the center (2D plot)
[Antenna.Rmicro_sort, Antenna.index ]= sort(Antenna.R);

%% Center in the middle of the index (2D plot)
for ii=1:nbrx_sca*nbry_sca;
    Antenna.Rmicro_sc(ii) = Antenna.Rmicro_sort(ii);
    Antenna.index_sc(ii) = Antenna.index(ii);
    if mod(ii,2) == 0
        Antenna.Rmicro_sc = circshift(Antenna.Rmicro_sc',[1,0]);Antenna.Rmicro_sc = Antenna.Rmicro_sc';
        Antenna.index_sc = circshift(Antenna.index_sc',[1 0]);Antenna.index_sc = Antenna.index_sc';    
        
    end
end

end