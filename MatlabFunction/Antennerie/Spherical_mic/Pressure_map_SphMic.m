function [ Pressure ] = Pressure_map_SphMic(Bmn,Sphmic,ct,N,var  )
% Calculate from Bmn coefficient the pressure map 

%% creation antenne
ct.pas_m = 2e-2; % pas de l'antenne
N.nbrx_sca =50; % nombre de micros par ligne
N.nbry_sca = 50; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca) ;

%% Calculate pressure from Bmn coefficient
 Pressure = Decoding_pressure_field(Bmn,Sphmic,Antenna,ct,var,N ) ;
%% Reconstruction

Pmes_mat = reshape(Pressure	,size(Antenna.X_mat));
figure

pcolor(Antenna.y,Antenna.x,real(Pmes_mat));
cax=caxis;
shading interp
r=ones(1,200)*5/(ct.k); theta=linspace(0,2*pi,200);
[x ,y ]=pol2cart(theta,r);
hold on
plot(x,y,'--r')
axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)

end

