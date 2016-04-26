function [ ] = Pressure_map_(field,ct,error)
% Calculate from Bmn coefficient the pressure map 

if exist('error','var')
    opt=1;
else
    opt=0;
end

%% creation antenne
ct.pas_m = 2e-2; % pas de l'antenne
N.nbrx_sca =50; % nombre de micros par ligne
N.nbry_sca = 50; % nombre de micros par ligne
ct.N_mic=N.nbrx_sca*N.nbry_sca;
[ Antenna ] = AntennArray( ct.pas_m,N.nbrx_sca,N.nbry_sca) ;
%% Reconstruction

Pmes_mat = reshape(field ,size(Antenna.X_mat)) ;
figure

pcolor(Antenna.y,Antenna.x,real(Pmes_mat)) ;
shading interp
r=ones(1,200)*5/(ct.k); theta=linspace(0,2*pi,200) ;
[x ,y ]=pol2cart(theta,r);
hold on
plot(x,y,'--r')
axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)
colorbar
if opt==1
caxis([0 100])
end
end

