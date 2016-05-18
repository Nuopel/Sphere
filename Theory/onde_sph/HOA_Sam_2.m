% Synthese d'une onde plane par HOA en base d'harmonique spherique
% Auteur : Dupont Samuel
% Version : 1.0 Fevrier 2016


close all; clear all; clc

%% SELECTION
M=5; % ordre de la decomposition

%% Choix de la source (Ae^j(wt-kx))

A=1;
f=1400; % frequency of the wave
w=2*pi*f;
c=343;k=w/c; % sound speed, wave number

%% POSITION DE LA SOURCE
rs=20; % distance de la source
pos_theta_deg=30; pos_theta_rad=deg2rad(pos_theta_deg); % Theta angle of the source
pos_phi_deg=0; pos_phi_rad=deg2rad(pos_phi_deg); % Phi angle of the source

[x_s,y_s,z_s] = sph2cart(pos_theta_rad,pos_phi_rad,rs);




%% SET UP
% N_hp = 2*M_sca+1; % Nombre de hp
N_hp=50; %lebedev grid 6,14,28,38,50...
r_hp_sca = 7; % rayon du cercle
[ x, y, z, w ] = ld_by_order(N_hp);
pos_hp_vect=[ x, y, z]*r_hp_sca;
[ theta_hp2, phi_hp2,~] = cart2sph(pos_hp_vect(:,1),pos_hp_vect(:,2),pos_hp_vect(:,3)); % conversion  to pherical coord possible problem with ref of cart2sph

H=zeros(1,M);
for ii=1:N_hp
    sigma=cos(phi_hp2(ii))*cos(pos_phi_rad)*cos(pos_theta_rad-theta_hp2(ii))+sin(phi_hp2(ii))*sin(pos_phi_rad);
    
    for jj=0:M
        S_int(jj+1)=HFm(jj,1,k,rs)/HFm(jj,1,k,r_hp_sca)*(2*jj+1)*legendreP(jj,sigma);
    end
    
    Sl(ii)=A*w(ii)*sum(S_int);
end

%% Antenne de mesure
pas_m = 0.01; % pas de l'antenne
nbrx_sca =100; % nombre de micros par ligne
nbry_sca =100; % nombre de micros par ligne
[ Antenna ] = AntennArray(pas_m,nbrx_sca,nbry_sca);


%% Pression de la source sur l'antenne
r_rs=sqrt(sum([Antenna.coord_vect(1,:)-x_s;Antenna.coord_vect(2,:)-y_s;zeros(1,length(Antenna.coord_vect))-z_s].^2));
Pmes_mat = reshape((exp(-1j*k*(r_rs))./(4*pi*r_rs)),size(Antenna.X_mat));


%% Pression reconstruite sur l'antenne on spherique
% k_mat = k.*[cos(azimut_hp_rad);sin(azimut_hp_rad)]'; % Les differents vecteurs d'ondes de chaques HP.
r_rhp2=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),pos_hp_vect(:,1)),bsxfun(@minus,Antenna.coord_vect(2,:),pos_hp_vect(:,2)),bsxfun(@minus,zeros(1,length(Antenna.coord_vect(1,:))),pos_hp_vect(:,3))).^2,3));
% concatene the 3 matrix |r-rs| --->[position_mic - coord_hp] with
% 3dimension being the different dimension x,y,z of r and rs for each positions
% of hp --->A=(n_hp,coord_mic,dimension) then do the sqrt(x^2+y^2+z^2) at all mic positon qnd for all hp
%
Precons_mat = sum(repmat(Sl',1,numel(Antenna.X_mat)).*exp(-1j*k*(r_rhp2))./(4*pi*r_rhp2));
Precons_mat =reshape(Precons_mat,size(Antenna.X_mat));

%% Erreur de reconstruction
erreur_mat=abs(Pmes_mat-Precons_mat)./abs(Pmes_mat)*100;


%% Affichage

figure(1)
subplot(1,2,1)
h=pcolor(Antenna.x,Antenna.x,real(Pmes_mat));
hcb = colorbar('location','EastOutside');
set(h, 'EdgeColor', 'none');
xlabel('Position x (m)')
ylabel('Position y (m)')
shading flat
cax = caxis;
% axis equal
% axis tight
hold on

subplot(1,2,2)

h=pcolor(Antenna.x,Antenna.x,real(Precons_mat));
shading flat
hcb = colorbar('location','EastOutside');
xlabel('Position x (m)')
ylabel('Position y (m)')
hold on
contour(Antenna.x,Antenna.x,erreur_mat,[0 7])
caxis(cax)
% axis equal
% axis tight
