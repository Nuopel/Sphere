% Synthese d'une onde spherique par HOA en base d'harmonique spherique
% influence de l'enlevement d'une source de reproduction puis ajout d'une source ambisoniaue pour compenser
% Auteur : Dupont Samuel
% Version : 1.0 Mai 2016

%% Define constants
ct.M = 5; % ordre de la decomposition
ct.Hanke_k_sca=2; % kind of hankel funtion, normally 2

% Choix de la source (Ae^j(wt-kx))
ct.A = 1;
ct.f = 1000; % frequency of the wave
ct.w = 2*pi*ct.f;
ct.c_air = 340;
ct.k=ct.w/ct.c_air; % sound speed, wave number

var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% POSITION DE LA SOURCE
source.rs = 5; % distance de la source
source.theta_deg=0; source.theta_rad = deg2rad(source.theta_deg); % Theta angle of the source
source.phi_deg=0; source.phi_rad = deg2rad(source.phi_deg); % Phi angle of the source
[source.x,source.y,source.z] = sph2cart(source.theta_rad,source.phi_rad,source.rs);% convertion cartesienne

%% SET UP
StrucRepro.N_hp_sca = 50; %lebedev grid 6,14,26,38,50...
StrucRepro.r_hp_sca = 1.07; % rayon du cercle
[ StrucRepro.x, StrucRepro.y, StrucRepro.z, StrucRepro.w ] = ld_by_order( StrucRepro.N_hp_sca );
StrucRepro.pos_hp_vect=[ StrucRepro.x, StrucRepro.y, StrucRepro.z ]*StrucRepro.r_hp_sca;
[ StrucRepro.theta, StrucRepro.phi,StrucRepro.r ] = cart2sph(StrucRepro.pos_hp_vect(:,1),StrucRepro.pos_hp_vect(:,2),StrucRepro.pos_hp_vect(:,3)); % conversion  to pherical coord possible problem with ref of cart2sph

%% Encodage: Decomposition en harmonique spherique de la source spherique
for ii=0:ct.M
    Encode.Fm((ii)^2+1:(ii+1)^2,1)=-HFm(ii,ct.Hanke_k_sca,ct.k,source.rs);
end
Encode.Ymn=sph_harmonic( ct.M,1,source.theta_rad,source.phi_rad);
Encode.Bmn=Encode.Ymn.*Encode.Fm;

%% Decodage
% Fabrication Matrice "C" spherical harmonic from position of the speaker
Decodage.Ymn_n3d_mat = sph_harmonic(ct.M,StrucRepro.N_hp_sca,StrucRepro.theta,StrucRepro.phi);
% Near field compensation (speakers acting as monopole)
for ii=0:ct.M
    Decodage.Fm((ii)^2+1:(ii+1)^2,1)=1./HFm(ii,2,ct.k,StrucRepro.r_hp_sca);
end
Decodage.S1=diag(StrucRepro.w)*Decodage.Ymn_n3d_mat.'*diag(Decodage.Fm)*Encode.Bmn; % decodage spherical wave

%% Decodage source 6
source2.theta = StrucRepro.theta(50);
source2.phi = StrucRepro.phi(50);
ct.r_hp_sca=StrucRepro.r_hp_sca;
Encode.Bmn2 = Bmn_monopole_encodage(ct.M, source2,ct,var ).';
% Fabrication Matrice "C" spherical harmonic from position of the speaker
Decodage.Ymn_n3d_mat_2 = sph_harmonic(ct.M,49,StrucRepro.theta(1:49),StrucRepro.phi(1:49));
% Near field compensation (speakers acting as monopole)

Decodage.S2=Decodage.S1(50)*diag(StrucRepro.w(1:49))*Decodage.Ymn_n3d_mat_2.'*diag(Decodage.Fm)*Encode.Bmn2; % decodage spherical wave

%%
Decodage.S(1:49,1)=Decodage.S1(1:49)+Decodage.S2;





%% __________________________Affichage______________________________________ 
%___________________________________________________________________________

%% Antenne de mesure
Antenna.pas_m = 0.01; % pas de l'antenne
Antenna.nbrx_sca =100; % nombre de micros par ligne
Antenna.nbry_sca =100; % nombre de micros par ligne
[ Antenna ] = AntennArray(Antenna.pas_m,Antenna.nbrx_sca,Antenna.nbry_sca);

%% Pression reconstruite sur l'antenne on spherique
figure
Affichage.r_rhp=sqrt(sum(cat(3,bsxfun(@minus,Antenna.coord_vect(1,:),StrucRepro.pos_hp_vect(1:49,1)),bsxfun(@minus,Antenna.coord_vect(2,:),StrucRepro.pos_hp_vect(1:49,2)),bsxfun(@minus,zeros(1,length(Antenna.coord_vect(1,:))),StrucRepro.pos_hp_vect(1:49,3))).^2,3));

% concatene the 3 matrix |r-rs| --->[position_mic - coord_hp] with
% 3dimension being the different dimension x,y,z of r and rs for each positions
% of hp --->A=(n_hp,coord_mic,dimension) then do the sqrt(x^2+y^2+z^2) at all mic positon qnd for all hp

Affichage.Precons_mat = sum(repmat(Decodage.S(1:49,1),1,numel(Antenna.X_mat)).*exp(-1j*ct.k*(Affichage.r_rhp))./(4*pi*Affichage.r_rhp));
Affichage.Precons_mat = reshape(Affichage.Precons_mat,size(Antenna.X_mat));
% Affichage.Precons_mat=Precons_mat2(end:-1:1,end:-1:1);


h=pcolor(Antenna.x,Antenna.x,real(Affichage.Precons_mat));
shading flat
hcb = colorbar('location','EastOutside');
xlabel('Position x (m)')
ylabel('Position y (m)')

 axis equal
 axis tight

