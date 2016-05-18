% Synthese d'une onde plane par HOA utilisant la decomposition en
% harmoniques cylindriques
% Auteur : Pierre Lecomte
% Version : 1.0 Juillet 2013
% D'après [Nicol, 1999. p. 86-88] Nicol. R. (1999). Restitution sonore
% spatialisee sur une zone etendue : application a la telepresence. These
% de doctorat


clear all
close all
clc

M_sca = 5; % ordre de la décomposition

%% Constante
c_m_s = 340; % celerite du son
t_sec = 0; % instant du calcul

%% Source 1
f_hz=500; % fréquence de l'onde monochromatique
a_sca = 1; % amplitude de l'onde monochromatique
k_sca =2*pi*f_hz/c_m_s; % vecteur d'onde

%% Direction de provenance de l'onde
azimut_deg = 20; % direction de provenance de l'onde en degré
azimut_rad = deg2rad(azimut_deg);
[k_vec(1),k_vec(2)] = pol2cart(deg2rad(azimut_deg),1);
k_vec = k_vec./norm(k_vec).*k_sca;


%% Dispositif de restitution
N_hp_sca = 2*M_sca+1; % nombre de hauts-parleurs pour la restitution
r_m = 0.5; % rayon du cercle
azimut_hp_rad = linspace(0,2*pi,N_hp_sca+1);
azimut_hp_rad = azimut_hp_rad(1:end-1);
[x_hp_m, y_hp_m] = pol2cart(azimut_hp_rad,r_m);

%% Décomposition en harmonique cylindrique [Nicol, 1999. p. 86]
ordre_vec = (1:M_sca);

% Vecteur u [Nicol, 1999. p. 86]
terme_cos_vec = sqrt(2).*(cos(ordre_vec*azimut_rad)');
terme_sin_vec = sqrt(2).*(sin(ordre_vec*azimut_rad)');
u_vec = ones(1,2*M_sca+1);
u_vec(2:2:end) = terme_cos_vec;
u_vec(3:2:end) = terme_sin_vec;
clear terme_cos_vec
clear terme_sin_vec

% Matrice U [Nicol, 1999. p. 86]
terme_cos_mat = sqrt(2).*cos(ordre_vec'*azimut_hp_rad)';
terme_sin_mat = sqrt(2).*sin(ordre_vec'*azimut_hp_rad)';
u_mat = ones(N_hp_sca,2*M_sca+1);
u_mat(:,(2:2:end)) = terme_cos_mat;
u_mat(:,(3:2:end)) = terme_sin_mat;
clear terme_cos_mat
clear terme_sin_mat

%% Gains de chaque HP [Nicol, 1999. p. 88]
G_vec = a_sca/N_hp_sca.*u_mat*u_vec';

%% Antenne de mesure
pasx_m = 0.01; % pas de l'antenne
nbm_sca = 101; % nombre de micros par ligne
x1_vec=-(pasx_m*(nbm_sca-1))/2:pasx_m:(pasx_m*(nbm_sca-1))/2; % coordonnees des micros
[X_mat,Y_mat]=meshgrid(x1_vec,x1_vec);
coord_vec = [reshape(X_mat,1,numel(X_mat));reshape(Y_mat,1,numel(Y_mat))];

%% Pression de la source sur l'antenne
Pmes_mat = a_sca*cos(k_vec*coord_vec - 2*pi*f_hz*t_sec);
Pmes_mat = reshape(Pmes_mat,size(X_mat));

%% Pression reconstruite sur l'antenne
k_mat = k_sca.*[cos(azimut_hp_rad);sin(azimut_hp_rad)]'; % Les différents vecteurs d'ondes de chaques HP.
Precons_mat = sum(repmat(G_vec,1,numel(X_mat)).*cos(k_mat*coord_vec - 2*pi*f_hz*t_sec),1);
Precons_mat = reshape(Precons_mat,size(X_mat));

%% Erreur de reconstruction [Trouver le bon critère, ici norme L2]
Erreur_mat = abs(Precons_mat - Pmes_mat) ./ abs(Pmes_mat)*100;
Erreur_mat = round(Erreur_mat/10)*10;

%% Vecteur vélocité
r_v_vec = [G_vec'./sum(G_vec)*cos(azimut_hp_rad)' G_vec'./sum(G_vec)*sin(azimut_hp_rad)'];


%% Vecteur énergie
r_e_vec = [abs(G_vec).^2'./sum(abs(G_vec).^2)*cos(azimut_hp_rad)' abs(G_vec).^2'./sum(abs(G_vec).^2)*sin(azimut_hp_rad)'];


%% Affichage

figure(1)
subplot(1,2,1)
pcolor(x1_vec,x1_vec,Pmes_mat)
xlabel('Position x (m)')
ylabel('Position y (m)')
shading interp
hold on
quiver(0,0,k_vec(1)/k_sca,k_vec(2)/k_sca,max(max(X_mat)))
subplot(1,2,2)
%plot(x_hp_m,y_hp_m,'o','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
hold on
pcolor(x1_vec,x1_vec,Precons_mat)
%pcolor(x1_vec,x1_vec,Erreur_mat)
xlabel('Position x (m)')
ylabel('Position y (m)')
shading interp
hold on
quiver(0,0,r_v_vec(1),r_v_vec(2),max(max(X_mat)))
quiver(0,0,r_e_vec(1),r_e_vec(2),max(max(X_mat)))
hold on
%contour(X_mat,Y_mat,Erreur_mat)%,'Color',[1,0,0]);
%h = findobj('Type','patch');
%set(h,'LineWidth',2)
hold off
