% Compute the orthogonality error affilied to the lebedev quadrature
% Auteur : Dupont Samuel
% Version : 1.0 June 2016
clear variables; %close all; clc ;

% This program aim to show the orthogonality error affilied to the lebedev quadrature
% depending on the number of nodes
% D= Id - YtWY

%% Initialisation of the variables
ct.M=8;
ct.N_hp=50;%% number of nodes

%% Default variables 
ct.r_hp=1.07;
var.m_vect=0:ct.M;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;
%% Calculate  speakers positions
[ LoudspArray.x, LoudspArray.y, LoudspArray.z, LoudspArray.w ] = ld_by_order(ct.N_hp);
LoudspArray.pos_hp = [ LoudspArray.x, LoudspArray.y, LoudspArray.z ]*ct.r_hp;
[ LoudspArray.theta, LoudspArray.phi ,LoudspArray.r ] = cart2sph( LoudspArray.pos_hp(:,1),LoudspArray.pos_hp(:,2),LoudspArray.pos_hp(:,3));

%% Calculate Spherical harmonics
Ymn = sph_harmonic( ct.M,ct.N_hp,LoudspArray.theta,LoudspArray.phi );

%% Orthogonality matrix
D=diag(ones((ct.M+1)^2,1))-(Ymn*diag(LoudspArray.w)*Ymn.');


h=pcolor(abs(D));
 
set(h, 'EdgeColor', 'none');
hold on
for k = var.m_vect+1;
    x = [var.m_sum_vect(k)+1 var.m_sum_vect(k)+1];
    y = [0 var.m_sum_vect(end)];
    plot(x,y,'Color','k','LineStyle','-');
    plot(y,x,'Color','k','LineStyle','-');
end
colorbar

set(gca,'XTick',var.m_sum_vect)
set(gca,'XTickLabel',[var.m_vect])

set(gca,'YTick',var.m_sum_vect)
set(gca,'YTickLabel',[var.m_vect])
xlabel('Spherical harmonic order')
ylabel('Spherical harmonic order')
