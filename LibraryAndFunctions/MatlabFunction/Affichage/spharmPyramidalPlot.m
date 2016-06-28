function spharmPyramidalPlot(L1, resolution)
% This function plots base spherical harmonic functions in their real forms
%
% Syntax:
% spharmPlot(L, resolution)
% spharmPlot(L)
% spharmPlot()
%
% Inputs:%
% L: degree of the functions, default value = 2;
% resolution: resolution of sphere surface, default value = 500.
%
% Results:
% red regions represent negative values
% red regions represent positive values
%
% Reference: http://en.wikipedia.org/wiki/Spherical_harmonics
%
% Modified version by Samuel Dupont from Mengliu Zhao, School of Computing Science, Simon Fraser University
% Update Date: 2016/may

if nargin < 1
    L1 = 2;
    resolution = 500;
end
if nargin < 2
    resolution = 500;
end

var.m_vect=0:L1;
var.m_sum_vect=(var.m_vect+1).^2;
var.nbr_m=(2.*var.m_vect)+1;

% discretize sphere surface
delta = pi/resolution;
theta = 0:delta:pi; % altitude
phi = 0:2*delta:2*pi; % azimuth
[phi,theta] = meshgrid(phi,theta);

% set figure background to white
figure('Color',[1 1 1])
L=0:L1;
for jj=1:length(L)
    M = L(jj):-1:-L(jj);
    for ii=1:length(M)
        % Legendre polynomials
        P_LM = legendre(L(jj),cos(theta(:,1)));
        P_LM = P_LM(abs(M(ii))+1,:)';
        P_LM = repmat(P_LM, [1, size(theta, 1)]);
        
        % normalization constant
        N_LM = sqrt((2*L(jj)+1)/4/pi*factorial(L(jj)-abs(M(ii)))/factorial(L(jj)+abs(M(ii))));
        
        % base spherical harmonic function
        if M(ii)>=0
            Y_LM = sqrt(2) * N_LM * P_LM .* cos(M(ii)*phi);
        else
            Y_LM = sqrt(2) * N_LM * P_LM .* sin(abs(M(ii))*phi);
        end
        
        % map to sphere surface
        r = Y_LM;
        x = abs(r).*sin(theta).*cos(phi);
        y = abs(r).*sin(theta).*sin(phi);
        z = abs(r).*cos(theta);
        
        % visualization
        subplot(L1+1,var.nbr_m(end),var.nbr_m(end)*(jj-1)+ round(var.nbr_m(end)/2)-floor(var.nbr_m(jj)/2)+ii-1);
        h = surf(x,y,z,double(r>=0));
        
        % adjust camera view
        view(40,30)
        camzoom(1.5)
        camlight left
        camlight right
        % 	lighting phong
        %
        axis([-1 1 -1 1 -1 1])
        
        % map positive regions to red, negative regions to green
        colormap(redbluecmap(2))
        
        % hide edges
        set(h, 'LineStyle','none')
        set(findobj(gcf, 'type','axes'), 'Visible','off')
        grid off
%         text(0,-.5,-1,['L = ', num2str(L(jj)), ', M = ', num2str(M(ii))])
    end
end
end