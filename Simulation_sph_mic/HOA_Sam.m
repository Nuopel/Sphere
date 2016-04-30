% Synthese d'une onde plane par HOA en base d'harmonique spherique
% Auteur : Dupont Samuel
% Version : 1.0 Fevrier 2016


close all; clear all; clc

%% SELECTION
M=5; % ordre de la decomposition
kk=1;

%% Choix de la source (Ae^j(wt-kx))

A=1;
f=1000; % frequency of the wave
w=2*pi*f;
c=340;k=w/c; % sound speed, wave number

%% POSITION DE LA SOURCE
rs=3.07; % distance de la source
pos_theta_deg=-45; pos_theta_rad=deg2rad(pos_theta_deg); % Theta angle of the source
pos_phi_deg=0; pos_phi_rad=deg2rad(pos_phi_deg); % Phi angle of the source

[x_s,y_s,z_s] = sph2cart(pos_theta_rad,pos_phi_rad,rs);

%% SET UP
N_hp_sca=50; %lebedev grid 6,14,26,38,50...
for M= [3 4 5 ]
    kk=kk+1;
    r_hp_sca =1.07; % rayon du cercle
    [ x, y, z, w ] = ld_by_order(N_hp_sca);
    pos_hp_vect=[ x, y, z]*r_hp_sca;
    
    
    %% Decomposition en harmonique spherique
    Hanke_k_sca=1;
    for ii=0:M
        Fm((ii)^2+1:(ii+1)^2)=FM_sph(ii,Hanke_k_sca,k,rs);
    end
    Ymn_mat=sph_harmonic( M,1,pos_theta_rad,pos_phi_rad);
    Bmn_mat=Ymn_mat.*Fm';
    clear terme_cos_vec; clear terme_sin_vec;
    %% Fabrication Matrice "C" spherical harmonic
    
    % [ theta_hp, phi_hp ] = xyz_to_tp ( x, y, z );
    [ theta_hp2_vect, phi_hp2_vect,r ] = cart2sph(pos_hp_vect(:,1),pos_hp_vect(:,2),pos_hp_vect(:,3)); % conversion  to pherical coord possible problem with ref of cart2sph
    Ymn_n3d_mat=sph_harmonic(M,N_hp_sca,theta_hp2_vect,phi_hp2_vect);
    
    %%  Decodage signal hp "C.S=Bmn"
    for ii=0:M
        
        Fm1_vect((ii)^2+1:(ii+1)^2,1)=1./HFm(ii,2,k,r_hp_sca);
    end
    
    Fm1_mat= diag(Fm1_vect);
    w= diag(w);
    % S2=w*Ymn_n3d_mat'*Bmn_mat; % decodage onde plane
    S2=w*Ymn_n3d_mat.'*Fm1_mat*Bmn_mat; % decodage spherical wave
    
    %% Antenne de mesure
    pasx_m = 0.006; % pas de l'antenne
    nbm_sca =200; % nombre de micros par ligne
    x1_vec=-(pasx_m*(nbm_sca-1))/2:pasx_m:(pasx_m*(nbm_sca-1))/2; % coordonnees des micros
    [X_mat,Y_mat]=meshgrid(x1_vec,x1_vec);% create 2 matrix of coordinate associated to the grid
    coord_vec = [reshape(X_mat,1,numel(X_mat));reshape(Y_mat,1,numel(Y_mat))]; % transform the grid in a single row of pair coordinate of the mic
    
    
    %% Pression de la source sur l'antenne
    r_rs=sqrt(sum([coord_vec(1,:)-x_s;coord_vec(2,:)-y_s;zeros(1,length(coord_vec))-z_s].^2));
    r_rs_slice=sqrt(sum([coord_vec(1,:)-x_s;zeros(1,length(coord_vec))-y_s;coord_vec(2,:)-z_s].^2));
    
    Pmes_mat = reshape((exp(-1j*k*(r_rs))./(4*pi*r_rs)),size(X_mat));
    Pmes_mat_slice = reshape((exp(-1j*k*(r_rs_slice))./(4*pi*r_rs_slice)),size(X_mat));
    
    
    %% Pression reconstruite sur l'antenne on spherique
    % k_mat = k.*[cos(azimut_hp_rad);sin(azimut_hp_rad)]'; % Les differents vecteurs d'ondes de chaques HP.
    r_rhp2=sqrt(sum(cat(3,bsxfun(@minus,coord_vec(1,:),pos_hp_vect(:,1)),bsxfun(@minus,coord_vec(2,:),pos_hp_vect(:,2)),bsxfun(@minus,zeros(1,length(coord_vec(1,:))),pos_hp_vect(:,3))).^2,3));
    r_rhp2_slice=sqrt(sum(cat(3,bsxfun(@minus,coord_vec(1,:),pos_hp_vect(:,1)),bsxfun(@minus,zeros(1,length(coord_vec(1,:))),pos_hp_vect(:,2)),bsxfun(@minus,coord_vec(2,:),pos_hp_vect(:,3))).^2,3));
    
    % concatene the 3 matrix |r-rs| --->[position_mic - coord_hp] with
    % 3dimension being the different dimension x,y,z of r and rs for each positions
    % of hp --->A=(n_hp,coord_mic,dimension) then do the sqrt(x^2+y^2+z^2) at all mic positon qnd for all hp
    
    Precons_mat2 = sum(repmat(S2,1,numel(X_mat)).*exp(-1j*k*(r_rhp2))./(4*pi*r_rhp2));
    Precons_mat2 =reshape(Precons_mat2,size(X_mat));
    Precons_mat2=Precons_mat2(end:-1:1,end:-1:1);
    
    Precons_mat2_slice = sum(repmat(S2,1,numel(X_mat)).*exp(-1j*k*(r_rhp2_slice))./(4*pi*r_rhp2_slice));
    Precons_mat2_slice =reshape(Precons_mat2_slice,size(X_mat));
    Precons_mat2_slice=Precons_mat2_slice(end:-1:1,end:-1:1);
    
    %% Pression reconstruite sur l'antenne onde plane
    % uhp_xmes_vect=sum(cat(3,bsxfun(@times,coord_vec(1,:),pos_hp_vect(:,1)/norm(pos_hp_vect(:,1))),bsxfun(@times,coord_vec(2,:),pos_hp_vect(:,2)/norm(pos_hp_vect(:,2))),bsxfun(@times,zeros(1,length(coord_vec(1,:))),pos_hp_vect(:,3)/norm(pos_hp_vect(:,3)))).^2,3);
    % %
    % for ii=1:N_hp_sca
    % r_rhp=sqrt(sum([coord_vec(1,:)*pos_hp_vect(ii,1)/norm(pos_hp_vect(:,1));coord_vec(2,:)*pos_hp_vect(ii,2)/norm(pos_hp_vect(:,2));zeros(1,length(coord_vec))*pos_hp_vect(ii,3)/norm(pos_hp_vect(:,3))].^2));
    %
    % if ii==1
    %     Precons_mat2 = S2(ii).*exp(-1j*k*(r_rhp));
    % else
    %     Precons_mat2 = Precons_mat2 + S2(ii).*exp(-1j*k*(r_rhp));
    % end
    %
    % end
    % Precons_mat2 = sum(repmat(S2,1,numel(X_mat)).*exp(-1j*k*(uhp_xmes_vect)),1);
    % Precons_mat2 =reshape(Precons_mat2,size(X_mat));
    
    %% Erreur de reconstruction
    erreur_mat=sqrt(abs(Pmes_mat-Precons_mat2).^2./abs(Pmes_mat).^2)*100;
    erreur_mat_slice=sqrt(abs(Pmes_mat_slice-Precons_mat2_slice).^2./abs(Pmes_mat_slice).^2)*100;

    %% Affichage
    if kk==2
        figure(1)
        subplot(2,2,1)
        pcolor(x1_vec,x1_vec,real(Pmes_mat));
        cax = caxis;

        hold on
        h=pcolor(x1_vec,x1_vec,real(Pmes_mat_slice));
        caxis(cax)
         rotate(h,[1,0,0],90)


        xlabel('Position x (m)')
        ylabel('Position y (m)')
        shading interp
        view(3)

    end
   
    subplot(2,2,kk)
    pcolor(x1_vec,x1_vec,real(Precons_mat2));
        caxis(cax)

    hold on
    h=pcolor(x1_vec,x1_vec,real(Precons_mat2_slice));
    rotate(h,[-1,0,0],90)
     caxis(cax)

    shading flat
    xlabel('Position x (m)')
    ylabel('Position y (m)')
    contour(x1_vec,x1_vec,erreur_mat,[0 15]);
    hold on
    [C,h2]=contour(x1_vec,x1_vec,erreur_mat_slice,[0 15]);
   
    xfm=hgtransform ;set(xfm,'parent',get(h2,'parent'));
    set(h2,'parent',xfm);
    R = makehgtform('xrotate',-pi/2);
    set(xfm,'Matrix',R)
    caxis(cax)
   view(3)
   
end
