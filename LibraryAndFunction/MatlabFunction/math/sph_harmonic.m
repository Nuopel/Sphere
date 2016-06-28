function [ Ymn_n3d_mat ] = sph_harmonic( M,N_hp_sca,theta_hp2_vect,phi_hp2_vect )
% [ Ymn_n3d_mat ] = sph_harmonic( M,N_hp_sca,theta_hp2_vect,phi_hp2_vect )
% Calculate the spherical harmonics : ouput Orderxangle
% M: maximun order of the spherical harmonics, calculate all Ymn up to M
% N_hp_sca: number of points for which spherical harmonics are calculated
% phi_hp2_vect,theta_hp2_vect: coordinate of the points to calculates

Ymn_n3d_mat = cell(M+1,1); %declaration

terme_cos_mat=zeros(M+1,N_hp_sca);terme_sin_mat=zeros(M,N_hp_sca);% declaration

for ii=1:M+1
    
    terme_cos_mat(ii,:) = cos((ii-1)*theta_hp2_vect)';
    Pmn_mat=bsxfun(@times,(-1).^(-abs(0:ii-1)).',legendre(ii-1,sin(phi_hp2_vect)));
    
    if ii<2
        Ymn_int =Pmn_mat(ii,:);
    else
        terme_sin_mat(ii-1,:) = sin((ii-1)*theta_hp2_vect);
        N_norm=sqrt(diag(factorial(ii-1-(0:ii-1))./factorial(ii-1+(0:ii-1))).*(diag([1,2*ones(1,ii-1)])).*(2*(ii-1)+1));
        Pmn_mat=N_norm*Pmn_mat;
        Ymn_int = [ flipud( bsxfun(@times,Pmn_mat(2:end,:),terme_sin_mat(1:ii-1,:))) ; bsxfun(@times,Pmn_mat,terme_cos_mat(1:ii,:))];
        
    end
    
    Ymn_n3d_mat{ii}=Ymn_int;
    clear Ymn_int;
end
Ymn_n3d_mat=cell2mat(Ymn_n3d_mat);
end

