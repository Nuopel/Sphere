function [ Sphmic,Pressure ] = Decoding_pressure_microphone( Bmn,Sphmic,N,ct,var )

[a, b ]=size(Bmn);
if b>a
    Bmn=Bmn.';
end

[a, b ]=size(ct.k);
if b>a
    ct.k=ct.k.';
end


[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Sphmic.x, Sphmic.y, Sphmic.z ) ;

% spherical harmonics of the microphone
Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta', Sphmic.phi' ) ;

Pressure.difract=zeros(N.N_sweep,ct.N_mic);%init
Fm=zeros(N.N_sweep,var.m_sum_vect(end));%init

for ii=0:ct.M_th
    var.Hprim(:,(ii)^2+1:(ii+1)^2) = repmat(1i^(ii-1)./((ct.k.*ct.r_micsph).^2.*Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)),1,var.nbr_m(ii+1)) ;%---> put out loop
end
    var.Hprim(isnan(var.Hprim(:)))=0+1i*0;

for jj=1:ct.N_mic
    clc;fprintf('Microphone %i',jj);
    var.Bmn_Ymn = bsxfun(@times,Bmn.source,Ymn.Mic(:,jj).') ;    
    
    %     for ii=0:ct.M_th
    %         test(ii+1)=var.sum_Bmn_Ymn;
    %         var.Hprim = 1i^(ii)*(Bessel_sph(ii,k.*ct.r_micsph)-Bessel_sph_1_deriv(ii,k.*ct.r_micsph)/Hankel_sph_1_deriv(ii,ct.hankel_order,k.*ct.r_micsph)*Hankel_sph(ii,ct.hankel_order,k.*ct.r_micsph));
    %         var.Hprim = 1i^(ii-1)/((ct.k.*ct.r_micsph).^2*Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)) ;%---> put out loop
    %         test4=var.Hprim.*var.Bmn_Ymn(var.m_sum_vect(ii):var.m_sum_vect(ii+1))
    %         test2(ii+1)=var.Hprim;
    Pressure.difract(:,jj)=sum(bsxfun(@times,var.Bmn_Ymn,var.Hprim),2);
    %         if ii==0;
    %             var.sum_Bmn_Ymn = var.Bmn_Ymn(1,:) ;
    %             var.pressure_diffr = permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Hprim ;
    %         else
    %             var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(var.m_sum_vect(ii)+1:var.m_sum_vect(ii+1),:),1) ;
    %             var.pressure_diffr = sum([var.pressure_diffr, permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Hprim],2) ;
    %
    %         end
    %     end
    %     Pressure.difract(:,jj)=var.pressure_diffr;
    %     Pressure.direct(:,jj)=var.pressure_direc;
    
    
end


end
