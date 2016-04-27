function [ Sphmic,Pressure ] = Decoding_pressure_microphone( Bmn,Sphmic,N,ct,var )

[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Sphmic.x, Sphmic.y, Sphmic.z ) ;

Ymn.Mic = sph_harmonic( ct.M_th, ct.N_mic, Sphmic.theta', Sphmic.phi' ) ;

Pressure.difract=zeros(N.N_sweep,ct.N_mic);
Pressure.direct=zeros(N.N_sweep,ct.N_mic);


for jj=1:ct.N_mic
    var.Bmn_Ymn = bsxfun(@times,permute(Bmn.source,[2, 1, 3]),Ymn.Mic(:,jj)) ;
%     
%     var.Hprim = 1i^(var.m_vect-1)./((ct.k.*ct.r_micsph).^2*Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)) ;%---> put out loop      
%     var.Hprim(isnan(var.Hprim(:)))=0+1i*0;
    % var.pressure = zeros(N.N_sweep,ct.M_th) ;------» tic toc method choice
    for ii=0:ct.M_th
        var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(1:var.m_sum_vect(ii+1),:),1) ;
        test(ii+1)=var.sum_Bmn_Ymn;
%         var.Hprim = 1i^(ii)*(Bessel_sph(ii,k.*ct.r_micsph)-Bessel_sph_1_deriv(ii,k.*ct.r_micsph)/Hankel_sph_1_deriv(ii,ct.hankel_order,k.*ct.r_micsph)*Hankel_sph(ii,ct.hankel_order,k.*ct.r_micsph));
        var.Hprim = 1i^(ii-1)/((ct.k.*ct.r_micsph).^2*Hankel_sph_1_deriv(ii,ct.hankel_order,ct.k.*ct.r_micsph)) ;%---> put out loop      
        test2(ii+1)=var.Hprim;
        var.Hprim(isnan(var.Hprim(:)))=0+1i*0;

        var.Bessel_int = 1i^(ii)*(Bessel_sph(ii,ct.k.*ct.r_micsph)) ;
        if ii==0;
            var.pressure_diffr = permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Hprim ;
            var.pressure_direc = permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Bessel_int ;
            test3=var.pressure_diffr;
        else
            var.pressure_diffr = sum([var.pressure_diffr, permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Hprim],2) ;
            test3(ii+1)=permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Hprim;
            var.pressure_direc = sum([var.pressure_direc, permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Bessel_int],2) ;

        end
    end
    Pressure.difract(:,jj)=var.pressure_diffr;
    Pressure.direct(:,jj)=var.pressure_direc;

    
end


end

