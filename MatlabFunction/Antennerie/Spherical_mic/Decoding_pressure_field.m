function [ Pressure ] = Decoding_pressure_field( Bmn,Sphmic,Antenna,ct,var,N ) 



[ Sphmic.theta, Sphmic.phi, Sphmic.r ] = cart2sph( Antenna.coord_vect(1,:)', Antenna.coord_vect(2,:)',zeros(length(Antenna.coord_vect(2,:)),1) ) ;
Ymn.Mic = sph_harmonic( ct.M, ct.N_mic, Sphmic.theta, Sphmic.phi ) ;

Pressure=zeros(N.N_sweep,ct.N_mic);
for jj=1:ct.N_mic
    var.Bmn_Ymn = bsxfun(@times,Bmn,Ymn.Mic(:,jj)) ;
    
    % var.pressure = zeros(N.N_sweep,ct.M_th) ;------� tic toc method choice
    for ii=0:ct.M
        var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(1:var.m_sum_vect(ii+1),:),1) ;
        var.Bessel_int = 1i^(ii)*(Bessel_sph(ii,ct.k.*Antenna.Rmicro(jj)));
        if ii==0;
            var.pressure_direc = permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Bessel_int ;

        else
            var.pressure_direc = sum([var.pressure_direc, permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Bessel_int],2) ;

        end
    end
    Pressure(:,jj)=var.pressure_direc;

    
end
% var.pressure(1,:)=0;% does it r


end

