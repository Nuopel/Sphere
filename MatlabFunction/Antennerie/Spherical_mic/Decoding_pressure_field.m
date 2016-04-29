function [ Pressure ] = Decoding_pressure_field(M,Bmn,Antenna,ct,var,N )



Ymn.Mic = sph_harmonic( M, ct.N_mic, Antenna.theta, Antenna.phi ) ;

Pressure=zeros(N.N_sweep,ct.N_mic);
for jj=1:ct.N_mic
    var.Bmn_Ymn = bsxfun(@times,Bmn,Ymn.Mic(:,jj)) ;
    
    % var.pressure = zeros(N.N_sweep,ct.M_th) ;------� tic toc method choice
    for ii=0:M
<<<<<<< HEAD
        var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(1:var.m_sum_vect(ii+1),:),1) ;
        var.Bessel_int = 1i^(ii)*(Bessel_sph(ii,ct.k.*Antenna.R(jj)));
=======
        var.Bessel_int = 1i^(ii)*(Bessel_sph(ii,ct.k.*Antenna.Rmicro(jj)));
>>>>>>> refs/remotes/origin/master
        if ii==0;
            
            
            var.sum_Bmn_Ymn = var.Bmn_Ymn(1,1) ;
            var.pressure_direc = permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Bessel_int ;
            
        else
            var.sum_Bmn_Ymn = sum(var.Bmn_Ymn(var.m_sum_vect(ii)+1:var.m_sum_vect(ii+1),:),1) ;
            var.pressure_direc = sum([var.pressure_direc, permute(var.sum_Bmn_Ymn,[2, 1, 3]).*var.Bessel_int],2) ;
            
        end
    end
<<<<<<< HEAD
    Pressure(jj)=var.pressure_direc;

=======
    Pressure(:,jj)=var.pressure_direc;
    
>>>>>>> refs/remotes/origin/master
    
end
% var.pressure(1,:)=0;% does it r


end

% test