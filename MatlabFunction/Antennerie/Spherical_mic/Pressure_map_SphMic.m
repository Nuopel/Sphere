function [ Pressure, var] = Pressure_map_SphMic(M,Bmn,ct,N,var,Antenna  )
% Calculate from Bmn coefficient the pressure map
% Samuel Dupont  may 2016

%% Calculate pressure from Bmn coefficient
Pressure = Decoding_pressure_field(M,Bmn,Antenna,ct,var,N ) ;
%% Reconstruction

Pmes_mat = reshape(Pressure	,size(Antenna.X_mat));
pcolor(Antenna.y,Antenna.x,real(Pmes_mat));
% cax=caxis;
shading interp
r=ones(1,200)*M/(ct.k); theta=linspace(0,2*pi,200);
[x ,y ]=pol2cart(theta,r);
hold on
plot(x,y,'--r')
xlabel('x')
axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)
colorbar
if isfield(var,'cax')==1
    caxis(var.cax);
else
    var.cax=caxis;
end
end
