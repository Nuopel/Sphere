function [var] = Pressure_map_(field,ct,Antenna,var,error)
% Calculate from Bmn coefficient the pressure map

if exist('error','var')
    opt=1;
else
    opt=0;
end

%% Reconstruction

Pmes_mat = reshape(field ,size(Antenna.X_mat)) ;


pcolor(Antenna.y,Antenna.x,real(Pmes_mat)) ;
shading interp

axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)
colorbar

if isfield(var, 'cax')==1
    caxis(var.cax);
else
    var.cax=caxis;
end
if opt==1
    caxis([0 100])
    r=ones(1,200)*5/(ct.k); theta=linspace(0,2*pi,200) ;
    [x ,y ] = pol2cart(theta,r);
    hold on
    plot(x,y,'--r')
end
end

