function [var] = Pressure_map_gui(axe,field,ct,Antenna,var,error)
% Calculate from Bmn coefficient the pressure map

if exist('error','var')
    opt=1;
else
    opt=0;
end

%% Reconstruction

Pmes_mat = reshape(field ,size(Antenna.X_mat)) ;


pcolor(axe,Antenna.y,Antenna.x,real(Pmes_mat)) ;
shading(axe,'interp')
axis(axe, 'equal')
axis(axe,'tight')
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)
colorbar(axe)
if opt==1
    caxis(axe,[0 100])
    r=ones(1,200)*5/(ct.k); theta=linspace(0,2*pi,200) ;
    [x ,y ] = pol2cart(theta,r);
    hold on
    plot(axe,x,y,'--r')
end
if isfield(var, 'cax')==1
    caxis(axe,var.cax);
else
    var.cax=caxis(axe);
end
end

