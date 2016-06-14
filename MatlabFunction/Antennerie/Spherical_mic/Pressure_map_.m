function [var] = Pressure_map_(field,ct,var,error,Antenna)
% Calculate from Bmn coefficient the pressure map

if exist('error','var')
   if error==1
       opt=1;
   else
       opt=0;
   end
else
    opt=0;
end

if nargin==4 || nargin==3 
r=ones(1,200)*ct.M/(ct.k); 
Size=r(1)*2;
Antenna = AntennArray_defined_size(40,Size);
end

%% Reconstruction

Pmes_mat = reshape(field ,size(Antenna.X_mat)) ;


pcolor(Antenna.y,Antenna.x,real(Pmes_mat)) ;
shading interp
     r=ones(1,200)*ct.M/(ct.k); theta=linspace(0,2*pi,200);
    [x ,y ] = pol2cart(theta,r);
    hold on
    plot(x,y,'--r')
axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)
colorbar
xlabel('Position x [m]')
ylabel('Position y [m]')

if isfield(var, 'cax')==1
    caxis(var.cax);
else
    var.cax=caxis;
end
if opt==1
    hold on
    [~,hfigc] = contour(Antenna.y,Antenna.x,real(Pmes_mat),[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1 ]);
        [~,hfigc] = contour(Antenna.y,Antenna.x,real(Pmes_mat),[0 30]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
    caxis([0 100])

end
end

