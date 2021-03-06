function [ TH ] = wavefrontarrow( Pressure,Antenna,opt )
%% wavefrontarrow(Pressure,Antenna,opt) compute an arrow representing the direction of the mean intensity
%
%INPUT: the pressure field PRESSURE, 
%        the structure ANTENNA associated (see AntennArray_defined_size.m),
%        An option for the arrow color, string like 'r' red, 'b' blue ...
%
% OUTPUT: The angle theta TH, of the mean arrow
%
% EXAMPLE:
%               wavefrontarrow( Pressure.p_recons_target,Antenna );
%               ct.theta_simu=wavefrontarrow( Pressure.test,Antenna,'r');
%               

if nargin<3
    opt=0;
end

Pressure=reshape(Pressure,size(Antenna.X_mat));

[dx,dy]=gradient(Pressure,abs(Antenna.x(1)-Antenna.x(end))/length(Antenna.x),abs(Antenna.x(1)-Antenna.x(end))/length(Antenna.x));
dx=-1j*conj(dx);
dy=-1j*conj(dy);
P = bsxfun(@times,[dx(:), dy(:)],Pressure(:));
Px=reshape(P(:,1),size(Antenna.X_mat));
Py=reshape(P(:,2),size(Antenna.X_mat));
% quiver(Antenna.x,Antenna.y,Px,Py)
% contour(Antenna.x,Antenna.y,real(Pressure.test))
hold on

ct.arrownorm=Antenna.x(end)/norm([mean((Px(:))) mean((Py(:)))])*0.5;
arrow=real([mean((Px(:))) mean((Py(:)))].*ct.arrownorm);
if opt==0
quiver(0,0,arrow(1),arrow(2),0, 'MaxHeadSize',1.90,'color','k')
else
quiver(0,0,arrow(1),arrow(2),0, 'MaxHeadSize',1.90,'color',opt)
end    
[TH,~] = cart2pol(arrow(1),arrow(2)) ;
% quiver(Antenna.x,Antenna.y,dx,dy)


end

