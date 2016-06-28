function [] = plot_reel_virt_erreur(h_sig,h_sig2,h_sig3,Antenna,ii,text,d )

if exist('text','var')~=1
    opt=0;
else
    opt=1;
end

if exist('text','var')~=1
    d=1;
end


%% check size
% data
[a, b ]=size(h_sig);
if b>a
    sprintf('The matrix data is transposed  ')
    h_sig=permute(h_sig,[2 1 3]);
end

[a, b ]=size(h_sig2);
if b>a
    sprintf('The matrix data is transposed ')
    h_sig2=permute(h_sig2,[2 1 3]);
end

[a, b ]=size(h_sig3);
if b>a
    sprintf('The matrix data is transposed ')
    h_sig3=permute(h_sig3,[2 1 3]);
end
%% Plot part
r=ones(1,200)*5/(2*pi*ii/340); theta=linspace(0,2*pi,200);
[x ,y ]=pol2cart(theta,r);


hfig=figure(d);
set(hfig,'Position',[0 0 300 1000],'name',sprintf('%i',text.freq(ii)));
subplot(313)
grid_mat_erreur=reshape(h_sig3(ii,:),size(Antenna.X_mat));
pcolor(Antenna.y,Antenna.x,grid_mat_erreur); % plot of the frame
caxis([0 100])
shading interp
zlim([0 100]);xlabel('x [m]');ylabel('y [m]')
hold on
plot(x,y,'--r')
axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)

[~,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 15]);
set(hfigc, 'LineWidth',1.0,'Color', [1 1 1]);
[~,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 30]);
set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);

hold off
if opt==1
    title(text.c)
end


subplot(311)
grid_mat=reshape(h_sig(ii,:),size(Antenna.Y_mat));
% surf(Antenna.y,Antenna.x,grid_mat); % plot of the frame
pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
shading interp;xlabel('x [m]');ylabel('y [m]')
axis equal
axis tight
if opt==1
%     a=sprintf('%s frequence %i',text.a,ii+100);
    a=sprintf(text.a);
    title(a)
end

subplot(312)
grid_mat=reshape(h_sig2(ii,:),size(Antenna.Y_mat));
% surf(Antenna.y,Antenna.x,grid_mat); % plot of the frame
pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
cax=caxis;
shading interp
hold on

[~,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 15]);
set(hfigc, 'LineWidth',1.0,'Color', [1 1 1]);
[~,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 30]);
set(hfigc, 'LineWidth',1.0,'Color', [1 1 0]);
xlabel('x [m]');ylabel('y [m]')
caxis(cax)
plot(x,y,'--r')
axis equal
axis tight
axis([Antenna.y(1) Antenna.y(end) Antenna.x(1) Antenna.x(end)]	)
hold off
if opt==1
    title(text.b)
end





end

