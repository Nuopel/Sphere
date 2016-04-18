function [] = MovieMaker_double_erreur( h_sig,h_sig2,h_sig3,Antenna,begining,ending ,ct,text)
% MovieMaker take the signal and the Antenna set up to make a movie of the
% signal sent

% hH_SIG : signal used
% ANTENNA : structure of the antenna set up
% BEGINNING : first sample of the video
% ENDING : end sample of the video

if  nargin == 8
    opt = 1;
else
    opt = 0;
end



%% check size
% data
[a, b ]=size(h_sig);
if b>a
    sprintf('The matrix is transposed for fft used which work only with column of data')
    h_sig=permute(h_sig,[2 1 3]);
end

[a, b ]=size(h_sig2);
if b>a
    sprintf('The matrix is transposed for fft used which work only with column of data')
    h_sig2=permute(h_sig2,[2 1 3]);
end

[a, b ]=size(h_sig3);
if b>a
    sprintf('The matrix is transposed for fft used which work only with column of data')
    h_sig3=permute(h_sig3,[2 1 3]);
end

%% initialisation of the setting

figure(21)

subplot(311)
[~,pos]=max(h_sig(:,1)); % max position for scaling
grid_mat=reshape(h_sig(pos,:),size(Antenna.Y_mat));
pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
shading interp
cax = caxis;% register color setting to uniformise the movie

subplot(312)
[~,pos]=max(h_sig2(:,1)); % max position for scaling
grid_mat=reshape(h_sig2(pos,:),size(Antenna.Y_mat));
pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
shading interp

cax2 = caxis;% register color setting to uniformise the movie


%% initialisation of the setting

writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 96000/10000; % How many frames per second.
open(writerObj);
for ii=begining:ending
    clc;
    
    subplot(313)
    grid_mat_erreur=reshape(real(h_sig3(ii,:)),size(Antenna.X_mat));
    %        surf(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    pcolor(Antenna.y,Antenna.x,grid_mat_erreur); % plot of the frame
    shading interp
    hold on
    
    [C,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1]);
    hold off;    
    axis equal
    axis tight
    zlim([0 100])
    hold off
    
    
    subplot(311)
    grid_mat=reshape(real(h_sig(ii,:)),size(Antenna.X_mat));
    %        pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    surf(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    title(text.freq(ii))
    caxis(cax)
    shading interp
    axis equal
    axis tight
    set(gca,'zlim',[min(real(h_sig(:))) max(real(h_sig(:)))] )
    
    
    if opt==1
        xlabel(text.a)
    end
    
    
    subplot(312)
    grid_mat=reshape(real(h_sig2(ii,:)),size(Antenna.X_mat));
    %        pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    surf(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    caxis(cax2)
    shading interp

    set(gca,'zlim',[min(real(h_sig2(:))) max(real(h_sig2(:)))] )
    hold on
   [C,hfigc] = contour(Antenna.y,Antenna.x,grid_mat_erreur,[0 15]);
    set(hfigc, 'LineWidth',1.0,'Color', [1 1 1]);
    axis equal
    axis tight
    hold off
    if opt==1
        xlabel(text.b)
    end
    
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
    
end
close(writerObj); % Saves the movie.


end

