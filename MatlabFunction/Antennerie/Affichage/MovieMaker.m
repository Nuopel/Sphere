function [] = MovieMaker( h_sig,Antenna,begining,ending ,ct)
% MovieMaker take the signal and the Antenna set up to make a movie of the
% signal sent

% hH_SIG : signal used
% ANTENNA : structure of the antenna set up
% BEGINNING : first sample of the video
% ENDING : end sample of the video

%% check size
% data
[a, b ]=size(h_sig);
if b>a
    sprintf('The matrix is transposed for fft used which work only with column of data')
    h_sig=permute(h_sig,[2 1 3]);
end

%% initialisation of the setting


[~,pos]=max(h_sig(:,1)); % max position for scaling
grid_mat=reshape(h_sig(pos,:),size(Antenna.Y_mat));
pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
shading interp

cax = caxis;% register color setting to uniformise the movie


%% initialisation of the setting
figure(21)

writerObj = VideoWriter('out.avi'); % Name it.
writerObj.FrameRate = 96000/1000; % How many frames per second.
open(writerObj);
for ii=begining:ending
    clc;
    
    grid_mat=reshape(real(h_sig(ii,:)),size(Antenna.X_mat));
    %        pcolor(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    surf(Antenna.y,Antenna.x,grid_mat); % plot of the frame
    title(ii)
%     caxis(cax)
    shading interp
%     axis equal
%     axis tight
    set(gca,'zlim',[min(real(h_sig(:))) max(real(h_sig(:)))] )
    frame = getframe(gcf);
    writeVideo(writerObj, frame);
end
close(writerObj); % Saves the movie.


end

