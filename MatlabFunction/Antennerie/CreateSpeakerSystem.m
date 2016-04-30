function [ ArraySpeaker, N ] = CreateSpeakerSystem(r_hp_sca)
%UCreateSpeakerSystem : 
% - Import the data of the position of the speakers 
%   on an ambisonics set up of 50 loudspeakers,microphone
% - Rotate the coordinates of 45 degree in order to 
%   be in in the same way than the physical set up 
% - Set the sphere radius at r_hp_sca
% - Use CreateSpeakerStructure function

<<<<<<< HEAD
imp=load('data/coords.mat');imp.coords(:,1:3)=imp.coords(:,1:3)*r_hp_sca;%imp.coords(:,2)=-imp.coords(:,2) ;
=======
imp=load('data/coords.mat');imp.coords(:,1:3)=imp.coords(:,1:3)*r_hp_sca;imp.coords(:,2)=imp.coords(:,2) ;
>>>>>>> refs/remotes/origin/master
% A=rotz(45*pi/180) ;
% for ii=1:length(imp.coords)
% imp.coords(ii,1:3)=imp.coords(ii,1:3)*A' ;
% end
[ArraySpeaker] = CreateSpeakerStructure( imp.coords(:,1), imp.coords(:,2), imp.coords(:,3), imp.coords(:,4) );
N=length(imp.coords(:,3));

end

