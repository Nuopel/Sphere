function [ ] = comp3( h_sig1,h_sig2,h_sig3,t,a )
% comparison between 3 signals
if nargin==3
    a=1;
end
    
    

figure(a);subplot(311)
plot(t,h_sig1,'r')
hold on
plot(t,h_sig2,'k')
legend('Simu','Reel');xlim([0 0.015])
xlabel('Time [s]');ylabel('Amplitude ')
hold off

subplot(312)
plot(t,h_sig2,'k')
hold on
plot(t,h_sig3,'b')
legend('Reel','Virt');xlim([0 0.015])
xlabel('Time [s]');ylabel('Amplitude ')
hold off


subplot(313)
plot(t,h_sig3,'b')
hold on
plot(t,h_sig1,'r')
legend('Virt','Simu');xlim([0 0.015])
xlabel('Time [s]');ylabel('Amplitude ')
hold off

end

