function [ ] = plot_micro( a,b,c,d )
if nargin==3
    d=1;
end
    

figure(d)
plot(a,1:56,'r+') 
hold on
plot(b,1:56,'+k') 
hold on
plot(c,1:56,'+b') 
legend('Simu','Reel','Virt')
hold off
xlabel('Delay in sample')
ylabel('Micropone number')

end
