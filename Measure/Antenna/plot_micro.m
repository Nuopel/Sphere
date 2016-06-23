function [ ] = plot_micro( a,b,c,d )
if nargin==3
    d=1;
end
    

figure(d)
plot(a,1:55,'r+') 
hold on
plot(b,1:55,'+k') 
hold on
plot(c,1:55,'+b') 
legend('Target','Decode','Simu')
hold off
xlabel('Delay in sample')
ylabel('Micropone number')

end