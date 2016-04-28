function [ ] = plot_micro( a,d )
if nargin==3
    d=1;
end
    

figure(d)
plot(a,1:56,'r+') 
xlabel('Delay in sample')
ylabel('Micropone number')

end