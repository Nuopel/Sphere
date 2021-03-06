function h = Fractional_delay_lagrange(N, delay)
%LAGRANGE  h=lagrange(N,delay) returns order N FIR
%          filter h which implements given delay
%          (in samples).  For best results,
%          delay should be near N/2 +/- 1.
% Implement matrix filter : Fractional_delay_lagrange_matrix(N, delay)
%
%Samuel Dupont
%May 2016
n = 0:N;
h = ones(1,N+1);
for k = 0:N
    index = find(n ~= k);
    h(index) = h(index) *  (delay-k)./ (n(index)-k);
end
end