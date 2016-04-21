function h = Fractional_delay_lagrange_matrix(N, delay)
%LAGRANGE  h=lagrange(N,delay) returns order N FIR
%          filter h which implements given delay
%          (in samples).  For best results,
%          delay should be near N/2 +/- 1.

% TODO implement matrix filter
% Need a delay column of vector
[a, b ]=size(delay);

if b>a
    delay=permute(delay,[2 1 3]);
end
[a, ~ ]=size(delay);

n = 0:N;
h = ones(a,N+1);
for k = 0:N
    index = find(n ~= k);
    h(:,index) = bsxfun(@times, h(:,index), bsxfun(@rdivide,(delay-k),(n(index)-k)));
end
end