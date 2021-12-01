function h = h(A, F, kn, ks, eps)

if nargin == 2
    kn = 1;
    ks = .25;
    eps = .1;
end

h = eps*(kn*A - ks*F);

end