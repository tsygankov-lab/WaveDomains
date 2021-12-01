function f = f(A, I, F, k0, s1, s2, A0, F0, gamma, delta)

if nargin == 6
    A0 = .4;
    F0 = .5;
    gamma = 1;
    delta = 1;
end

f = (k0 + gamma*A.^3./(A0^3+A.^3)).*I - delta*(s1+s2*F./(F0+F)).*A;
end