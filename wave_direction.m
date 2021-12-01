function [im,theta,rho] = wave_direction(A1, A2, d, vel_tol)
%WAVE_DIRECTION Summary of this function goes here
%   Detailed explanation goes here

if (nargin == 3)
    vel_tol = 1e-13;
end

Im1 = double(A1>0);
Im2 = double(A2>0);

s = size(A1);
i_ids = repmat((1:s(1))',1,s(2));
j_ids = repmat(1:s(2),s(1),1);

%calculate gradient
[G1mag, G1dir] = imgradient(A1);
%components of the vector
G1i = -G1mag.*sin(G1dir/180*pi);
G1j = G1mag.*cos(G1dir/180*pi);

[G2mag, G2dir] = imgradient(A2);
%components of the vector
G2i = -G2mag.*sin(G2dir/180*pi);
G2j = G2mag.*cos(G2dir/180*pi);

vel_i = zeros(s);
vel_j = zeros(s);

for i_shift = -d:d
    for j_shift = -d:d

        C_err = abs(A2(d+1+i_shift:end-d+i_shift,d+1+j_shift:end-d+j_shift) - A1(d+1:end-d,d+1:end-d));
        G_err = abs(G2i(d+1+i_shift:end-d+i_shift,d+1+j_shift:end-d+j_shift) - G1i(d+1:end-d,d+1:end-d)) + ...
            abs(G2j(d+1+i_shift:end-d+i_shift,d+1+j_shift:end-d+j_shift) - G1j(d+1:end-d,d+1:end-d));

        M = 1./(C_err+G_err+1).^5;
        
        vel_i(d+1:end-d,d+1:end-d) = vel_i(d+1:end-d,d+1:end-d) - M.*(i_ids(d+1:end-d,d+1:end-d) - i_ids(d+1+i_shift:end-d+i_shift,d+1+j_shift:end-d+j_shift));
        vel_j(d+1:end-d,d+1:end-d) = vel_j(d+1:end-d,d+1:end-d) - M.*(j_ids(d+1:end-d,d+1:end-d) - j_ids(d+1+i_shift:end-d+i_shift,d+1+j_shift:end-d+j_shift));
        
    end
end

vel_i(abs(vel_i)<vel_tol) = 0;
vel_j(abs(vel_j)<vel_tol) = 0;

vel_i = vel_i.*Im1;
vel_j = vel_j.*Im1;

vel_mag = sqrt(vel_i.^2+vel_j.^2);
vel_i = vel_i./vel_mag;
vel_j = vel_j./vel_mag;

[theta,rho] = cart2pol(vel_j,-vel_i);

im = theta;
im = (im+pi)/(2*pi);
im(Im1 ==0) = 0;

end

