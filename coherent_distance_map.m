function r = coherent_distance_map(theta, Ccr)

%compute coherent distance map
%the function computes for each pixel the largest radius of the circle
%where standard deviation of theta is less then Ccr

%INPUT
% theta - 2d array with values of theta
% Ccr - critical value of standard deviation

%OUTPUT
%r - 2d array, coherent distance map

r = zeros(size(theta));
[I,J] = ndgrid(1:size(theta,1),1:size(theta,2));
for x = 1:size(theta,2)
    for y = 1:size(theta,1)
        x1 = max(1,x-1); x2 = min(size(theta,2),x+1);
        y1 = max(1,y-1); y2 = min(size(theta,1),y+1);
        r_max = max(r(y1:y2, x1:x2),[],'all');
        D = (I-y).^2 + (J-x).^2;
        if r_max == 0
            s0 = 0;
            c = 0;
            while s0 <= Ccr
                c = c + 1;
                [~,s0] = circ_std(theta(D <= c^2));
            end
            r(y,x) = c;
        else
            c = r_max;
            [~,s0] = circ_std(theta(D <= c^2));
            if s0 <= Ccr
                while s0 <= Ccr
                    c = c + 1;
                    [~,s0] = circ_std(theta(D <= c^2));
                end
                r(y,x) = c;
            else
                while s0 > Ccr
                    c = c - 1;
                    [~,s0] = circ_std(theta(D <= c^2));
                end
                r(y,x) = c + 1;
            end
        end
    end
end
r(1,:)=0; r(end,:)=0; r(:,1)=0; r(:,end)=0;
end