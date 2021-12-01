clear; clc;

rng shuffle;

N = 3;
model = createpde(N);

% Coordinates
lowerLeft  = [0, 0];
lowerRight = [4, 0];
upperRight = [4, 4];
upperLeft =  [0, 4];
% Geometry matrix
S = [3,4 lowerLeft(1), lowerRight(1), upperRight(1), upperLeft(1), ...
         lowerLeft(2), lowerRight(2), upperRight(2), upperLeft(2)];                     
gdm = S';
% Names
ns = 'S';
% Set formula 
sf = 'S';
% Invoke decsg
g = decsg(gdm,ns,sf');
geometryFromEdges(model,g);

precision = 0.05;
%precision = 0.2;
generateMesh(model,'Hmax',precision);

[p,e,t] = meshToPet(model.Mesh);
 
DA = 0.001/3;
DI = 0.1/3;
DF = 0;

gamma = 1;
A0 = 0.4;
delta = 1;
s1 = 0.5;
F0 = 0.5;
kn = 1;
ks = 0.25;
eps = 0.1;

k0 = 0.2;
s2 = 0.8;

c = [DA;DA;DI;DI;DF;DF];
f = @(region,state) f_fun(region, state, N, k0, gamma, A0, delta, s1, s2, F0, kn, ks, eps);
ic = @(location) ic_rand_fun(location,N);
%ic = @(location) ic_spike_fun(location,N);

specifyCoefficients(model, 'm', 0, 'd', 1, 'c', c, 'a', 0, 'f', f);

clampedBC1 = applyBoundaryCondition(model,'neumann','Edge',1,'g',[0 0 0]);
clampedBC2 = applyBoundaryCondition(model,'neumann','Edge',2,'g',[0 0 0]);
clampedBC3 = applyBoundaryCondition(model,'neumann','Edge',3,'g',[0 0 0]);
clampedBC4 = applyBoundaryCondition(model,'neumann','Edge',4,'g',[0 0 0]);
 
setInitialConditions(model,ic);

time = 0:1:4000;
%time = 0:10:1000;
%time = 0:1:100;

tic;
result = solvepde(model,time);
toc;

fold = 'Z:\Siarhei Hladyshau\Integrative_Model_of_Cell_Morphodynamics\actin_cortical_waves\FEM_Matlab';
save_file = strcat('k0_', num2str(k0), '_s2_', num2str(s2), '.mat');

save(fullfile(fold, save_file), 'model', 'precision', 'DA', 'DI', 'DF', ...
    'k0', 'gamma', 'A0', 'delta', 's1', 's2', 'F0', 'kn', 'ks', 'eps', 'time', 'result', '-v7.3');

u = result.NodalSolution;

fig = figure('Position', [50 500 1000 200]);
for i = 1:length(time)
    disp(i);
    clf;
    subplot(1,3,1);
    hold on;
    pdeplot(model,'XYData',u(:,1,i),'ColorMap','gray');
    title('A');
    
    subplot(1,3,2);
    pdeplot(model,'XYData',u(:,2,i),'ColorMap','gray');
    title('I');
    
    subplot(1,3,3);
    pdeplot(model,'XYData',u(:,3,i),'ColorMap','gray');
    title('F');
    drawnow;
end

disp('done');


function f = f_fun(region, state, N, k0, gamma, A0, delta, s1, s2, F0, kn, ks, eps)

    nr = length(region.x); % Number of columns
    f = zeros(N,nr); % Allocate f
    
    u = state.u;
    
    F1 = (k0 + gamma*u(1,:).^3./(A0^3+u(1,:).^3)).*u(2,:) - delta*(s1+s2*u(3,:)./(F0+u(3,:))).*u(1,:); 
    F2 = eps*(kn*u(1,:)-ks*u(3,:));
    
    f(1,:) = F1;
    f(2,:) = -F1;
    f(3,:) = F2;
    
end

function ic = ic_rand_fun(location,N)
    M = length(location.x);
    noise = 0.000001*abs(randn(1,M));
    ic = zeros(N,M);
    ic(1,:) = noise;
    ic(2,:) = 1-noise;
    ic(3,:) = 0.000001*abs(randn(1,M));
end

function ic = ic_spike_fun(location,N)
    M = length(location.x);
    ic = zeros(N,M);
    ids = (location.x<2) & (location.y<2);
    ic(1,:) = 0;
    ic(1,ids) = 0.1;
    ic(2,:) = 1;
    ic(2,ids) = ic(2,ids) - 0.1;
    ic(3,:) = 0;
end