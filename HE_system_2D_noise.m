function HE_system_2D_noise(k0, s1, s2, DA, DI, DF, dx, dt, Im, A0, I0, F0, ...
        alf_AI, alf_F, A_ex_mcs, ex_t_vals, K, save_T, save_fold, func_fold)
% k0, s1, s2 [DOUBLE] - kinetic parameters of the model
% DA, DI, DF [DOUBLE] - diffusion coefficients
% dx, dt [DOUBLE] - spatial and demporal steps for the Euler difference
%   scheme
% Im [INT] - binary matrix, specifying the geometry of the domain, 
%   where the RD-system should be simmulated
% A0, I0, F0 [DOUBLE] - initial values of the concentrations
% alf_AI - magnitude of noise for A and I
% alf_F - magnitude of noise for F
% A_ex_mcs [CELL, DOUBLE] - cell array of A preturbations, that should be 
%   applied over the course of the simulations. The values of A_ex_mcs{i} 
%   will be added to the A variable. In order to have the mass conservation, 
%   I values will be corrected (distributed evenly): 
%   I = I - sum(A_ex(:))/sum(Im(:)); I(Im == 0) = 0;
% ex_t_vals [INT] - number of iteration, when the specified (from A_ex_mcs)
%   should be applied
% K [INT] - number of iterations to simmulate
% save_T [INT] - save the state of the system every save_T iteration 
%   results are saves in the folder 'A_tif' - images of A concentrations and
%   in the 'mat', where the matrices of all components (A, I and F) are
%   saved. Also the script saves at the beginnong the parameters of the 
%   simmulation in the file 'parameters.mat', located in the same directory
%   as 'A_tif' and 'mat'
% save_fold [STRING] - folder to save the output
% func_fold [STRING] - folder to add to the path, that contains supporting
%   functions (f, h and laplacian)

rng('shuffle');


addpath(func_fold);

s = size(Im);

% helper matrix to compute the laplacian with no flux boundary conditions
U = 0*Im;
U(2:(end-1),2:(end-1))=Im(3:end,2:(end-1))+Im(1:(end-2),2:(end-1))+...
    Im(2:(end-1),3:end)+Im(2:(end-1),1:(end-2));
U(Im==0)=0;

%specify initial conditions
A = A0; I = I0; F = F0;

%check perturbations
i = 0; %current iteration
for ii = find(ex_t_vals == i)
    disp(ii);
    A_ex = A_ex_mcs{ii};
    A = A + A_ex_mcs{ii};
    I = I - sum(A_ex(:))/sum(Im(:)); I(Im == 0) = 0;
end

mkdir(fullfile(save_fold, 'mat'));
mkdir(fullfile(save_fold, 'A_tif'));
mkdir(fullfile(save_fold, 'F_tif'));

save(fullfile(save_fold, 'parameters.mat'), 'k0', 's1', 's2', ...
        'DA', 'DI', 'DF', 'dx', 'dt', 'Im', 'A0', 'I0', 'F0', ...
        'A_ex_mcs','ex_t_vals', 'K', 'save_T', 'save_fold', 'func_fold');

for i=1:K
    %check perturbations
    for ii = find(ex_t_vals == i)
        A_ex = A_ex_mcs{ii};
        A = A + A_ex_mcs{ii};
        I = I - sum(A_ex(:))/sum(Im(:)); I(Im == 0) = 0;
    end
    
    fi = f(A,I,F,k0,s1,s2);
    noise = randn(s)*alf_AI;
    
    A_new = A + (fi + DA*laplacian(A,dx,U) + noise)*dt;
    A_new(Im==0)=0;
    A_neg = A_new;
    A_neg(A_neg > 0) = 0;
    A_new(A_new < 0) = 0;
    
    I_new = I + (-fi + DI*laplacian(I,dx,U) - noise)*dt;
    I_new(Im==0)=0;
    I_neg = I_new;
    I_neg(I_neg > 0) = 0;
    I_new(I_new < 0) = 0;
    
    F_new = F + (h(A,F) + randn(s)*alf_F)*dt;
    F_new(Im==0)=0;
    F_new(F_new < 0) = 0;
    
    A = A_new + I_neg;
    I = I_new + A_neg;
    F = F_new;

    if mod(i,save_T) == 0
        save(fullfile(save_fold, 'mat', strcat(num2str(i/save_T), '.mat')), 'A', 'I', 'F');
        imwrite(A, fullfile(save_fold, 'A_tif', strcat(num2str(i/save_T), '.tif')));
        imwrite(F, fullfile(save_fold, 'F_tif', strcat(num2str(i/save_T), '.tif')));
    end
end

end