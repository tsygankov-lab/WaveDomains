clear; clc;

s = [202, 202];

k0 = 0.2; s1 = 0.5; s2 = 1.2;
DA = .001/3; DI = .1/3; DF = 0;
dx = 0.02; dt = .001;
alf_AI = 3; alf_F = 0.1;
Im = ones(s); Im(1,:) = 0; Im(end,:) = 0; Im(:,1) = 0; Im(:,end) = 0;
A0 = zeros(s); A0 = A0.*Im;
I0 = ones(s); I0 = I0.*Im;
F0 = zeros(s);
K = 4000000;
save_T = 1000;

ex_t_vals = [];
A_ex_mcs = {};

save_root_fold = 'Z:\Siarhei Hladyshau\Integrative_Model_of_Cell_Morphodynamics\actin_cortical_waves\noise_triggering\AI_correct_mass_concervation\square_geom_different_noise_values_for_IA_and_F';
str_pars = strcat('k0_', num2str(k0), '__s1_', num2str(s1), '__s2_', num2str(s2), ... 
    '_alf_AI_', num2str(alf_AI), '_alf_F_', num2str(alf_F), '__size_', num2str(s(1)), 'x', num2str(s(2)));
save_fold = fullfile(save_root_fold, str_pars);
func_fold = 'Z:\Siarhei Hladyshau\Integrative_Model_of_Cell_Morphodynamics\actin_cortical_waves\noise_triggering\AI_correct_mass_concervation\square_geom_different_noise_values_for_IA_and_F';

tic;
HE_system_2D_noise(k0, s1, s2, DA, DI, DF, dx, dt, Im, A0, I0, F0, ...
                   alf_AI, alf_F, A_ex_mcs, ex_t_vals, K, save_T, save_fold, func_fold);
toc;
