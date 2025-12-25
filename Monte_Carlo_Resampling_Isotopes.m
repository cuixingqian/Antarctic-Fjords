% Monte Carlo Resampling used in He et al., Antarctic fjords paper 
% d13C, Fm, Frac_d13C, and Frac_Fm are bulk δ13C, bulk F14C, RPO fractional δ13C, and RPO fractional F14C, respectively.
clear; clc;
close all;

% import data
BD = readtable('file path of Antarctic RPO-14C dataset.xlsx','Sheet','bulk & RPO');
FD = readtable('file path of Antarctic RPO-14C dataset.xlsx','Sheet','RPO fracs');

d13C = table2array(BD(:,15));
Fm = table2array(BD(:,16));
Frac_d13C = table2array(FD(1:16,9));
Frac_Fm = table2array(FD(1:16,11));

% Estimate end-member contribution
% 10000 iterations
% initialization
iter = 10000;
Mar_d13C = -25+2*rand(iter,1);
Mar_Fm = 0.84+0.01*rand(iter,1);
Terr_d13C = -26.5+2.8*rand(iter,1);
Terr_Fm = 0.71+0.04*rand(iter,1);
Petro_d13C = -21.5+rand(iter,1);
Petro_Fm = zeros(iter,1);

p_Mar = zeros(28,iter);
p_Terr = zeros(28,iter);
p_Petro = zeros(28,iter);

% normalization by max-min approach
d13C_max = -22.2;
d13C_min = -24.71;
Fm_max = 0.837;
Fm_min = 0.287;

d13C_EM = [d13C; Frac_d13C];
Fm_EM = [Fm; Frac_Fm];

d13C_norm = (d13C_EM-d13C_min)./(d13C_max-d13C_min);
Fm_norm = (Fm_EM-Fm_min)./(Fm_max-Fm_min);
isotope_norm = [d13C_norm'; Fm_norm'];

% get the standard deviation of parameters
sigma1 = std(d13C_norm);
sigma2 = std(Fm_norm);
sigma = [sigma1; sigma2];

% normalization
Mar_d13C_norm = (Mar_d13C-d13C_min)./(d13C_max-d13C_min);
Mar_Fm_norm = (Mar_Fm-Fm_min)./(Fm_max-Fm_min);
Terr_d13C_norm = (Terr_d13C-d13C_min)./(d13C_max-d13C_min);
Terr_Fm_norm = (Terr_Fm-Fm_min)./(Fm_max-Fm_min);
Petro_d13C_norm = (Petro_d13C-d13C_min)./(d13C_max-d13C_min);
Petro_Fm_norm = (Petro_Fm-Fm_min)./(Fm_max-Fm_min);

Mar_isotope = [Mar_d13C_norm'; Mar_Fm_norm'];
Terr_isotope = [Terr_d13C_norm'; Terr_Fm_norm'];
Petro_isotope = [Petro_d13C_norm'; Petro_Fm_norm'];

% estimate endmember proportion by 10000 iterations
parfor i = 1:iter
    for j = 1:28
        objective = @(p) sqrt(sum(((p(1)*Mar_isotope(:,i)+p(2)*Terr_isotope(:,i)+p(3)*Petro_isotope(:,i)-isotope_norm(:,j))./sigma).^2));
        A_eq = [1, 1, 1]; b_eq = 1; c0 = [1/3, 1/3, 1/3];
        lb = [0, 0, 0]; ub = [1, 1, 1];
        options = optimoptions('fmincon', 'Display', 'final', 'Algorithm', 'sqp');
        [p_opt, D_min] = fmincon(objective, c0, [], [], A_eq, b_eq, lb, ub, [], options);

        p_Mar(j,i) = p_opt(1);
        p_Terr(j,i) = p_opt(2);
        p_Petro(j,i) = p_opt(3);
    end
end

% calculate the average values and standard deviations of three endmember proportions
p_Mar_mean = mean(p_Mar,2);
p_Mar_std = std(p_Mar')';
p_Terr_mean = mean(p_Terr,2);
p_Terr_std = std(p_Terr')';
p_Petro_mean = mean(p_Petro,2);
p_Petro_std = std(p_Petro')';