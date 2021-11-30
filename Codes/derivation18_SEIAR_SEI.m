close all; clear all; clc;

%% SEIAR_SEI
syms S_p E_p I_p A_p R_p S_m E_m I_m;
syms br_p br_m dr_p dr_m beta_pp beta_pm beta_mp beta_mm;
syms kappa omega_18pPrime omega_18p omega_18m f_p gamma gammap;
syms a c n q f_p f_m;


variables = [E_p, I_p, A_p, E_m, I_m]; % all infected cells

% Define Newly Infections and Transitions
dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta_pp * S_p * (I_p + kappa*A_p) + beta_mp * S_p * I_m; 
F(4) = beta_mm * S_m * I_m + beta_pm * S_m *(I_p + A_p);
V(1) = dr_p * E_p + q * omega_18pPrime * E_p + (1-q) * omega_18p * E_p;
V(2) = -(1-q) * omega_18p * E_p + (gamma + dr_p + f_p) * I_p;
V(3) = -q * omega_18pPrime * E_p + (dr_p + gammap) * A_p;
V(4) = (omega_18m + dr_m) * E_m;
V(5) = -a*c*n*I_m - omega_18m * E_m + (dr_m + f_m) * I_m;



%% construct the next generation matrix 
% compute the jacobian matrix 
JF = sym(zeros(dim)); 
JV = sym(zeros(dim)); 
for i = 1:dim
    for j = 1:dim
        JF(i,j) = diff(F(i),variables(j));
        JV(i,j) = diff(V(i),variables(j));
    end
end
JF
JV

% the matrix F*V^(-1)
% the matrix is in fact the reproduction matrix for interactive
% transmissibility
Mat = JF*inv(JV)

M11 = latex(Mat(1:3,1:3))
M22 = latex(Mat(4:5,4:5))
M12 = latex(Mat(1:3,4:5))
M21 = latex(Mat(4:5,1:3))




