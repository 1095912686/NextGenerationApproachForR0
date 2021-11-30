close all; clear all; clc;
format long;

%% SEIARV  model
syms S E I A R N V_1 V_2;
syms b_r d_r beta kappa omegap omega p f;
syms gamma gammap;

variables = [E I A]; % all infected cells

% Define Newly Infections and Transitions
dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta*S*(I+kappa*A) + beta*V_1*(I+kappa*A);
V(1) = d_r*E + p*omegap*E + (1-p)*omega*E;
V(2) = -(1-p)*omega*E + (d_r+f+gamma)*I;
V(3) = -p*omegap*E + (d_r+gammap)*A;




%% construct the next generation matrix and compute its leading eigenvalue.
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
invJV = inv(JV)
Mat = JF*invJV


latexF = latex(F)
latexV = latex(V)
latexJF = latex(JF)
latexJV = latex(JV)

latexInvJV = latex(invJV)

% eigenvalues
eigMat = eig(Mat)
latexEigMat = latex(eigMat(end))