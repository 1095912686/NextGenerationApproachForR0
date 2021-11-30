close all; clear all; clc;
format long;

%% SEIS model
syms S E I A N;
syms br dr beta kappa omegap omega p f;


variables = [E I A]; % all infected cells

% Define Newly Infections and Transitions
dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta*S*(I+kappa*A);
V(1) = dr*E + p*omegap*E + (1-p)*omega*E;
V(2) = -(1-p)*omega*E + (dr+f)*I;
V(3) = -p*omegap*E + dr*A;




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