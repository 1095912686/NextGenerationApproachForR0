close all; clear all; clc;

%% SEIAR model
syms S E I A R N;
syms br dr beta kappa omegap omega p f;
syms gamma gammap;
 

%% Define Newly Infections and Transitions：

% all infected cells
variables = [E I A];

dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta*S*(I+kappa*A);
V(1) = dr*E + p*omegap*E + (1-p)*omega*E;
V(2) = -(1-p)*omega*E + (dr+f+gamma)*I;
V(3) = -p*omegap*E + (dr+gammap)*A;



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

%% compute its leading eigenvalue
% the matrix F*V^(-1)
invJV = inv(JV)
Mat = JF*invJV


% latexF = latex(F)
% latexV = latex(V)
% latexJF = latex(JF)
% latexJV = latex(JV)
% 
% latexInvJV = latex(invJV)

% eigenvalues
eigMat = eig(Mat)
%latexEigMat = latex(eigMat(end))


%% numerical experiment
S = 1e7;
beta = 9e-8;
omega = 1/5;
omegap = 1/7;
gamma = 1/3;
gammap = 1/5;
kappa = 0.8;
p = 0.3;
br = 0;
dr = 0;
f = 0;

R0 = eval(eigMat(end))