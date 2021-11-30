close all; clear all; clc;
format long;

%% SIRC model
syms S I R C;
syms Gamma a p gamma mu beta;

variables = [I C]; % all infected cells

% Define Newly Infections and Transitions
dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta*S*(I+C);
V(1) = (a + gamma)*I;
V(2) = -p*gamma*I + (mu + a)*C;


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

% the matrix F*V^(-1)
invJV = inv(JV)
Mat = JF*invJV

% eigenvalues
eig(Mat)