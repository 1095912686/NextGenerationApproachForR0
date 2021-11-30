close all; clear all; clc;
format long;

%% SEIS model
syms S E I N;
syms br dr beta omega f gamma;


variables = [E I]; % all infected cells

% Define Newly Infections and Transitions
dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta*S*I;
V(1) = omega*E + dr*E;
V(2) = -omega * E + (dr + f + gamma) * I;



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

% eigenvalues
eigMat = eig(Mat)
latexEigMat = latex(eigMat(2))