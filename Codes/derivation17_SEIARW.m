close all; clear all; clc;

%% SEIARW
syms beta beta_w p omega omegap mu mup gamma gammap epsilon kappa dr f;
syms S E I A R W N;


variables = [E, I, A, W]; % all infected cells

% Define Newly Infections and Transitions
dim = numel(variables);  
F = sym(zeros(dim,1)); 
V = sym(zeros(dim,1));

F(1) = beta*S*(I+kappa*A) + beta_w*S*W;
V(1) = (1-p)*omega*E + p*omegap*E + dr*E;
V(2) = -(1-p)*omega*E + gamma*I + dr*I + f*I;
V(3) = -p*omegap*E + gammap*A + dr*A;
V(4) = -mu*I - mup*A + epsilon*W;



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
% the matrix is in fact the reproduction matrix for interactive
% transmissibility
Mat = JF*inv(JV)

% Compute the eigenvalues of the matrix F*V^(-1)
[EigenVectors,EigenValues] = eig(Mat);
EigenValues = diag(EigenValues)


% take the maximum EigenValue and simplify
R0 = EigenValues(EigenValues ~= 0);   % given the possible maximum eigenvalues\
R0 = subs(R0,S,N);
latex(R0)

latexInvJV = latex(simplify(inv(JV)))

simplified_R0 = simplify(R0)
latexForSimplifiedR0 = latex(simplified_R0)


