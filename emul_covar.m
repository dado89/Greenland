% This function runs the emulator on a (small) set of inputs 'Test_Param'
% and returns both mean and covariance function for the emulated vector of
% length equal to the number of inputs given.

% INPUTS:
% - Test_Param: Nx8 matrix, in each row one input parameter for the
%               emulator (list of PC coefficients)
% - Design_Par: nx8 matrix, in each row one of the input parameters on
%               which the simulator has been run (list of PC coeffients)
% - Simul_Outputs: nx1 vector of known outputs corresponding to the inputs
%                  in Input_Runs
% - cor_fun: one of the strings 'exp2', 'matern32', 'matern52', 'abs_exp',
%            to specify which correlation function to use

% Typical Inputs
%{
Location=3;
T=readtable('Data/outputs.csv'); 
no_nan= ~isnan( table2array(T(:, end)) );
T=T(no_nan, 2:end);
Simul_Outputs = table2array(T(:, Location));
Design_par = csvread('Data/Coeff.csv');
Design_par = Design_par(no_nan,:);
clear no_nan
%}


function [mu,Sigma] = emul_covar(Test_Param, Design_par, Simul_Outputs, cor_fun)

n = size(Design_par,1); 
q = size(Design_par,2) + 1;
N = size(Test_Param,1);

y = Simul_Outputs;
d = 40*ones(1,q-1);
A = Corr_fun(Design_par, Design_par, d, cor_fun);    % nxn
H = [ones(n,1), Design_par]; % nxq
h = [ones(N,1), Test_Param]; % Nxq
K = H'/A;                    % qxn, K = H'*(A^-1)
B = K*H;                     % qxq, B = H'*(A^-1)*H
b = B\(K*y);                 % qx1, b = (B^-1)*Ky
f = y - H*b;                 % nx1
e = A\f;                     % nx1, e = A^-1 (y - Hb)
s2 = (f'*e)/(n-q-2);    % scalar

t = Corr_fun(Test_Param, Design_par, d, cor_fun); % Nxn
p = h - t*K'; % Nxq

mu = (h*b) + (t*e);
C = Corr_fun(Test_Param, Test_Param, d, cor_fun);
v1= t*(A\t');
v2= p*(B\p');
Sigma = s2*(C - v1 + v2);

% If round-off errors are present, overwrite them
S2=(Sigma+Sigma')/2;
if max(max(abs(Sigma-S2)))<1.e-15
    Sigma=S2;
end

end