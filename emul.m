% Loc = {'NEEM' 'NGRIP' 'GRIP' 'Renland' 'camp' 'DYE3' 'GISP2'};

% This function implements emulation on the data for Greenland orographies.
% Design_par is the set of design points (8-dimensional), Location is a
% string for the location where emulator is run, New_par is the set of new
% points where to emulate. The string retreat specifies whether retreat has
% to be taken into account during emulation or not.

% INPUTS:
% - Design_par: Nx8 matrix. In each row, an 8D input parameter for the
%               simulator/emulator (coefficients of PCs).
% - y:          Nx1 vector of simulator outputs corresponding to inputs in 
%               Design_par.
% - Total_par: Nx9 matrix. In the first 8 columns, each row is an 8D
%                  input where to run emulator (predictions). 
%                  9th column represents
%                  percentage of heat flow to be taken into
%                  consideration, expressed as a number in [0,1].
% - retreat: a string, either 'with_retreat' or 'no_retreat'. If
%            'no_retreat' is specified, the last column of Total_par will
%            not be used.
% - cor_fun: one of the strings 'exp2', 'matern32', 'matern52', 'abs_exp',
%            to specify which correlation function to use

% OUTPUTS:
% - M: Nx1 vector of emulated temperature mean for the parameters in
%      Total_par and the location 'Location'
% - S: Nx1 vector of emulated standard deviation

function [M, S, s2] = emul(Design_par, y, Total_par, retreat, cor_fun, cor_len, nu)

%% Check that input string retreat is OK
if ~strcmp(retreat, 'with_retreat') && ~strcmp(retreat, 'no_retreat')
    error('Third input must either the string ''with_retreat'' or the string ''no_retreat''.');
end

New_par = Total_par(:,1:(end-1));

%% Read the data for design points and corresponding response values

n=size(Design_par, 1); N=size(New_par,1);
q=size(Design_par,2)+1;

%% Code for emulator mean and variance

d = cor_len*ones(1,q-1);
A = Corr_fun(Design_par, Design_par, d, 0, cor_fun);    % nxn
A = A + nu*eye(size(A));

H = [ones(n,1), Design_par]; % nxq
h = [ones(N,1), New_par];    % Nxq
K = H'/A;                    % qxn, K = H'*(A^-1)
B = K*H;                     % qxq, B = H'*(A^-1)*H
b = B\(K*y);                 % qx1, b = (B^-1)*Ky
f = y - H*b;                 % nx1
e = A\f;                     % nx1, e = A^-1 (y - Hb)
s2 = (f'*e)/(n-q-2);         % scalar

%% Complete computation for emulator by dividing into for loops to save memory

N_block=10000;
N_loop=ceil(N/N_block);
t=zeros( N, size(Design_par,1) );

% Computation of t and p
for i=1:N_loop
    ind1 = (i-1)*N_block + 1;
    ind2 = min(i*N_block, N);
  %  t(ind1:ind2,:) = Corr_fun(New_par(ind1:ind2,:), Design_par, d, nu, cor_fun);     % Nxn
    t(ind1:ind2,:) = Corr_fun(New_par(ind1:ind2,:), Design_par, d, 0, cor_fun);     % Nxn
end


p = h - t*K'; % Nxq, p = h(x) - H'*A^-1*t


% Computation of final standard deviation:
% store v1=diag(t*(A\t')) and v2=diag(p'*(B\p)) through for loop

v1=zeros(N,1); v2=zeros(N,1); v3=zeros(N,1);
for i=1:N
    v1(i) = Corr_fun(New_par(i,:), New_par(i,:), d, 0, cor_fun) + nu;
    %v1(i) = 1;
    v2(i) = t(i,:)*(A\t(i,:)');
    v3(i) = p(i,:)*(B\p(i,:)');
end

M = (h*b) + (t*e);     % Nx1
S = s2*(v1 - v2 + v3);  % Nx1
S = sqrt(abs(S));


%% Now, if 'with_retreat' was specified, compute the correction factor

if strcmp(retreat, 'with_retreat')
    warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames');
    x = readtable('Data/data_sea_ice.xlsx', 'Range', 'B:B');
    x = table2array(x(2:end,1));
    y = readtable('Data/data_sea_ice.xlsx', 'Range', 'E:K');
    y = table2array(y(2:end,Location));
    shift = smoothing(x,y,25,300*Total_par(:,end)) - smoothing(x,y,25,0);
    M = M + shift;
end

end


