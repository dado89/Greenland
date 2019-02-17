% This function first computes the emulators for the parameters specified in
% Input_par and the locations corresponding to index_loc. It then measures
% the compatibility between each emulator output and the data stored in
% data, which come with error bars too, through the following measure:
%
% I(x, loc) = (emul_mean(x,loc) - data(loc))/sqrt(emul_var(x,loc) + data_error(loc))
%
% where x is the input parameter (8-dim) and loc is the location being
% considered. These values are stored in the matrix X. The output
% index_compat is a logical vector detecting input parameters x for which
% I(x,loc)<c for all locations in index.
%
%
% INPUTS:
% - Design_par: Nx8 matrix with parameters of i-th design point at row i.
% - Outputs:    Nx7 vector of simulator outputs corresponding to inputs in 
%               Design_par and different locations (one for each column).
% - cor_fun: one of the strings 'exp2', 'matern32', 'matern52', 'abs_exp',
%            to specify which correlation function to use
% - d:       7x1 vector with correlation lenghts for each location
% - nu:      7x1 vector with nugget term for each location
% - Input_par: Nx9 matrix. In the first 8 columns, each row is an 8D
%              input where to run the emulator (predictions). 
%              9th column represents fraction (in [0,1]) of heat flow.
% - range: 3x7 matrix. In 2nd row, d18O values for all locations. In 1st
%          and 3rd row, minimum and maximum estimates of d18O values, 
%          respectively.
% - index_loc: vector of integers of length at most 7, specifying which 
%              locations should be emulated (index_loc=[1,2,3,5,6,7];)
% - thr_meas: vector of same length of index_loc, with positive values 
%             (typical 1 or 2) used as threshold to measure compatibility 
%             between data and emulator response. thr_meas(i) is the value
%             to be used for location index_loc(i).
% - retreat: a string, either 'with_retreat' or 'no_retreat'
%
% OUTPUTS
% - X: N x length(index_loc) matrix. X(i,j)= compatibility measure between
%      data and emulator prediction for input i, at location index_loc(j)
% - index_compat: logical column vector of length N with 1s at position i
%                 iff the compatibility measures for that input are less
%                 than c for all locations

% ORDER OF LOCATIONS
% 1: NEEM
% 2: NGRIP
% 3: GRIP
% 4: Renland
% 5: Camp
% 6: DYE3
% 7: GISP2

function [X, index_compat]=data_match(Design_par, Outputs, cor_fun, d, nu, Input_par, range, index_loc, thr_meas, retreat)

% From most input variables, select only locations that are specified in index_loc 
data=range(2,index_loc);
error_top= (range(3,index_loc) - range(2, index_loc))/sqrt(3);
error_bottom= (range(2,index_loc) - range(1, index_loc))/sqrt(3);
Outputs = Outputs(:, index_loc);
d = d(index_loc);
nu = nu(index_loc);

N=size(Input_par,1);
L=length(index_loc); % number of actual locations

% Carry out emulation on locations of interest
M=zeros(N,L);
S=zeros(N,L);
for i=1:L
    % Next line allows to keep visual track of progress
    disp(['Start of emulation for Location ' , num2str(i), '.']);
    y = Outputs(:,i);
    [M(:,i), S(:,i)] = emul(Design_par, y, Input_par, retreat, cor_fun, d(i), nu(i));
end

% Compute the compatibility quantity I(x, loc) specified in the preface
X = M - (ones(N,1)*data) ;

% Create a matrix with errors for the data, according to whether values in
% X are positive or negative (ie, atw emulator prediction is higher or
% lower than data value)

Err_data = zeros(size(S));
for i=1:L
    Err_data(:,i) = error_bottom(i);
    higher_prediction = X(:,i)>0;
    Err_data(higher_prediction, i) = error_top(i);
end

Var=S.^2 + Err_data.^2;
for i=1:L
    X(:,i)=X(:,i)./sqrt(Var(:,i));
end

% Detect which inputs give rise to max{I(x, loc_1), ..., I(x, loc_L)} < threshold
index_compat = true(N,1);
for i=1:L
    index_compat = index_compat & (abs(X(:,i)) < thr_meas(i));
end

end
