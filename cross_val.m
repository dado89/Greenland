% This function computes a cross-validation measure of the goodness of the
% emulation fit carried out with correlation lengths d, on a set of n
% p-dim inputs and 1-dim outputs.
% Emulation is carried out on the dataset where one point, in turn, is 
% neglected and the prediction for that point is computed. The goodness of 
% the prediction is measured by evaluating a normal density with the mean 
% and variance provided by the emulator, on the known output y(j) that had 
% been neglected while building the emulator.
% The average of the log of such values is returned as global measure
% (thus, the higher this value, the better the fit).

% INPUTS
% - d:  p-dimensional vector with positive entries representing correlation
%       lengths along the p dimensions.
% - nu: nugget term, to introduce uncertainty in observed values
% - X:  nxp design matrix
% - y:  nx1 response vector
% - cor_fun: one of the strings 'exp2', 'matern32', 'matern52', 'abs_exp',
%            to specify which correlation function to use
% - method: one of the strings 'dens' or 'std', according to whether it is 
%           returned the pdf of the prediction evaluated at the left-out
%           point, or the number of standard deviations that the known
%           value of the left-out point is from the predicted mean.

% OUTPUT
% - resp: average of logarithm of normal densities, computed as explained
%         in the description above

% NOTE: uncomment the lines starting with '%%%' if just the difference
% between emulated prediction and real value is to be considered as
% goodness of fit measure.

function resp = cross_val(d, nu, X_full, y_full, cor_fun, method)

n=size(X_full,1);
H_full=[ones(n,1), X_full];
q=size(H_full,2);
dens=zeros(n,1);
val=zeros(n,1);

for j=1:n
    ind=[1:j-1, j+1:n];
    X=X_full(ind,:);
    H=H_full(ind,:);
    y=y_full(ind);
    h=H_full(j,:);  
    
    % Formulas for emulation
    A = Corr_fun(X, X, d, 0, cor_fun);
    A = A + nu*eye(size(A));
    K = H'/A; B = K*H; b = B\(K*y);
    f = y - H*b;  e = A\f;
    t = Corr_fun(X_full(j,:), X, d, 0, cor_fun);
    y_pred = h*b + t*e;
    
    s2 = (f'*e)/(n-1-q-2);  % there are n-1 design points
    p = h' - K*t';
    x = X_full(j,:);
    V = s2*( 1+nu - (t*(A\t')) + (p'*(B\p)) );
    S = sqrt(V);
    
    if strcmp(method, 'dens')
       % resp(j) = normpdf(y_full(j), y_pred, sqrt(V));
       deg=n-1-q;
       resp(j) = tpdf((y_full(j)-y_pred)/S, deg)/S;
    elseif strcmp(method, 'std')
        resp(j) = (y_pred-y_full(j))/S;
    else
        error('Last input must be one of ''dens'' or ''std''.');
    end
end

%resp = sum(log(dens))/n;
%resp=sqrt(sum(val.^2)/n);

end