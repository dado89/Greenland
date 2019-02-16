% This function finds Principal Components from the shapes stored in the
% file 'Data/Dataset.txt'. Weights are used to account for spherical
% geometry of Earth.

% OUTPUTS
% - PC: (Nlat)x(Nlon)x n matrix with principal components. PC(:,:,j) is
%       the j-th Principal Component.
% - M:  (Nlat)x(Nlon) matrix, with cell-by-cell average of original shapes.
% - Std: column vector with standard deviations of Principal Components


function [PC, M, Std]=pca_greenland()

mask=ncread('Data/Mask.nc', 'Mask'); % needed to store values for Nlat and Nlon 
Nlat=size(mask,1);
Nlon=size(mask,2);
Ntot=numel(mask);

land=mask<1.5; % logical vector storing land positions
X1=dlmread('Data/Dataset.txt', '');
w=X1(1,:);
w=w/sum(w);
X=X1(2:end,:); % rows of X contain the original shapes, as long 1D vectors
n=size(X,1);  % n is the number of original shapes
clear X1;

% p is here 17540 (only land cells considered)
Mn = mean(X,1); % 1xp vector
Xbar = X - ones(n,1)*Mn; % nxp

%Now compute Y=X*(sqrt(W))
Y=Xbar;
for j=1:size(Xbar,2)
    Y(:,j)=Xbar(:,j)*sqrt(w(j));
end

[~,S,PC_old] = svd(Y, 'econ');

PC_mid=PC_old; % PC_mid = W^(-1/2)*PC_old
for i=1:size(PC_mid,1)
    PC_mid(i,:)=PC_old(i,:)/sqrt(w(i));
end
clear PC_old

PC = zeros(Ntot,n-1);
PC(land,:)=PC_mid(:,1:n-1); % the last PC is meaningless and corresponds to zero eigenvalue
PC = reshape(PC, [Nlat, Nlon, n-1]);

M  = zeros(Nlat, Nlon);
M(land)=Mn;

Std = diag(S/sqrt(n-1));
Std = Std(1:end-1);

% This is in order to adjust the mask by leaving current values
PC=remask(PC, nan, 0);
M=remask(M, nan, 0);

end
