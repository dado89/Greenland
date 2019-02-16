% Given PCs and coefficients, the function builds N new shapes as linear
% combinations of the PCs, and then adjusts result according to bedrock.

% INPUTS
% PC:    (Nlat)x(Nlon)x n matrix with n Principal Components
% M:     (Nlat)x(Nlon) matrix of mean to add to linear combination
% coeff: Nxr matrix of coefficients. In each row, the r coefficients
%        corresponding to the PCs for that input are specified.

% OUTPUT
% Z: (Nlat)x(Nlon)x N, with j-th new shape in Z(:,:,j).


function Z = build_shapes(PC, M, coeff)

Nlat=size(PC,1);
Nlon=size(PC,2);
p=Nlat*Nlon;

N = size(coeff,1);
r = size(coeff,2);
if r>size(PC,3)
    error('The coefficients are provided for a number of PCs greater than the PCs available');
end

% Retain only needed components and reshape all quantities as vectors to 
% allow multiplications
PC=PC(:, :, 1:r);
PC=reshape(PC, [p,r]);
M=reshape(M, [p,1]);

% Generate shapes as linear combination, and then adjust bedrock
Z = M*ones(1,N) +  PC*(coeff'); % p x N
Z = reshape(Z, [Nlat, Nlon,N]);

str='Original Morphologies/Regridded Morphologies/nc files/Stone_123.5_Regrid.nc';
Bed=ncread(str, 'Bedrock');
%Bed=remask(Bed, nan, 0);  % puts nan outside physical_mask, and 0s in cells being inside the
                          % physical mask but outside the original "larger" one
Z(Z<0)=0;
Z(isnan(Z))=0;
for j=1:N
    Z(:,:,j) = max(Z(:,:,j), Bed);
end

end
