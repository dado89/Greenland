%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program is structures into two main parts. In the first one, all   %%%
%%% the subplots needed in the final main figure are produced, and saved as %%%
%%% .fig images. In the second part, these images are loaded within a for   %%%
%%% loop, where they are inserted as subpanels of the main figure; on the   %%% 
%%% diagonal, the PCs are drawn.                                            %%%
%%% The script is particularly long, since it is in control of a high       %%%
%%% number of graphical details of the final picture.                       %%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GENERAL VARIABLES

load('Inputs.mat');
Z=Input_par; clear Input_par;
load('Compatibility_Variables_tight', 'fully_comp_tight');
load('Compatibility_Variables_middle', 'fully_comp_middle');
load('Compatibility_Variables_loose', 'fully_comp_loose');

tight = fully_comp_tight; clear fully_comp_tight
middle = fully_comp_middle; clear fully_comp_middle
loose = fully_comp_loose; clear fully_comp_loose

lat = ncread('Data/Mask.nc', 'lat');
lon = ncread('Data/Mask.nc', 'lon');
[lat2, lon2]=adjust_latlon(lat,lon);
load coastlines.mat
Greenland=3977:4210;

[PC, ~, ~] = pca_greenland;

n=3;
ind=1:30000;
MaxVal=zeros(n,1);
for c=1:n
MaxVal(c)=max(abs(Z(ind,c)));
