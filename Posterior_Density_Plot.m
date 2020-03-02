%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This program is structures into two main parts. In the first one, all   %%%
%%% the subplots needed in the final main figure are produced, and saved as %%%
%%% .fig images. In the second part, these images are loaded within a for   %%%
%%% loop, where they are inserted as subpanels of the main figure; on the   %%% 
%%% diagonal, the PCs are drawn.                                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% GENERAL VARIABLES

%{
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
MaxVal(c)=500*floor(MaxVal(c)/500);
end
%}


%% Color, Line and other settings

granata=[0.64 0.08 0.18];
gold=[0.9 0.75 0];
N=60;
str='../Pictures/NxN_Plots/';
%s=[str, 'Different Colours/', func2str(col_map_fun)];
str = '/home/dario/Desktop/JGR - Earth Surface paper/Plots/3x3_Plots/Posterior_Densities/';

line_col_1='r'; text_col_1='r';
line_col_2='b';  text_col_2='b';
dot_sz = 35;
ln_width =1.5;
ln_style_1 = '-.';
ln_style_2 = '-.';
font_sz = 14;
font = 'Century Schoolbook L';

col_map_fun=@gray;
cm=col_map_fun();
cm=flipud(col_map_fun());
cm=cm(2:end,:);

%{
%% Figure Cx=2, Cy=1
v1=[1,4]/100;
v2=[2,5,8]/100;

figure;
Cx=2; Cy=1;
xbins=linspace(-MaxVal(Cx), MaxVal(Cx), N+1);
ybins=linspace(-MaxVal(Cy), MaxVal(Cy), N+1);
Count_tot = histcounts2(Z(:,Cx), Z(:,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC1  =  histcounts2(Z(tight,Cx), Z(tight,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC3  =  histcounts2(Z(loose,Cx), Z(loose,Cy), xbins, ybins, 'Normalization', 'count');
ind=Count_tot>0;
Ratio1 = zeros(size(Count_tot));
Ratio1(ind) = Count_DC1(ind)./Count_tot(ind);
Ratio1 = Gaussian_Filtering(Ratio1', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);
Ratio3 = zeros(size(Count_tot));
Ratio3(ind) = Count_DC3(ind)./Count_tot(ind);
Ratio3 = Gaussian_Filtering(Ratio3', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);

%Plots
pcolor(xbins, ybins, add_nan(Count_tot')); %Background
axis square
colormap(cm); shading flat; hold on;
[~,h]=contour(xbins, ybins, add_nan(Ratio1), v1); %, 'ShowText', 'on');
set(h, 'LineColor', line_col_1, 'LineStyle', ln_style_1, 'LineWidth', ln_width);
[~,h]=contour(xbins, ybins, add_nan(Ratio3), v2); % 'ShowText', 'on'); 
set(h, 'LineColor', line_col_2, 'LineStyle', ln_style_2, 'LineWidth', ln_width);
scatter(335, 165, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_1);
scatter(240, 100, dot_sz,  'o', 'filled', 'MarkerFaceColor', line_col_2);

%text(620, -200, '\fontsize{11}1', 'Color', text_col_1);
t=text(630, -210, '1', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
%t=text(470, -150, '3', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(350,  -90, '4', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(-270, -620, '2', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text(-30, -370, '5', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text(-240, -180, '8', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;

savefig([str, num2str(Cx), '_' num2str(Cy)]);
close
%}

%{
%% Figure Cx=1, Cy=2
v1=[2,4]/100;
v2=[3,6,9]/100;

figure;
Cx=1; Cy=2;
xbins=linspace(-MaxVal(Cx), MaxVal(Cx), N+1);
ybins=linspace(-MaxVal(Cy), MaxVal(Cy), N+1);
Count_tot = histcounts2(Z(:,Cx), Z(:,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC1  =  histcounts2(Z(tight,Cx), Z(tight,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC3  =  histcounts2(Z(loose,Cx), Z(loose,Cy), xbins, ybins, 'Normalization', 'count');
ind=Count_tot>0;
Ratio1 = zeros(size(Count_tot));
Ratio1(ind) = Count_DC1(ind)./Count_tot(ind);
Ratio1 = Gaussian_Filtering(Ratio1', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);
Ratio3 = zeros(size(Count_tot));
Ratio3(ind) = Count_DC3(ind)./Count_tot(ind);
Ratio3 = Gaussian_Filtering(Ratio3', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);

pcolor(xbins, ybins, add_nan(Count_tot'));
axis square
colormap(cm); shading flat; hold on;

[~,h]=contour(xbins, ybins, add_nan(Ratio1), v1); %, 'ShowText', 'on');
set(h, 'LineColor', line_col_1, 'LineStyle', ln_style_1, 'LineWidth', ln_width);
[~,h]=contour(xbins, ybins, add_nan(Ratio3), v2);  % 'ShowText', 'on'); 
set(h, 'LineColor', line_col_2, 'LineStyle', ln_style_2, 'LineWidth', ln_width);
scatter(165, 335, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_1);
scatter(100, 240, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_2);

t=text(-100, 620, '2', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(-110, 350, '4', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(-580, -230, '3', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text(-490, 0, '6', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text(-200, 200, '9', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;

savefig([str, num2str(Cx), '_' num2str(Cy)]);
close
%}

%{
%% Figure Cx=3, Cy=1
v1=[1,2]/100;
v2=[2,6,10]/100;

figure;
Cx=3; Cy=1;
xbins=linspace(-MaxVal(Cx), MaxVal(Cx), N+1);
ybins=linspace(-MaxVal(Cy), MaxVal(Cy), N+1);
Count_tot = histcounts2(Z(:,Cx), Z(:,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC1  =  histcounts2(Z(tight,Cx), Z(tight,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC3  =  histcounts2(Z(loose,Cx), Z(loose,Cy), xbins, ybins, 'Normalization', 'count');
ind=Count_tot>0;
Ratio1 = zeros(size(Count_tot));
Ratio1(ind) = Count_DC1(ind)./Count_tot(ind);
Ratio1 = Gaussian_Filtering(Ratio1', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);
Ratio3 = zeros(size(Count_tot));
Ratio3(ind) = Count_DC3(ind)./Count_tot(ind);
Ratio3 = Gaussian_Filtering(Ratio3', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);

pcolor(xbins, ybins, add_nan(Count_tot'));
axis square
colormap(cm); shading flat; hold on;

[~,h]=contour(xbins, ybins, add_nan(Ratio1), v1); %'ShowText', 'on');
set(h, 'LineColor', line_col_1, 'LineStyle', ln_style_1, 'LineWidth', ln_width);
[~,h]=contour(xbins, ybins, add_nan(Ratio3), v2);%'ShowText', 'on'); 
set(h, 'LineColor', line_col_2, 'LineStyle', ln_style_2, 'LineWidth', ln_width);
scatter(190, 235, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_1);
scatter(145, 40, dot_sz, 'o',  'filled', 'MarkerFaceColor', line_col_2);

t=text( 260, 480, '1', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text( 220, 220, '2', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(-340, -450, '2', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text(-200, -420, '6', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text( -10, -260, '10', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;

savefig([str, num2str(Cx), '_' num2str(Cy)]);
close

%}

%{
%% Figure Cx=1, Cy=3
v1=[1,2]/100;
v2=[4,8,12]/100;

figure;
Cx=1; Cy=3;
xbins=linspace(-MaxVal(Cx), MaxVal(Cx), N+1);
ybins=linspace(-MaxVal(Cy), MaxVal(Cy), N+1);
Count_tot = histcounts2(Z(:,Cx), Z(:,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC1  =  histcounts2(Z(tight,Cx), Z(tight,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC3  =  histcounts2(Z(loose,Cx), Z(loose,Cy), xbins, ybins, 'Normalization', 'count');
ind=Count_tot>0;
Ratio1 = zeros(size(Count_tot));
Ratio1(ind) = Count_DC1(ind)./Count_tot(ind);
Ratio1 = Gaussian_Filtering(Ratio1', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);
Ratio3 = zeros(size(Count_tot));
Ratio3(ind) = Count_DC3(ind)./Count_tot(ind);
Ratio3 = Gaussian_Filtering(Ratio3', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);

pcolor(xbins, ybins, add_nan(Count_tot'));
axis square
colormap(cm); shading flat; hold on;

[~,h]=contour(xbins, ybins, add_nan(Ratio1), v1); %, 'ShowText', 'on');
set(h, 'LineColor', line_col_1, 'LineStyle', ln_style_1, 'LineWidth', ln_width);
[~,h]=contour(xbins, ybins, add_nan(Ratio3), v2); %, 'ShowText', 'on'); 
set(h, 'LineColor', line_col_2, 'LineStyle', ln_style_2, 'LineWidth', ln_width);
scatter(235, 190, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_1);
scatter(40, 145, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_2);

t=text( 420, 280, '1', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text( 40, 240, '2', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(-630, -140, '4', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text( -320, -100, '8', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text( -230, 100, '12', 'Color', text_col_2); t.FontSize=0.9*font_sz; t.FontName=font;

savefig([str, num2str(Cx), '_' num2str(Cy)]);
close
%}

%{
%% Figure Cx=3, Cy=2
v1=[1,3]/100;
v2=[3,7,11]/100;

figure;
Cx=3; Cy=2;
xbins=linspace(-MaxVal(Cx), MaxVal(Cx), N+1);
ybins=linspace(-MaxVal(Cy), MaxVal(Cy), N+1);
Count_tot = histcounts2(Z(:,Cx), Z(:,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC1  =  histcounts2(Z(tight,Cx), Z(tight,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC3  =  histcounts2(Z(loose,Cx), Z(loose,Cy), xbins, ybins, 'Normalization', 'count');
ind=Count_tot>0;
Ratio1 = zeros(size(Count_tot));
Ratio1(ind) = Count_DC1(ind)./Count_tot(ind);
Ratio1 = Gaussian_Filtering(Ratio1', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);
Ratio3 = zeros(size(Count_tot));
Ratio3(ind) = Count_DC3(ind)./Count_tot(ind);
Ratio3 = Gaussian_Filtering(Ratio3', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);

pcolor(xbins, ybins, add_nan(Count_tot'));
axis square
colormap(cm); shading flat; hold on;

[~,h]=contour(xbins, ybins, add_nan(Ratio1), v1); %'ShowText', 'on');
set(h, 'LineColor', line_col_1, 'LineStyle', ln_style_1, 'LineWidth', ln_width);
[~,h]=contour(xbins, ybins, add_nan(Ratio3), v2); %'ShowText', 'on'); 
set(h, 'LineColor', line_col_2, 'LineStyle', ln_style_2, 'LineWidth', ln_width);
scatter(-35, 370, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_1);
scatter(185, -235, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_2);

t=text( -280, 680, '1', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text( -220, 450, '3', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text( -280, 100, '3', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text(  205, 70, '7', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text( 200, -160, '11', 'Color', text_col_2); t.FontSize=0.9*font_sz; t.FontName=font;

savefig([str, num2str(Cx), '_' num2str(Cy)]);
close
%}


%%{
%% Figure Cx=2, Cy=3
v1=[2,4]/100;
v2=[1,5,9]/100;

figure;
Cx=2; Cy=3;
xbins=linspace(-MaxVal(Cx), MaxVal(Cx), N+1);
ybins=linspace(-MaxVal(Cy), MaxVal(Cy), N+1);
Count_tot = histcounts2(Z(:,Cx), Z(:,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC1  =  histcounts2(Z(tight,Cx), Z(tight,Cy), xbins, ybins, 'Normalization', 'count');
Count_DC3  =  histcounts2(Z(loose,Cx), Z(loose,Cy), xbins, ybins, 'Normalization', 'count');
ind=Count_tot>0;
Ratio1 = zeros(size(Count_tot));
Ratio1(ind) = Count_DC1(ind)./Count_tot(ind);
Ratio1 = Gaussian_Filtering(Ratio1', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);
Ratio3 = zeros(size(Count_tot));
Ratio3(ind) = Count_DC3(ind)./Count_tot(ind);
Ratio3 = Gaussian_Filtering(Ratio3', linspace(-0.1, 0.1, N), linspace(-0.1, 0.1, N), 7, 7, 50);

pcolor(xbins, ybins, add_nan(Count_tot'));
axis square
colormap(cm); shading flat; hold on;

[~,h]=contour(xbins, ybins, add_nan(Ratio1), v1); %, 'ShowText', 'on');
set(h, 'LineColor', line_col_1, 'LineStyle', ln_style_1, 'LineWidth', ln_width);
[~,h]=contour(xbins, ybins, add_nan(Ratio3), v2); %, 'ShowText', 'on'); 
set(h, 'LineColor', line_col_2, 'LineStyle', ln_style_2, 'LineWidth', ln_width);
scatter(370, -35, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_1);
scatter(-235, 185, dot_sz, 'o', 'filled', 'MarkerFaceColor', line_col_2);

t=text(  500, -260, '2', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text(  350, -145, '4', 'Color', text_col_1); t.FontSize=font_sz; t.FontName=font;
t=text( -580, 320, '1', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;
t=text( -180, 335, '5', 'Color', text_col_2); t.FontSize=0.95*font_sz; t.FontName=font;
t=text( -10, 195, '9', 'Color', text_col_2); t.FontSize=font_sz; t.FontName=font;

savefig([str, num2str(Cx), '_' num2str(Cy)]);
close
%}


%% FOR loop to produce final picture

set(groot, 'defaultAxesTickLabelInterpreter','latex');

fact_space=1/25; % the bigger, the more space between plots
scale = 1.2;     % scale factor for size of panels
transl = [-0.04, -0.03]; % [a,b]: moves each panel of a units to the right, 
                         % and b units to the top
font_sz = 18;   % font size for axis numbers in panels
tick_lngt=0.023;

hf=figure('units','normalized','outerposition',[0 0 1 1]);
x0=20; y0=20; 
width=770; height=750;
width=800; height=750;
set(gcf,'units','points','position',[x0,y0,width,height]);
        
for Cy=1:n  % row index, we consider starting from below
    %    subaxis(r,r, sub2ind([r,r], i, r+1-i), 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
    H=PC(:,:,Cy);
    index_nan=isnan(H); H(index_nan)=0;
    H=smooth_Greenland(H, 5, 9, 100); H(index_nan)=nan;
    p=subplot(n,n, sub2ind([n,n], Cy, n+1-Cy)); pos=get(p, 'Position');
    pos(1:2) = pos(1:2) + (2*pos(1:2)-1)*fact_space; % enlargen around center of picture (rescale around (0.5, 0.5))
    pos(1:2) = pos(1:2) + transl; % translate everything towards left and bottom
    pos(3:4) = scale*pos(3:4);    % enlargen width and height of each panel
    set(p, 'Position', pos);
    M = max( abs(nanmin(H(:))), abs(nanmax(H(:))) );
    ax=worldmap([59.5,84.1], [-74,-10.5]);
    setm(ax, 'grid', 'off', 'frame', 'off', 'meridianlabel', 'off', 'parallellabel', 'off');
    pcolorm(lat2, lon2, add_nan(H));
    plotm(coastlat(Greenland), coastlon(Greenland), 'k', 'LineWidth', 1.5);
    caxis([-M,M]);
    colormap(p,redblue); %ch=colorbar;
    %set(ch, 'TickLength', 0.025, 'TickDirection', 'both', 'LineWidth', 1.5);
    %set(ch, 'FontSize', 20); ch.Location='westoutside';
    th=title(['PC ', num2str(Cy)], 'FontSize', 20);
    pos=get(ax, 'Position'); 
    pos(3)=1.1*pos(3); pos(1)=pos(1)-0.01; set(ax, 'Position', pos);
    
    for Cx=Cy+1:n
        pause(0.2);
        p1=subplot(n,n, sub2ind([n,n], Cx, n+1-Cy));  % bottom right half
        axis square;
        pos=get(p1, 'Position');
        pos(1:2) = pos(1:2) + (2*pos(1:2)-1)*fact_space; % enlargen around center of picture (rescale around (0.5, 0.5))
        pos(1:2) = pos(1:2) + transl; % translate everything towards left and bottom
        pos(3:4) = scale*pos(3:4);        % enlargen width and height of each panel
        set(p1, 'Position', pos);
        h1 = openfig([str, num2str(Cx), '_' num2str(Cy) '.fig'],  'reuse');
        ax1 = gca; fig1 = get(ax1,'children');
        copyobj(fig1,p1, 'legacy');
        close(h1);
        box on; colormap(p1,cm);
        set(gca, 'XTickLabel', [], 'YTickLabel', [], 'LineWidth', 1.5, 'TickLength', 1.4*get(gca,'TickLength'));
        set(gca, 'TickDir', 'both'); set(gca, 'FontSize', font_sz);
        Mx=max(abs(xticks)); p1.XAxis.MinorTick='on';
        XVal=linspace(-Mx,Mx,5);
        p1.XAxis.TickLength(1)=tick_lngt; p1.XAxis.TickValues=XVal;
        p1.XAxis.MinorTickValues=linspace(-Mx,Mx,2*length(XVal)-1);
        My=max(abs(yticks)); p1.YAxis.MinorTick='on';
        YVal=linspace(-My,My,5);
        p1.YAxis.TickLength(1)=tick_lngt; p1.YAxis.TickValues=YVal;
        p1.YAxis.MinorTickValues=linspace(-My,My,2*length(YVal)-1);
        if mod(n,2)==0
            if Cy==1 && mod(Cx,2)==0 % bottom row, even cells
                set(gca, 'XTickLabel', num2cell(XVal));
            end
            if Cx==n && mod(Cy,2)==1 % last column, odd cells (starting counting from below)
                set(gca, 'YAxisLocation', 'right');
                set(gca, 'YTickLabel', num2cell(YVal));
            end
            
        else % n is odd
            if Cy==1 && mod(Cx,2)==1 % bottom row, odd cells
                set(gca, 'XTickLabel', num2cell(XVal));
            end
            if Cx==n && (mod(Cy,2)==0 || Cy==1)
                set(gca, 'YAxisLocation', 'right');
                set(gca, 'YTickLabel', num2cell(YVal));
            end
        end
        
        pause(0.2);
        p2=subplot(n,n, sub2ind([n,n], Cy, n+1-Cx)); % top left part
        axis square;
        pos=get(p2, 'Position');
        pos(1:2) = pos(1:2) + (2*pos(1:2)-1)*fact_space; % enlargen around center of picture (rescale around (0.5, 0.5))
        pos(1:2) = pos(1:2) + transl; % translate everything towards left and bottom
        pos(3:4) = scale*pos(3:4);    % enlargen width and height of each panel
        set(p2, 'Position', pos);
        h2 = openfig([str, num2str(Cy), '_' num2str(Cx) '.fig'],  'reuse');
        ax2 = gca; fig2 = get(ax2,'children');
        copyobj(fig2,p2, 'legacy');
        close(h2);
        box on; colormap(p2, cm);
        set(gca, 'XTickLabel', [], 'YTickLabel', [], 'LineWidth', 1.5, 'TickLength', 1.4*get(gca,'TickLength'));
        set(gca, 'TickDir', 'both'); set(gca, 'FontSize', font_sz);
        Mx=max(abs(xticks)); p2.XAxis.MinorTick='on'; 
        XVal=linspace(-Mx,Mx,5);
        p2.XAxis.TickLength(1)=tick_lngt; p2.XAxis.TickValues=XVal;
        p2.XAxis.MinorTickValues=linspace(-Mx,Mx,2*length(XVal)-1);
        My=max(abs(yticks)); p2.YAxis.MinorTick='on'; 
        YVal=linspace(-My,My,5);
        p2.YAxis.TickLength(1)=tick_lngt; p2.YAxis.TickValues=YVal;
        p2.YAxis.MinorTickValues=linspace(-My,My,2*length(YVal)-1);       
        if mod(n,2)==0
            if Cy==1 && mod(Cx,2)==0 % left column, even cells
                set(gca, 'YTickLabel', num2cell(YVal));
            end
            if Cx==n && mod(Cy,2)==1 % j==t: top row, odd cells (starting counting from below)
                set(gca, 'XAxisLocation', 'top');
                set(gca, 'XTickLabel', num2cell(XVal));
            end
        else % n is odd
            if Cy==1 && mod(Cx,2)==1 % left column
                set(gca, 'YTickLabel', num2cell(YVal));
            end
            if Cx==n &&(mod(Cy,2)==0 || Cy==1)    % j==t: top row
                set(gca, 'XAxisLocation', 'top');
                set(gca, 'XTickLabel', num2cell(XVal));
            end
        end % end of bigger if...else.
        
    end % end for in j
end


print([str, '3x3_gray'], '-dpng'); close(hf);
%print([str, '3x3_gray'], '-dsvg'); close(hf);
