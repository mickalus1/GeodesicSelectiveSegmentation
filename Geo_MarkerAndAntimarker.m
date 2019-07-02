function [ SSF, cols, rows, R_min,R_max,C_min,C_max,z] = Geo_MarkerAndAntimarker( z, cols, rows,...
                                                    cols_anti, rows_anti )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

global use_antimarkers edit_image use_cropping

if nargin == 1
figure('units','normalized','outerposition',[0 0 1 1]);...
    imagesc(z);colormap(gray) % Me
title('Marker Set')
[~,cols,rows] = freehanddraw(gca,'color','g','linewidth',3);
Mask = roipoly(z,cols,rows);
close
if use_antimarkers == 1
figure('units','normalized','outerposition',[0 0 1 1]);...
    imagesc(z);colormap(gray);hold on;contour(Mask,[0.5 0.5],'g','LineWidth',3);
title('Anti-Marker Set')
[~,cols_anti,rows_anti] = freehanddraw(gca,'color','r','linewidth',3);
close
Mask_Anti = roipoly(z,cols_anti,rows_anti);
end

else
    Mask = roipoly(z,cols,rows);
end

if use_cropping
    boundary = 40;
    C_row = round(mean(rows));
    C_col = round(mean(cols));
    R_min = max(C_row - boundary,1);
    R_max = min(C_row + boundary,size(Mask,1));
    C_min = max(C_col - boundary,1);
    C_max = min(C_col + boundary,size(Mask,2));
    Mask = Mask(R_min:R_max,C_min:C_max);
    if use_antimarkers
    Mask_Anti = Mask_Anti(R_min:R_max,C_min:C_max);
    end
    z = z(R_min:R_max,C_min:C_max);
end

[rows2,cols2] = find(Mask==1);
rows = rows2;
cols = cols2;

if (size(cols,1)*size(cols,2) == 1) &&... 
        (size(rows,1)*size(rows,2) == 1)
       Mask(rows:rows+1,cols:cols+1) = 1;
end


% load('HeartCT_Seg_Markers')

%% Anti-Markers
if use_antimarkers == 1
    [rows2_anti,cols2_anti] = find(Mask_Anti==1);
rows_anti = rows2_anti;
cols_anti = cols2_anti;
Pd_anti = geodesic_map(z,(Mask_Anti>0.5)-0.5, rows_anti, cols_anti);
% SSF_anti = 1 - Pd_anti/max(Pd_anti(:));
SSF_anti2 = Pd_anti/max(Pd_anti(:));
alpha = 200;
SSF_anti = (exp(-alpha*SSF_anti2)-exp(-alpha))/(1-exp(-alpha));
end

z_new = z;

if use_antimarkers ==1
BW = roipoly(z,round(cols_anti),round(rows_anti));
E = edge(BW);
se = strel('disk',2);
E = imdilate(E,se);
[rows_anti2,cols_anti2] = find(E>0.5);
for i = 1:size(rows_anti2,1)
z_new(round(rows_anti2(i)),round(cols_anti2(i))) = 1-z(round(rows_anti2(i)),round(cols_anti2(i)));
end
end
if edit_image == 0
z_new = z;
end

%% Markers
Pd1 = geodesic_map(z_new,(Mask>0.5)-0.5, rows, cols);
SSF_markers = Pd1/max(Pd1(:));

%% Combine

if use_antimarkers ==1
    if edit_image ==0
        SSF = SSF_markers + SSF_anti;
    else
        SSF = SSF_markers;
    end
else
 SSF = SSF_markers;% + SSF_anti;   
end

SSF = (SSF - min(SSF(:)))/(max(SSF(:)) - min(SSF(:)));

end

