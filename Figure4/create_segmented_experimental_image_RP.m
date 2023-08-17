close all; clear all; clc

%% Preprocessing Experimental data
image = "Experimental_data_RP.png";

[I] = imread(image);
greyI = rgb2gray(I);
% rotate image to correct 
greyI = imrotate(greyI,-1,'bilinear','crop');
BW = imbinarize(greyI);
[B,image_mask] =bwboundaries(BW,8);
thresh_im=(image_mask<1);
[bdry] = edge(thresh_im,'Sobel');

%%get the center row
ind = find(BW ==1);
[row,col] = ind2sub(size(bdry),ind);
%find the midpoint of the first col
first_col = find(col==1);
first_row = row(first_col);
first_mid = min(first_row)+round((max(first_row)-min(first_row))/2);
%find the midpoint of the last col
last_col = find(col==size(bdry,2));
last_row = row(last_col);
last_mid = min(last_row)+round((max(last_row)-min(last_row))/2);

% get the fiber thickness
fiber_image = edge(BW,'Canny');
ind_fiber = find(fiber_image ==1);
[row,col] = ind2sub(size(bdry),ind_fiber);
selected_col = 32; %for default
thresh = multithresh(row(find(col==selected_col)),6);
top_fiber_bdry =  round(thresh(3));
fiber_thickness = (first_mid-top_fiber_bdry);
figure;imshow(bdry); 
hold on;
plot([0 size(bdry,2)],[first_mid+fiber_thickness  last_mid+fiber_thickness],'b','LineWidth',1); hold on;
plot([0 size(bdry,2)],[first_mid-fiber_thickness  last_mid-fiber_thickness],'b','LineWidth',1); hold on;

