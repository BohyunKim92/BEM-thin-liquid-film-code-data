close all; clear all; clc
image = "Experimental_data_IS.jpg";
[I] = imread(image);
greyI = im2gray(I);
greyI = imrotate(greyI,-0.1,'bilinear','crop');
BW = imbinarize(greyI);
[B,image_mask] =bwboundaries(BW,8);
thresh_im=(image_mask<1);
[bdry] = edge(thresh_im,'Sobel');
%%get the medium row
ind = find(bdry ==1);
[row,col] = ind2sub(size(bdry),ind);
%find the midpoint of the first col
first_col = find(col==860);
first_row = row(first_col);
first_mid = 113;
%find the midpoint of the last col
last_col = find(col==size(bdry,2));
last_row = row(last_col);
last_mid = first_mid;


%% get the fiber thickness
fiber_image = edge(BW,'Canny');
ind_fiber = find(fiber_image ==1);
[row,col] = ind2sub(size(bdry),ind_fiber);
selected_col = 692; 
thresh = multithresh(row(find(col==selected_col)),6);
top_fiber_bdry =  110; 
fiber_thickness = (first_mid-top_fiber_bdry);

%cropping from 603- 2006
cropped_bdry = bdry(:,603:2006);

%figure; 
imshow(cropped_bdry);
hold on;
plot([0 size(bdry,2)],[first_mid+fiber_thickness  last_mid+fiber_thickness],'b','LineWidth',1); hold on;
plot([0 size(bdry,2)],[first_mid-fiber_thickness  last_mid-fiber_thickness],'b','LineWidth',1); hold on;



