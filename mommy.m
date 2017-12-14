function [res, initno]=mommy_myelin(filein,mag)
addpath('/Users/alex/Documents/MATLAB/emanalysis/')
%sample call
% [res, initno]=mommy_myelin('/Users/alex/Desktop/andy/ematpathlab/CCALOSSUM/CVN/130828_3_CC_cvn/7100/130828_3_CC_cvn_030.TIF', 7100)

if mag==2650
    radius_mf=3;
    myarea=2500;
else
    radius_mf=6;
    myarea=4000;
end
    
im1=imread(filein); 
%'Axon1.TIF'
im1=im1(1:2624,1:2624);

figure(1)
subplot(1,3,1)
imagesc(im1);
title('original')

flag =1;

switch 1
    case flag==1
        %enhance contrast
        %flag=1
        im11=imadjust(im1,[0.3 0.7],[]);
        
        
        im3 = medfilt2(im11,[radius_mf radius_mf]);
        subplot(1,3,2)
        imagesc(im3);
        title('median filtered')
        
        level = graythresh(im3);
        BW = im2bw(im3,0.9*level);
        BW=uint8(BW);
        subplot(1,3,3)
        imagesc(BW);
        title('otsu thresholding')
        
        
        
    case flag==2
        %average filter
        %flag=2
        im2 = filter2(fspecial('average',3),im1)/255;
        subplot(2,2,3)
        imagesc(im2);
        BW = im2bw(im2,level);
        BW=uint8(BW);
        subplot(2,2,3)
        imagesc(BW);
        
    case flag==3
        %flag=3
        im3 = medfilt2(im2,[3 3]);
        subplot(2,2,4)
        imagesc(im3);
        BW = im2bw(im3,level);
        BW=uint8(BW);
        subplot(2,2,4)
        imagesc(BW);
    case flag==4
        %flag=4;
        nhoodSize = 3;
        smoothValue  = 0.01*diff(getrangefromclass(im3)).^2;%was 0.001
        B = imguidedfilter(im1, im3, 'NeighborhoodSize',nhoodSize, 'DegreeOfSmoothing',smoothValue);
        figure(2)
        subplot(1,2,1)
        imagesc(B);
        sqrt(500
        BW = im2bw(B,level);
        BW=uint8(BW);
        subplot(1,2,2)
        imagesc(BW);
end

colormap(gray)

% colormap(gray)
% figure(1)
% 
% subplot(2,2,1)
% imagesc(im1)
% title('original');
% 
% subplot(2,2,2)
% imagesc(BW)
% title('thresh 1');





figure(3)
title('erode')
mydisk=strel('disk',2); %mydisk=strel('disk',4);
% im_dil=imdilate(imdilate(imdilate(BW,mydisk),mydisk),mydisk);
% im_erode=imerode(imerode(imerode(im_dil,mydisk),mydisk),mydisk);

if mag==2650
    imdil=BW;
else
im_dil=imdilate(imdilate(BW,mydisk),mydisk);
end

im_erode=imerode(imerode(im_dil,mydisk),mydisk);


subplot(1,2,1);
imagesc(im_erode);
title('dilx3,errx3')

cc=bwconncomp(~im_erode);
stats_erode=regionprops(cc,'Area');
ind1=find([stats_erode.Area] > myarea);%was 2500
bw2=ismember(labelmatrix(cc),ind1);
imshow(bw2)

%serode  = regionprops(~im_erode, 'centroid');
serode  = regionprops(~bw2, 'centroid');
centroids_erode = cat(1, serode.Centroid);
  
   hold on
   plot(centroids_erode(:,1), centroids_erode(:,2), 'r*')
   hold off
   title('eroded')
   

  initno=size(centroids_erode);
  initno=initno(1);
  
  
  
im_open=imopen(BW,mydisk);
subplot(1,2,2);
imagesc(im_open);
title('open');



cc=bwconncomp(~im_open);
stats_open=regionprops(cc,'Area');
ind1=find([stats_open.Area] > myarea);
bw2=ismember(labelmatrix(cc),ind1);
imshow(bw2)
title('open');

res=bw2;

sopen  = regionprops(~bw2, 'centroid');
centroids_open = cat(1, sopen.Centroid);
  

   hold on
   plot(centroids_open(:,1), centroids_open(:,2), 'r*')
   hold off
   
initno2=size(centroids_open);
  initno2=initno2(1)
  
%keep circular areas only

%make sure center is on black , not white

%active snakes?
 figure(5)
   [D,IDX] = bwdist(~bw2,'euclidean');
   
    mymap=jet(24);
   mymap(1,:)=[0 0 0];
   
    %IMGRADIENT
    
  % BW = imregionalmax(D,4);
   D2=D*0;
   D2(find(D))=D(find(D));
   [counts,x] = imhist(D2,100); stem(x,counts)
   imagesc(D2)
   colormap(mymap);
   %ultimateErosion = bwulterode(originalBW);-use on myelin
   