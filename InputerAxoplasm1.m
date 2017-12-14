clear all
train=0;%1 to generate new inputs for training NN
pos = 1;
dog=0;
prompt = 'Axons:';

mag=2650;
filein='/Users/alex/Desktop/andy/train_140228_1_NOS_2650.TIF'; %'/Users/alex/Documents/MATLAB/Axon1.TIF';
 
mag=7100;
%nos1
runno='140313';
index={'08', '09', '10', '11', '12', '18', '19', '20', '22', '23' , '27', '28', '35', '39', '40', '41'}
filein='/Users/alex/Desktop/andy/ematpathlab/FORNIX/NOS/140313_1_NOS_f/7100/140313_NOS_008.TIF'; %'/Users/alex/Desktop/andy/train_140228_1_NOS_7100.TIF'; %
%nos2
% runno='140228_1'
% index={ '09', '10', '11', '18', '19', '20', '21', '22'};
% filein=['/Users/alex/Desktop/andy/ematpathlab/FORNIX/NOS/' runno '_NOS_f/7100/' runno '_NOS_009.TIF'];

%cvn1
% runno='130828_3';
% index={'05' '06' '18' '19' '21' '24' '28' '29' '31' '35' '40' '41' '42' '48' '49' '50' '54' '65' '66' '67' '68'}; 
% filein=['/Users/alex/Desktop/andy/ematpathlab/FORNIX/CVN/' runno '_CVN_f/7100/' runno '_CVN_005.TIF'];
% 
% runno='130930_2';
% index={'06' '15' '16' '17' '18' '19' '20' '21' '24' '36' '37' '38' '41' '42' '43' '44' '45' '46' '47' '48' '59' '60' '61' '62' '63'}; 

%/Users/alex/Desktop/andy/ematpathlab/FORNIX/CVN/130930_2_CVN_f/7100/130930_2cvn_006.TIF





root='FORNIX/NOS';
%root='FORNIX/CVN'
fiber='FORNIX';
genotype='NOS'; %'CVN';
genotype2='nos'; %'cvn';

pathin=['/Users/alex/Desktop/andy/ematpathlab/' fiber '/' genotype '/' runno '_' genotype '_f/7100/'];
flist=dir(fullfile(pathin, '*TIF'))

nn=numel(flist);

for imindex=1:nn %numel(index)
    imindex

  

filein=[pathin flist(imindex).name];

rawimg = imread(filein); 
%prompt = 'Axons:';
% rawimg = imgaussfilt(rawimg,2);
% [lowthresh, c]=graythresh(rawimg);
% lowthresh=lowthresh*255;
% %img = rawimg > 0.95*lowthresh;
% img = rawimg > lowthresh*0.90;
% ker = ones(1);
% for i = 1:20
%     img = imerode(img,ker);
%     img = imdilate(img,ker);
% end
% img = imerode(img,ker);
close all
[img, rawimg_small,init_no] = mommy_axoplasm_neg(filein,mag);
 

%img = imcomplement(img);
props = regionprops(img,rawimg_small,'Area','Eccentricity','EquivDiameter','EulerNumber',...
    'Extent','FilledArea','MajorAxisLength','MinorAxisLength',...
    'Orientation','Perimeter','Solidity','MeanIntensity', 'ConvexArea');



showaxons = regionprops(img, 'centroid');%was ~bw2
centroids_axons = cat(1, showaxons.Centroid);
  
figure(1)
imagesc(rawimg)
hold on
plot(centroids_axons(:,1), centroids_axons(:,2), 'r*')
hold off
title('innitial centroids')

%generate secondary features for NN

clearvars FilledPercent HoleArea Geomean PerArea AxesRatio Convexity

for ii=1:(size(props,1))
    FilledPercent(ii)= props(ii).Area./props(ii).FilledArea;

    HoleArea(ii)= props(ii).FilledArea-props(ii).Area;

    Geomean(ii)= geomean([props(ii).MajorAxisLength,props(ii).MinorAxisLength]);

    PerArea(ii)= (props(ii).Perimeter./props(ii).FilledArea);
    
    AxesRatio(ii)=(props(ii).MajorAxisLength./props(ii).MinorAxisLength);
    
    Convexity(ii)=props(ii).Area./props(ii).ConvexArea;
    
    Elongation(ii)=props(ii).EquivDiameter./props(ii).MinorAxisLength
    
    
end



%saved as magAxoplasmSet13
% Inputs1 = [props.Eccentricity; ...
%     props.EulerNumber; ...
%     props.Solidity; ...
%     props.Extent; ...
%     props.Orientation;...
%     props.Perimeter;...
%     FilledPercent;...
%     HoleArea;...
%     Geomean;...
%     AxesRatio;...
%     Convexity;...
%     Elongation;...
%     PerArea]';

%saved as magAxoplasm
Inputs1 = [props.Area; props.MajorAxisLength; props.MinorAxisLength; props.Eccentricity;...
    props.Orientation; props.FilledArea;props.EulerNumber; props.EquivDiameter;...
    props.Solidity; props.Extent; props.Perimeter;FilledPercent;...
    HoleArea;Geomean;PerArea].';


%get classified image from NN
Axons=round(abs(myNNAxoplasm1(Inputs1')),1);

ind_axons=find(Axons(1,:)>=0.5);
ind_spaces=find(Axons(2,:)>=0.5);

%andy
Asum = sum(Axons);

%alex
Asum=(sum(Axons,2));
Asum=Asum(1)


figure(2)
imagesc(rawimg_small)
colormap gray
hold on
plot(centroids_axons(ind_axons,1), centroids_axons(ind_axons,2), 'g*')
hold off
title('final centroids after NN');
saveas(gcf, ['/Users/alex/Desktop/andy/ematpathlab/stats_em/' 'fig_testcentroids' num2str(imindex) '.png']);
%% 

%make function to return g ratio
mybwimage=img;
mycentroids=centroids_axons(ind_axons,:);
mydiameters=cat(1, props(ind_axons).EquivDiameter);
myorientation=cat(1,props(ind_axons).Orientation);
mygratios=gratio_noshow(mybwimage, rawimg_small,mycentroids, mydiameters, myorientation);

%% make function to write stats

res=write_emstats(filein,root, imindex, Asum,mygratios, mydiameters); 
%% 


j=1;
g=0;
if dog ~=0
for ii=1:size(props,1)
    HoleArea(ii)= (props(ii).Area-props(ii).FilledArea)*(-1);
end
for i = 1:size(Axons,2)
    if HoleArea(i) > 1000
        if Axons(i) > 1
            if Axons(i) < 4
                R1 = sqrt(HoleArea(i));%axoplasm
                R2 = sqrt(Inputs1(i,6));%fill area; myelin + axoplasm
                %Thickness = R1-R2;
                %Thinkness(j) = Thickness/R2;
                %The ratio of the inner axonal diameter to the total outer diameter or g-ratio i
              g(j)=R1/R2;
              j=j+1;
            end
        else
            'no axons found'
            g(j)=0;
        end
        
    end
end

g
end
    %figure(10);
    %imshow(rawimg);
    %hold
%plot(allres(1,:),allres(2,:),'+r');
%hold


end
%end loopfor imindex
%% 

%train the network
if train==1;
for i = 1:1
% clear FilledPercent
% clear Geomean
% clear HoleArea
% clear PerArea
% clear Inputs
% clear Inputs2
% clear props
% rawimg = imread('/Users/alex/Documents/MATLAB/Axon2.TIF'); %'TestImg_CHT_a2.bmp');
% rawimg = imgaussfilt(rawimg,1);
% [lowthresh, c]=graythresh(rawimg);
% lowthresh=lowthresh*255;
% %img = rawimg > 0.95*lowthresh;
% img2 = rawimg > lowthresh*0.93;
% ker = ones(1);
% for i = 1:10
%     img2 = imerode(img2,ker);
%     img2 = imdilate(img2,ker);
% end
% img2 = imcomplement(img2);
% props = regionprops(img2,'Area','Eccentricity','EquivDiameter','EulerNumber',...
%     'Extent','FilledArea','MajorAxisLength','MinorAxisLength',...
%     'Orientation','Perimeter','Solidity')
% 
% for ii=1:(size(props,1))
%     FilledPercent(ii)= props(ii).Area./props(ii).FilledArea;
% end
% for ii=1:size(props,1)
%     HoleArea(ii)= props(ii).Area-props(ii).FilledArea*-1;
% end
% for ii=1:size(props,1)
%     Geomean(ii)= geomean([props(ii).MajorAxisLength,props(ii).MinorAxisLength]);
% end
% for ii=1:size(props,1)
%     PerArea(ii)= (props(ii).Perimeter./props(ii).FilledArea);
% end
% 
% Inputs2 = [props.Area; props.MajorAxisLength; props.MinorAxisLength; props.Eccentricity;...
%     props.Orientation; props.FilledArea;props.EulerNumber; props.EquivDiameter;...
%     props.Solidity; props.Extent; props.Perimeter;FilledPercent;...
%     HoleArea;Geomean;PerArea].';
% Inputs = [Inputs1;Inputs2];
% Cents = regionprops(img,'Centroid');
% Cents = Cents(1:size(Axons));
% for i = 1:size(Axons)
% if abs(Axons(i)) > 0.3
% a(i)=(Cents(i))
% end
% end
% 
% allres=zeros(2,numel(a));
% for i=1:numel(a)
%     res=floor(a(:,i).Centroid);
%     if numel(res)==2
%         allres(:,i)=res;
%     end
% end
%     figure(1)
%     imshow(rawimg);
%     hold
% plot(allres(1,:),allres(2,:),'+r')
Images1 = (regionprops(img,'Image'));
Indices=(regionprops(img,'PixelList'));
subimage=regionprops(img,'SubarrayIdx');


imshow(Images1(i).Image)
% Images2 = (regionprops(img2,'Image'));
% Images = [Images1;Images2];
for i = pos:size(Inputs1,1)
    if Inputs1(i,2) * Inputs1(i,3) > 1000
    i/size(Inputs1,1)
    %imshow(Images1(i).Image)
    twolists=subimage(i).SubarrayIdx;
    my_xs=[min(subimage(i).SubarrayIdx{1}):max(subimage(i).SubarrayIdx{1})];
    my_ys=[min(subimage(i).SubarrayIdx{2}):max(subimage(i).SubarrayIdx{2})];
    imshow(rawimg(my_xs,my_ys));
    Targets(i) = input(prompt);
    else
    Targets(i) = 0;
    end
    pos = pos + 1;
% save('TargetsAndyAxon1.mat','Targets')
% save('Inputs1AndyAxon1.mat','Inputs1')
Targets2=Targets';
%had set 13
%now set cvn
save(['Targets' num2str(mag) 'magAxoplasmSetCVN.mat'],'Targets2')
save(['Inputs' num2str(mag) 'magAxoplasmSetCVN.mat'],'Inputs1')
end
end
end