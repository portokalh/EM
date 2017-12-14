clear all
train=1;
pos = 1;
prompt = 'Axons:';
mag=2650;
filein='/Users/alex/Desktop/andy/train_140228_1_NOS_2650.TIF'; %'/Users/alex/Documents/MATLAB/Axon1.TIF';

mag=7100;
filein='/Users/alex/Desktop/andy/train_140228_1_NOS_7100.TIF'; %
rawimg = imread(filein); %'TestImg_CHT_a2.bmp');
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
[img, init_no] = mommy(filein,mag);
 
close all
%img = imcomplement(img);
props = regionprops(img,'Area','Eccentricity','EquivDiameter','EulerNumber',...
    'Extent','FilledArea','MajorAxisLength','MinorAxisLength',...
    'Orientation','Perimeter','Solidity')

for ii=1:(size(props,1))
    FilledPercent(ii)= props(ii).Area./props(ii).FilledArea;
end
for ii=1:size(props,1)
    HoleArea(ii)= props(ii).Area-props(ii).FilledArea*-1;
end
for ii=1:size(props,1)
    Geomean(ii)= geomean([props(ii).MajorAxisLength,props(ii).MinorAxisLength]);
end
for ii=1:size(props,1)
    PerArea(ii)= (props(ii).Perimeter./props(ii).FilledArea);
end

Inputs1 = [props.Area; props.MajorAxisLength; props.MinorAxisLength; props.Eccentricity;...
    props.Orientation; props.FilledArea;props.EulerNumber; props.EquivDiameter;...
    props.Solidity; props.Extent; props.Perimeter;FilledPercent;...
    HoleArea;Geomean;PerArea].';
Axons =  (round(abs(Axoner(Inputs1')),1));
Asum = sum(Axons)

j=1;
g=0;

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
    %figure(10);
    %imshow(rawimg);
    %hold
%plot(allres(1,:),allres(2,:),'+r');
%hold

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
save(['Targets' num2str(mag) 'magClean.mat'],'Targets')
save(['Inputs' num2str(mag) 'magClean.mat'],'Inputs1')
end
end
end