function mygratios=gratio_noshow(mybwimage, rawimg_small, mycentroids, mydiameters, myorientation)

mygratios=zeros(numel(mydiameters),1);
%calculate and return g ratio for all axons returned from NN 
%sample call
%mygratios=gratio(mybwimage, rawimg_small,mycentroids, mydiameters, myorientation);
simage=size(mybwimage);
sx=simage(1);
sy=simage(2);

parpool(4);
parfor i=1:numel(mydiameters)
   
    mycenter=mycentroids(i,:);
    mydiameter=mydiameters(i);
    myx=floor([max(1,mycenter(1)-1.2*mydiameter/2):min(sx,mycenter(1)+1.2*mydiameter/2)]);
    myy=floor([max(1,mycenter(2)-1.2*mydiameter/2):min(sy,mycenter(2)+1.2*mydiameter/2)]);
    
    myx_ys=repmat(mycenter(2), numel(myx),1);
    myy_xs=repmat(mycenter(1), numel(myy),1);
        
        %add more lines
        %theta=45;135 at least
        %addx=floor(max(1,[mycenter(1)-0.6*mydiameter*cosd(theta)]):min(sy,[mycenter(1)+0.6*mydiameter*cosd(theta)]));
        %addy=floor(max(1,[mycenter(2)-0.6*mydiameter*sind(theta)]):min(sy,[mycenter(2)+0.6*mydiameter*sind(theta)]));
        
        %find intersection with distance map
        [D,IDX] = bwdist(mybwimage,'euclidean');%was ~bw2
   
    %start temp
%     figure(2)
%       imagesc(rawimg_small)
%       colormap gray
%         hold on
%         plot(mycenter(1,1), mycenter(1,2), 'g*')
%         
%        plot(myx,myx_ys, 'c');
%        plot(myy_xs,myy, 'm');
        
        
        
        
%         hold off
%         title(['my ' num2str(i) 'centroid and ring after NN']);
        %end temp
        %find intersection with distance map
        
        
   
%    figure(2)
%    mymap=jet(24);
%    mymap(1,:)=[0 0 0];
%    colormap(mymap)
%    imagesc(D)
%    hold
%    plot(myx,myx_ys, 'c');
%    plot(myy_xs,myy, 'm');
%    title('intersect distance map (Euclidian)')
%    hold off
   
%    figure(3)
   myx=uint16(myx);
   myx_ys=uint16(myx_ys);
   myy=uint16(myy);
   myy_xs=uint16(myy_xs);
   
   ImXvals=zeros(numel(myx,1));ImYvals=zeros(numel(myy,1));
   xcoords=[myx',myx_ys];
  for j=1:numel(myx)
      ImXvals(j)=D(xcoords(j,2),xcoords(j,1));
  end
  
   ycoords=[myy_xs,myy'];
   
   for j=1:numel(myy)
       ImYvals(j)=D(ycoords(j,2),ycoords(j,1));
   end
   
%    figure(4)
   sImXvals = smooth(myx,ImXvals,0.05,'rloess');
   sImYvals=smooth(myy,ImYvals,0.05,'rloess');
%    subplot(1,2,1)
%    plot(myx,ImXvals, 'black');
%    hold
%    plot(myx,sImXvals, 'c');
%    hold
%    title(['x line profiles for ' num2str(i) 'centroid (Euclidian)']);
%    legend('original', 'X profile');
%    
%    subplot(1,2,2)
%    plot(myy,ImYvals, 'black');
%    hold
%    plot(myy,sImYvals, 'm');
%    title(['y line profiles for ' num2str(i) 'centroid (Euclidian)']);
%    hold
%    legend('original', 'Y profile');
   
   

   myfx=diff(smooth(ImXvals,3));
   myfxleft=circshift(myfx,-1,1);
   myfxright=circshift(myfx,1,1);
   
   myfy=diff(smooth(ImYvals,3));
   myfyleft=circshift(myfy,-1,1);
   myfyright=circshift(myfy,1,1);
   
   
   [xmyelin_ind,xmyelin_vals]=find(myfxleft.*myfxright<0);%finds max in distance map
   [ymyelin_ind,ymyelin_vals]=find(myfyleft.*myfyright<0);%finds max in distance map
  
   xmthick=ImXvals(xmyelin_ind);
   ymthick=ImYvals(ymyelin_ind);
   
   myelinthick=2*median([xmthick,ymthick]);
   gration1=mydiameter./(mydiameter+2*myelinthick);
   mygratios(i)=gration1;
   %[ymyelin_ind,ymyelin_vals]=find(diff(smooth(ImXvals,3)));
   
% plot(diff((diff(smooth(ImXvals,3)))))
% hold
%
% plot(diff(smooth(ImXvals,3)))
% legend('diff2', 'diff1')
% 
%    maxy=max(yy2);
%    miny=min(yy2);
        
        
       
    
    
end



delete(gcp('nocreate'));

