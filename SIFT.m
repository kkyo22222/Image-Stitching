function [ points_with_desc append_direction ] = SIFT( image,octaves,intervals )
%SIFT Summary of this function goes here
%   Detailed explanation goes here

debug =0;

%intervals = s 代表每個octaves中的間隔

    if ~exist('octaves')
       octaves = 4;
    end

    if ~exist('intervals')
       intervals = 3;
    end
    
    contrast_thre=0.02;
    edge_thre =10.0;
    edge_thre_ratio =(edge_thre+1)^2 /edge_thre;
    %在做偏微分時所用的kernel
    xx = [ 1 -2  1 ];
    yy = xx';
    xy = [ 1 0 -1; 0 0 0; -1 0 1 ]/4;
    
    %先判斷讀入的image的類型
    
    if isgray(image) ==0
        image = rgb2gray(image);
    end
    if isfloat(image) ==0
        image=double(image)./255;
    end

    %原圖先用0.5sigma去blur放大兩倍在做presmoothing sigma=1.6 can be robost
    height = size(image,1);
    width  = size(image,2);
    
    hsize = GKernelSize(0.5);
    gfilter =fspecial('gaussian',hsize,0.5);
    image =imfilter(image,gfilter);     

    
    imgtemp = imresize(image,[height*2 width*2],'nearest');    
    hsize = GKernelSize(1);
    gfilter =fspecial('gaussian',hsize,1);
    
    imgtemp =imfilter(imgtemp,gfilter); 
    
    ori_sigma =1;
%在同一個octave中套上不同的sigma~ 要有interval+3個圖
    gau_imgset=cell(octaves,intervals+3);
    dog_imgset=cell(octaves,intervals+2);
    truly_sigma=zeros(octaves,intervals+3);
    
    gau_imgset{1,1}=imgtemp;
    
    
    for i=1:octaves        
        %fprintf('octaves %d\n',i);
        nowsigma = ori_sigma;       
        nowsigma = nowsigma* 2^(-1/intervals);
        abs_sigma(i,1)=nowsigma*2^(i-1);
        %fprintf('\t\tinterval:%d \tsigma %f\n',1,nowsigma*2^(i-1));
        for j=2:intervals+3
            %gau_imgset{i,j} =zeros(height*2^(2-i),width*2^(2-i));
            nowsigma = nowsigma* 2^(1/intervals);
            truly_sigma(i,j)=nowsigma;
            abs_sigma(i,j)=nowsigma*2^(i-1);
            %fprintf('\t\tinterval:%d \tsigma %f\n',j,nowsigma*2^(i-1));           
            hsize = GKernelSize(nowsigma* 2^(1/intervals));
            gfilter =fspecial('gaussian',hsize,nowsigma* 2^(1/intervals));
            %其實每一張都可以從前一張套sigma得出
            gau_imgset{i,j}=imfilter(gau_imgset{i,j-1},gfilter);   
        end
        if(i<4)
            gau_imgset{i+1,1} = imresize(gau_imgset{i,intervals+1},[height*2^(1-i) width*2^(1-i)],'nearest');    
        end
       
    end
    
%abs_sigma    
    
%      for i=1:octaves
%          for j=1:intervals+3
%              subplot(octaves,6,(i-1)*6+j);
%              imshow(gau_imgset{i,j});
%          end
%      end
  %製造一張大圖片來show gau的結果
  if debug==1
    y=1;
    for i=1:octaves
        hh = height*2^(2-i);
         y = y+hh; 
    end
    
    showimg=zeros(y,(intervals+3)*width*2);
    y=1;
    for i=1:octaves        
         for j=1:intervals+3
             hh = height*2^(2-i);
             ww = width*2^(2-i);
             x=1+ (j-1)*ww;             
             showimg(y:y+hh-1,x:x+ww-1)=gau_imgset{i,j};
         end
        y=y+hh; 
    end   
    imshow(showimg);
    pause();
  end

    for i=1:octaves
        for j=1:intervals+2
            dog_imgset{i,j} = gau_imgset{i,j+1}-gau_imgset{i,j};
        end
    end
    
%     for i=1:octaves
%          for j=1:intervals+2
%              subplot(octaves,5,(i-1)*5+j);
%              imshow(dog_imgset{i,j});
%          end
%     end
   %製造一張大圖片來show dog的結果
   if debug==1
    y=1;
    for i=1:octaves
        hh = height*2^(2-i);
         y = y+hh; 
    end
    
    showimg=zeros(y,(intervals+2)*width*2);
    y=1;
    for i=1:octaves        
         for j=1:intervals+2
             hh = height*2^(2-i);
             ww = width*2^(2-i);
             x=1+ (j-1)*ww;             
             showimg(y:y+hh-1,x:x+ww-1)=dog_imgset{i,j};
         end
        y=y+hh; 
    end   
    imshow(showimg);  
    pause();
   end
   candidatelist =cell(octaves,1);
   after_contrast{i} =cell(octaves,1);
   after_edge=cell(octaves,intervals+1);
   for i=1:octaves 
    candidatelist{i}=zeros(10000,2);
    after_contrast{i} =zeros(10000,2);
    for j=2:intervals+1      
        after_edge{i,j} =zeros(10000,2);
    end
   end

   
 disp('find points'); 

   pointlist=zeros(3,9);
   for i=1:octaves 
       countc=1;
       countac=1;

       hh = height*2^(2-i);
       ww = width*2^(2-i);
       for j=2:intervals+1          
           countae=1;
           redge= ceil(GKernelSize(truly_sigma(i,j))/2)  ;                  
           for x=2+redge:ww-1-redge               
               for y=2+redge:hh-1-redge
                   %找26個點來比
                   %中心是dog_imgset{i,j}(y,x) pointlist(2,5)
                   pointlist(1:3,1:3)= dog_imgset{i,j-1}(y-1:y+1,x-1:x+1);
                   pointlist(1:3,4:6)= dog_imgset{i,j}(y-1:y+1,x-1:x+1);
                   pointlist(1:3,7:9)= dog_imgset{i,j+1}(y-1:y+1,x-1:x+1);                 
                   
                   if  (pointlist(2,5)==min(pointlist(:)) ) || (pointlist(2,5)==max(pointlist(:)))
                       %加到候選人中
                       candiexist=0;
%                        for Ncandi=2:size(candidatelist{i},1) 
%                             if candidatelist{i}(Ncandi,1)==x && candidatelist{i}(Ncandi,2) ==y   
%                                 candiexist=1;
%                             end
%                        end
                       if candiexist ==0
                         if debug==1
                             candidatelist{i}(countc,1) = x;
                             candidatelist{i}(countc,2) = y;
                             countc =countc+1;
                         end
                         if abs(dog_imgset{i,j}(y,x)) >= contrast_thre
                             if debug==1
                                 after_contrast{i}(countac,1) = x ;
                                 after_contrast{i}(countac,2) = y ;
                                 countac =countac+1;
                             end
                             DXX = sum(dog_imgset{i,j}(y,x-1:x+1).*xx );
                             DYY = sum(dog_imgset{i,j}(y-1:y+1,x).*yy );
                             DXY = sum(sum(dog_imgset{i,j}(y-1:y+1,x-1:x+1).*xy ));
                             
                             TrH =  DXX+DYY;
                             DetH = DXX*DYY-DXY^2;
                             if ((TrH^2/DetH) < edge_thre_ratio) 
                                 after_edge{i,j}(countae,1) = x;  
                                 after_edge{i,j}(countae,2) = y;
                                 countae =countae+1;
                             end
                         end
                       end
                   end                   
               end
           end
           after_edge{i,j} = after_edge{i,j}(1:countae-1,:);
       end
   end
 %{  
   candidatelist
   size(image)
   imshow(image)
   hold on
   for i=1:octaves 
       for Ncandi=2:size(candidatelist{i},1)     
           plotx = (candidatelist{i}(Ncandi,1)-1)*2^(i-2)+1;
           ploty = (candidatelist{i}(Ncandi,2)-1)*2^(i-2)+1;          
           plot( plotx,ploty,'y+','Markersize',5)
       end
   end
   hold off
   pause();

   after_contrast
   size(image)
   imshow(image)
   hold on
   for i=1:octaves 
       for Ncandi=2:size(after_contrast{i},1)     
           plotx = (after_contrast{i}(Ncandi,1)-1)*2^(i-2)+1;
           ploty = (after_contrast{i}(Ncandi,2)-1)*2^(i-2)+1;          
           plot( plotx,ploty,'y+','Markersize',5)
       end
   end
   hold off
   pause();

   
   after_edge
   size(image)
   imshow(image)
   hold on
   for i=1:octaves 
       for Ncandi=2:size(after_edge{i},1)     
           plotx = (after_edge{i}(Ncandi,1)-1)*2^(i-2)+1;
           ploty = (after_edge{i}(Ncandi,2)-1)*2^(i-2)+1;          
           plot( plotx,ploty,'y+','Markersize',5)
       end
   end
   hold off
   pause();
 %}  
  
   
   disp('calculate direction')
   %after_edge
   %對之前gau的圖片找其gradient與方向 這是預處理 節省時間
    for i=1:octaves
        for j=2:intervals+1
            [mag_x{i,j} mag_y{i,j}] =gradient (gau_imgset{i,j}) ; 
            mag_sum{i,j} = sqrt(mag_x{i,j}.^2+mag_y{i,j}.^2);
            dire{i,j} = atan2(mag_y{i,j},mag_x{i,j});
            dire{i,j}(find(dire{i,j} == pi)) = -pi;            
        end
    end
    
    %對每一個scale的keypoint去用一window(1.5該scale的sigma)計算裡面的方向累積圖 ,計算最高的方向
    num_bins = 36;
    hist_step = 2*pi/num_bins;
    hist_orient = [-pi:hist_step:(pi-hist_step)];
    
    append_direction=zeros(1,4);
    % x y scale 方向

    for i=1:octaves
        for j=2:intervals+1
            nowsigma =truly_sigma(i,j)*2;
            hsize = GKernelSize(1.5*nowsigma);
            gfilter =fspecial('gaussian',hsize,1.5*nowsigma);
            hrange = (hsize-1)/2;
            hh = height*2^(2-i);
            ww = width*2^(2-i);            
            for Ncandi=2:size(after_edge{i,j},1)
                %相對於原始圖的座標
                %nowx= after_edge{i,j}(Ncandi,1)*2^(i-2)+1;
                %nowy= after_edge{i,j}(Ncandi,2)*2^(i-2)+1;
                nowx= after_edge{i,j}(Ncandi,1);
                nowy= after_edge{i,j}(Ncandi,2);
                truex= (nowx-1)*2^(i-2)+1;
                truey= (nowy-1)*2^(i-2)+1;
              
                %對keypoint附近的做fliter
                now_magregion =zeros(hsize,hsize);
                now_dire = zeros(hsize,hsize);
                for m=-hrange:hrange
                    for n= -hrange:hrange
                         if nowy+m>=1 && nowy+m <= hh && nowx+n>=1 && nowx+n <= ww 
                            now_magregion(hrange+m+1,hrange+n+1) = mag_sum{i,j} (nowy+m,nowx+n);
                            now_dire(hrange+m+1,hrange+n+1) = dire{i,j}(nowy+m,nowx+n);
                         end
                    end
                end               
               
                each_point_vote_weight =gfilter.*now_magregion;               
                
                now_hist =zeros(num_bins,1);
                for m=1:hsize
                    for n= 1:hsize
                        %bin是-pi,-pi+step...,為了對應回去 所以+pi+1 以投到正確位置
                        tmpvalue = (now_dire(m,n)+pi)/hist_step;
                        tmpvalue =mod(tmpvalue,36)+1;
                        upb = mod(ceil(tmpvalue)-1,36)+1; 
                        lowb = floor(tmpvalue);
                        ratio = tmpvalue -lowb;
                        now_hist(upb) = now_hist(upb) + each_point_vote_weight(m,n)*(ratio);
                        now_hist(lowb) = now_hist(lowb) + each_point_vote_weight(m,n)*(1-ratio);  
                    end
                end
                now_hist;
                [maxvalue maxpos]=max(now_hist);
                nowpeak =maxvalue;
                while nowpeak >0.8*maxvalue
                    %把這個方向加到keypoint擁有的方向
                    
                    A = []; %a是現在跟旁邊peak的位置
                    b = []; %b對應到其值
                    for jj = -1:1
                       A = [A; (maxpos+jj).^2  maxpos+jj 1];
                        binv= mod(maxpos+jj-1,num_bins)+1;
                       b = [b; now_hist(binv)];
                    end
                    c = pinv(A)*b;  
                    max_orient = -c(2)/(2*c(1));
                    while( max_orient < 0 )
                        max_orient = max_orient + 36;
                    end
                    while( max_orient >= 36 )
                       max_orient = max_orient - 36;
                    end            
                               
                    
                    append_direction =[append_direction; truex truey truly_sigma(i,j)*2^(i-1) max_orient];
                    now_hist(maxpos)=0;
                    [nowpeak maxpos]=max(now_hist);
                end
                
                
                %看各個對應點的方向  來投票到對應的hist上
            end
        end
    end
   if debug ==2
    imshow(image)

    hold on
    for i=2:size(append_direction,1)
        %line([起??坐?，???坐?],[起??坐?，???坐?])，
         nscale = append_direction(i,3);
         nscale = 25/nscale;
         ndire =  append_direction(i,4);
         ndire = (ndire-1)*hist_step;
%         goal_x = append_direction(i,1)+ nscale*cos(ndire)
%         goal_y = append_direction(i,2)+ nscale*sin(ndire)        
%         line([append_direction(i,1) goal_x],[append_direction(i,2) goal_y],'Color','y');       
        aaa=quiver(append_direction(i,1),append_direction(i,2),nscale*cos(ndire),nscale*sin(ndire) );
        set(aaa,'color',[1 1 0]);
        plot(append_direction(i,1),append_direction(i,2),'r.');         
    end
    
    
    
    hold off
    pause();  
   end
     
   
   
   
    %%產生descriptor
    disp('make descriptor');
  
    
orient_bin_spacing = pi/4;
orient_angles = [-pi:orient_bin_spacing:(pi-orient_bin_spacing)];

now_coord= zeros(16,16,2);

   for ii=1:16
       for jj=1:16
           now_coord(ii,jj,1)= jj-8.5;
           now_coord(ii,jj,2)= ii-8.5;
       end
   end

gfilter =fspecial('gaussian',16,0.5*16); 

points_with_desc =[];
   
for i=2:size(append_direction,1)
   nscale = append_direction(i,3);
   ndire =  append_direction(i,4);
   ndire = (ndire-1)*hist_step;
   
   [octa inter] = find(abs_sigma==nscale);  
   if size(octa,1)>=2
       octa=octa(1);
       inter=inter(1);
   end
   
   x = append_direction(i,1)*2^(2-octa); %轉為最接近尺度的座標
   y = append_direction(i,2)*2^(2-octa);
   hh = round(height*2^(2-octa));
   ww = round(width*2^(2-octa));
   M = [cos(ndire) -sin(ndire); sin(ndire) cos(ndire) ];
   
   feat_desc = [];
   for ii=1:16
       for jj=1:16
           now_coord(ii,jj,1)= jj-8.5;
           now_coord(ii,jj,2)= ii-8.5;
       end
   end
   
   %將對於keypoint的offset做旋轉
   for ii=1:16
       for jj=1:16
           now_coord(ii,jj,:) = M*[now_coord(ii,jj,1);now_coord(ii,jj,2)]; 
       end
   end
   tmpimage =zeros(16,16);
   for ii=1:16
       for jj=1:16
           shiftx =now_coord(ii,jj,1);
           shifty =now_coord(ii,jj,2);
           tmpimage(jj,ii)=gau_imgset{octa,inter}(min(max(round(y+shifty),1),hh), min(max(round(x+shiftx),1),ww) ); 
       end
   end
   %imshow(tmpimage);
   
   [tmppx tmppy]= gradient(tmpimage);
   mag_tmp = sqrt(tmppx.^2+tmppy.^2);
   orit_tmp = atan2(tmppy,tmppx);
   
   %gfilter =fspecial('gaussian',16,0.5*16);
   true_mag = gfilter.*mag_tmp;
   
   for feature_h=1:4
       for feature_w=1:4
           now_orient=zeros(1,8);
           for hhh=1:4
               for www=1:4
                   yyy =4*(feature_h-1)+hhh;
                   xxx =4*(feature_w-1)+www;
                   tmpori = orit_tmp(yyy ,xxx);
                   if tmpori==pi
                       tmpori=-pi;
                   end
                   votebin = ((tmpori+pi)/orient_bin_spacing) ;
                   votebin = mod(votebin,8)+1;                   
                   upb = mod(ceil(votebin)-1,8)+1 ;               
                   lowb = floor(votebin);
                
                   
                   ratio = votebin -lowb;
                   now_orient(upb) = now_orient(upb) + ratio* true_mag(yyy,xxx)*100;
                   now_orient(lowb) = now_orient(lowb) + (1-ratio)*true_mag(yyy,xxx)*100; 
                   
               end
           end
           feat_desc =[feat_desc now_orient];           
       end
   end
   feat_desc= feat_desc/norm(feat_desc);
   feat_desc( find(feat_desc > 0.2) ) = 0.2;
   feat_desc = feat_desc / norm(feat_desc);
   points_with_desc=[points_with_desc; feat_desc];
   %gau_imgset{octa,inter}()
   %now_coord
   %pause
    

end

append_direction= append_direction(2:size(append_direction,1),:);
   

end

