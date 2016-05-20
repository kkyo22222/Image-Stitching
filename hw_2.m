% parameter region--------------------------------
reloadsift=1;
blendingmethod =2;
% 0 direct blending 1 poisson blending 2 linear blending
bundleadjustflag =1;
%0 不做 1做
clearflag =0;
imagespath = './images/grail/';
imagestype = 'jpg';
f=2000;
s=2000;
show_match_pair=0;
%-------------------------------------------------

if ~exist('results', 'dir')
        mkdir('results');
end


tic
%先選定要stiching的資料夾
list=dir([imagespath '*.' imagestype]);
%把每張圖片讀入並計算其descriptor
imset =cell(size(list,1),1);
if reloadsift==1
desc = cell(size(list,1),1);
pos = cell(size(list,1),1);
warppos =cell(size(list,1),1);
end

for i=1:size(list,1)
    imset{i}=imread([imagespath list(i).name]);
    fprintf('image %d calculate sift now\n',i);
    if reloadsift==1
    [desc{i} pos{i}]= SIFT(imset{i});
    end
end
%對每張圖片做warping

disp('warping image and points');
for i=1:size(list,1)
    [imset{i} cutx] =(cwarpimage(imset{i},f,s));
    imset{i} =uint8(imset{i});
%對pos也做warping
    centralx= round((1+size(imset{i},2))/2)+cutx ;
    centraly= round((1+size(imset{i},1))/2) ;
    for j=1:size(pos{i},1)
        [warppos{i}(j,1) warppos{i}(j,2)] = cwarp(pos{i}(j,1),pos{i}(j,2),centralx,centraly,f,s,0);
        warppos{i}(j,3)=pos{i}(j,3);
        warppos{i}(j,4)=pos{i}(j,4);
        warppos{i}(j,1)=warppos{i}(j,1)-cutx;       
    end

end

bigkdtreept =[];

for i=1:size(list,1)
    bigkdtreept=[bigkdtreept desc{i}'];
    eachkdtree{i} = ann(desc{i}');
end

bigkdtree =ann(bigkdtreept);
clear bigkdtreept;
bigkdindex =zeros(size(list,1),1);
bigkdindex(1) = size(desc{1},1);

for i=2:size(desc,1)
bigkdindex(i)=bigkdindex(i-1)+size(desc{i},1);
end


imgmatchscore=zeros(size(list,1),size(list,1));

%每張圖計算和其他圖feature的match程度
disp('calculate similarity');
for i=1:size(list,1)    
    for j=1:size(desc{i},1)    
        pt = desc{i}(j,:)';        
         [idx dst] =ksearch(bigkdtree,pt, 5, 1);
          %利用idx看是跟哪一張圖match
         for k=1:5
             for searchi=1:size(bigkdindex,1)
                if idx(k) < bigkdindex(searchi) 
                    break
                end
             end
         %idx找到後再做一次確認
             [subidx subdst] = ksearch(eachkdtree{searchi},pt, 2, 1);
             if subdst(1)/subdst(2) <0.7
                imgmatchscore(i,searchi) =  imgmatchscore(i,searchi)+1;
             end
         end
        
    end  
end

transmodel =cell(size(list,1),size(list,1));
disp('calculate the model');
for i=1:size(list,1)    
    disp(i);
    
    %找出除了自己以外6張相像的
    nowscore = imgmatchscore(i,:);
    [scorebar sortidx ]=sort(nowscore,'descend');
    mostlikeflag=0;
    for j=1:7
        if sortidx(j) ~=i && imgmatchscore(i,sortidx(j)) >6
         
            
            % ransac
            imgshiftx=0;
            imgshifty=0;    
            
            %現在是第i張圖 對第sortidx(j)圖
            pointlistdata =[];
            bestmodel=[];
            besterror =realmax('double');
            bestinlier =0;
            
            if show_match_pair==1 &&mostlikeflag==0
                 showshowimg =zeros(size(imset{i},1) ,size(imset{i},2)+size(imset{sortidx(j)},2),3);
                 showshowimg(1:size(imset{i},1),1:size(imset{i},2),:) =imset{i};
                 showshowimg(1:size(imset{i},1),size(imset{i},2)+1:size(imset{i},2)+size(imset{sortidx(j)},2),: ) =imset{sortidx(j)};
                 imshow(uint8(showshowimg));
                 hold on
            end
            
            for ii = 1:size(desc{i},1)                
                pt = desc{i}(ii,:)';
                [subidx subdst] = ksearch(eachkdtree{sortidx(j)},pt, 2, 1);
                if subdst(1)/subdst(2) <0.7
                    %儲存keypoint間的關係
                    ptbar = desc{sortidx(j)}(subidx(1),:)';
                    [subidxbat subdstbar] = ksearch(eachkdtree{i},ptbar, 2, 1);
                    if subdstbar(1)/subdstbar(2) <0.7
                        pointlistdata =[pointlistdata; warppos{i}(ii,1),warppos{sortidx(j)}(subidx(1),1)  ,warppos{i}(ii,2),warppos{sortidx(j)}(subidx(1),2) ];
                    end
                    if show_match_pair ==1 &&mostlikeflag==0
                    	line([warppos{i}(ii,1) warppos{sortidx(j)}(subidx(1),1)+size(imset{i},2) ],[warppos{i}(ii,2) warppos{sortidx(j)}(subidx(1),2)]);
                    end  
                end
            end
            
            if show_match_pair ==1 &&mostlikeflag==0
                mostlikeflag=1;
                outputname = ['imgpair_' num2str(i) '_' num2str(sortidx(j)) '.jpg'];
                print('-djpeg',outputname);
            end
            
            
            %pointlistdata 圖對圖對應的match的位移
            pointlistdata;
            bbb = size(pointlistdata,1);           
            if bbb==0
                continue;
            end
            for iter=1:300 
                %隨機挑6個點
                randidx = fix(random('uniform',1,bbb,6,1));                
                selectpoint = pointlistdata(randidx,:);                
                %nowmodel =mean(selectpoint);
                dx=0;
                dy=0;
                for kkp=1:6
                    dx = dx + (selectpoint(kkp,2)-selectpoint(kkp,1));
                    dy = dy + (selectpoint(kkp,4)-selectpoint(kkp,3));                   
                end
                %Ax=b  x是parameter

                nowmodel =[dx/6 dy/6];
                
                
                nowerror =0;
                countinlier=0;
                for jj=1:bbb
                    %判斷是否為outlier 不管outlier的error                    
                    %tmperror = sqrt(sum((pointlistdata(jj,:)-nowmodel).^2));
                    predictp = [pointlistdata(jj,1)+nowmodel(1) pointlistdata(jj,3)+nowmodel(2)];
                    tmperror = sqrt((predictp(1)-pointlistdata(jj,2))^2 + (predictp(2)-pointlistdata(jj,4))^2);                     
                    if tmperror >15
                        %即為outlier
                    else
                        nowerror = nowerror +tmperror;
                        countinlier = countinlier +1;
                    end
                    
                    
                end
                if countinlier >bestinlier
                    %取代成新的模型
                    bestmodel =nowmodel ;
                    besterror =nowerror ;
                    bestinlier=countinlier;                    
                end
            end
            bestinlier;
            if bestinlier >(5.9+0.22*bbb)
                transmodel{i,sortidx(j)}=bestmodel;
            end
           
        end
    end
end



disp('image stitching');
repeatflag =0;
%從圖片一開始往右邊找可以貼的
panoimage=imset{1};
nowsearchidx=1;

search_dire=1;
% 1向右 0向左
righttail=0;
lefttail=0;
right_totaldy=0;
left_totaldy=0;
right_totaldx=0;
left_totaldx=0;
while repeatflag~=1
    %找x大於0 且最近的
    nowcloset=realmax('double');
    stitichidx=0;
    if search_dire ==1
        for j=1:size(transmodel,1)
            if size(transmodel{nowsearchidx,j},1) >0
                if transmodel{nowsearchidx,j}(1) <0 && abs(transmodel{nowsearchidx,j}(1))< nowcloset
                    stitichidx =j;
                    nowcloset = abs(transmodel{nowsearchidx,j}(1)); 
                end
            end
        end
        stitichidx;
        if stitichidx ~=0
            h1=size(panoimage,1);  %原圖高
            h2=size(imset{stitichidx},1);  %要貼的圖的高
            w1=size(panoimage,2);
            w2=size(imset{stitichidx},2);
           
            parameter = round(transmodel{nowsearchidx,stitichidx}(:));            
            parameter(2) = parameter(2)+right_totaldy;
            parameter(1) = parameter(1)+right_totaldx;
            if parameter(2)<0
                right_totaldy= parameter(2);
            else
                left_totaldy =left_totaldy+parameter(2);
                right_totaldy=0;
            end            
            
            right_totaldx= parameter(1);
            tmpwidth =abs(parameter(1));
            if parameter(2)>0
                newh=h1+parameter(2);                
            else
                newh=max(h1,h2+abs(parameter(2)));
            end
            neww=max(w1,w2+tmpwidth);
         

            pasteimage =double(imset{stitichidx});
            if parameter(2)<0 %dy小於0 代表新圖在原圖下方
                tmppanoimage =zeros(newh,neww,3);
                tmppanoimage(1:h1,1:w1,:)=panoimage;
                panoimage = tmppanoimage;
            else
                tmppanoimage =zeros(newh,neww,3);                
                tmppanoimage(1+abs(parameter(2)):newh,1:w1,:) =panoimage;                
                panoimage = tmppanoimage;
            end
                        
%             figure(1)
%             imshow(uint8(panoimage));
%             figure(2)
%             imshow(uint8(pasteimage));
%             figure(3)
            
           if blendingmethod ==1 
                %重疊區 用poisson blending     
               if parameter(2)<0
                    panoimage(:,:,1) =Poisson(pasteimage(:,:,1),panoimage(:,:,1),[1 1 w1-tmpwidth h2],[tmpwidth+1 1+abs(parameter(2))]);
                    panoimage(:,:,2) =Poisson(pasteimage(:,:,2),panoimage(:,:,2),[1 1 w1-tmpwidth h2],[tmpwidth+1 1+abs(parameter(2))]);
                    panoimage(:,:,3) =Poisson(pasteimage(:,:,3),panoimage(:,:,3),[1 1 w1-tmpwidth h2],[tmpwidth+1 1+abs(parameter(2))]); 
                                    %Poisson(src, dst, srcPos, dstPos)
                    panoimage(1+abs(parameter(2)):abs(parameter(2))+h2,w1+1:neww,:) = pasteimage(:,w1-tmpwidth+1:w2,:);
                else
                    panoimage(:,:,1) =Poisson(pasteimage(:,:,1),panoimage(:,:,1),[1 1 w1-tmpwidth h2],[tmpwidth+1 1]);
                    panoimage(:,:,2) =Poisson(pasteimage(:,:,2),panoimage(:,:,2),[1 1 w1-tmpwidth h2],[tmpwidth+1 1]);
                    panoimage(:,:,3) =Poisson(pasteimage(:,:,3),panoimage(:,:,3),[1 1 w1-tmpwidth h2],[tmpwidth+1 1]);
                    panoimage(1:h2,w1+1:neww,:) = pasteimage(:,w1-tmpwidth+1:w2,:);
               end
           elseif blendingmethod==0 || blendingmethod==2
               if tmpwidth+w2 <w1
                    summother=w2;
               else
                    summother=w1-tmpwidth;
               end
               
               if parameter(2)<0
                   for ww=tmpwidth+1:neww
                       for hh=1+abs(parameter(2)):h2+abs(parameter(2))
                           if ww-tmpwidth <=w2 && ww <= w1
                               if blendingmethod ==0
                                panoimage(hh,ww,:) = 0.5*panoimage(hh,ww,:)+0.5*pasteimage(hh-abs(parameter(2)),ww-tmpwidth,:);
                               elseif blendingmethod ==2                                
                                blend_ratio= (ww-tmpwidth)/summother;
                                panoimage(hh,ww,:) = (1-blend_ratio)*panoimage(hh,ww,:)+(blend_ratio)*pasteimage(hh-abs(parameter(2)),ww-tmpwidth,:);   
                               end
                           elseif  ww-tmpwidth <=w2 && ww > w1 %重複區大
                            panoimage(hh,ww,:) = pasteimage(hh-abs(parameter(2)),ww-tmpwidth,:); 
                           elseif  ww-tmpwidth >w2 && ww < w1 %重複區小
                               % do nothing
                           end
                       end
                   end
               else
                   for ww=tmpwidth+1:neww
                       for hh=1:h2
                           if ww-tmpwidth <=w2 && ww <= w1
                              if blendingmethod ==0 
                                panoimage(hh,ww,:) = 0.5*panoimage(hh,ww,:)+0.5*pasteimage(hh,ww-tmpwidth,:);
                              elseif blendingmethod ==2
                                blend_ratio= (ww-tmpwidth)/summother;
                                panoimage(hh,ww,:) = (1-blend_ratio)*panoimage(hh,ww,:) + blend_ratio*pasteimage(hh,ww-tmpwidth,:);
                              end
                           elseif ww-tmpwidth <=w2 && ww > w1
                               panoimage(hh,ww,:) = pasteimage(hh,ww-tmpwidth,:); 
                           elseif ww-tmpwidth >w2 && ww < w1
                           end
                       end
                   end
               end               
           end
            
           %重疊邊界區套上gaussian filter 讓邊界變不明顯
%             hsize = GKernelSize(0.5);
%             gfilter =fspecial('gaussian',hsize,0.5);
%             if tmpwidth-4-hsize >1 && w1+5+hsize < neww
%                 tmpimg=imfilter(panoimage(:,tmpwidth-4-hsize:tmpwidth+5+hsize,:),gfilter);
%                 panoimage(:,tmpwidth-4:tmpwidth+5,:) =tmpimg(:,hsize:hsize+9,:);
%                 tmpimg =imfilter(panoimage(:,w1-4-hsize:w1+5+hsize,:),gfilter);  
%                 panoimage(:,w1-4:w1+5,:) = tmpimg(:,hsize:hsize+9,:);
%             end
%             imshow(uint8(panoimage));
            
            
            
            nowsearchidx =stitichidx;
%             imshow(uint8(panoimage))
%              pause
            righttail=stitichidx;
        else  %一直往右走卻遇到阻礙 改成從頭向左走            
            nowsearchidx=1;
            search_dire =0;
            left_totaldy = -left_totaldy;
            
        end
    else
        for j=1:size(transmodel,1)
            if size(transmodel{nowsearchidx,j},1) >0
                if transmodel{nowsearchidx,j}(1) >0 && abs(transmodel{nowsearchidx,j}(1))< nowcloset
                    stitichidx =j;
                    nowcloset = abs(transmodel{nowsearchidx,j}(1)); 
                end
            end
        end
        
        
        stitichidx;
        if stitichidx ~=0
            h1=size(panoimage,1);  %原圖高
            h2=size(imset{stitichidx},1);  %要貼的圖的高
            w1=size(panoimage,2);
            w2=size(imset{stitichidx},2);
           
            parameter = round(transmodel{nowsearchidx,stitichidx}(:));
            
            parameter(2) = parameter(2)+left_totaldy;
            left_totaldy = parameter(2);
            if parameter(2)>0
                left_totaldy=0;
                newh=h1+parameter(2);                
            else
                newh=max(h1,h2+abs(parameter(2)));
            end
            
            neww=w1+abs(parameter(1));
            tmpwidth =abs(parameter(1));

            pasteimage =double(imset{stitichidx});
            if parameter(2)<0 %dy小於0 代表新圖在原圖下方
                tmppanoimage =zeros(newh,neww,3);                
                tmppanoimage(1:h1,1+tmpwidth:neww,:)=panoimage;
                panoimage = tmppanoimage;
            else
                tmppanoimage =zeros(newh,neww,3);                
                tmppanoimage(1+abs(parameter(2)):newh,1+tmpwidth:neww,:) =panoimage;                
                panoimage = tmppanoimage;
            end
                        
%             figure(1)
%             imshow(uint8(panoimage));
%             figure(2)
%             imshow(uint8(pasteimage));
%             figure(3)
           if blendingmethod ==1 
               if parameter(2)<0
                    panoimage(:,:,1) =Poisson(pasteimage(:,:,1),panoimage(:,:,1),[tmpwidth+1 1 w2  h2],[tmpwidth+1 abs(parameter(2))]);
                    panoimage(:,:,2) =Poisson(pasteimage(:,:,2),panoimage(:,:,2),[tmpwidth+1 1 w2  h2],[tmpwidth+1 abs(parameter(2))]);
                    panoimage(:,:,3) =Poisson(pasteimage(:,:,3),panoimage(:,:,3),[tmpwidth+1 1 w2  h2],[tmpwidth+1 abs(parameter(2))]);
                                    %Poisson(src, dst, srcPos, dstPos)
                    panoimage(1+abs(parameter(2)):abs(parameter(2))+h2,1:tmpwidth,:) = pasteimage(:,1:tmpwidth,:);
                else
                    panoimage(:,:,1) =Poisson(pasteimage(:,:,1),panoimage(:,:,1),[tmpwidth+1 1 w2  h2],[tmpwidth+1 1]);
                    panoimage(:,:,2) =Poisson(pasteimage(:,:,2),panoimage(:,:,2),[tmpwidth+1 1 w2  h2],[tmpwidth+1 1]);
                    panoimage(:,:,3) =Poisson(pasteimage(:,:,3),panoimage(:,:,3),[tmpwidth+1 1 w2  h2],[tmpwidth+1 1]);

                    panoimage(1:h2,1:tmpwidth,:) = pasteimage(:,1:tmpwidth,:);
               end
           elseif blendingmethod==0 || blendingmethod==2
             
               summother=w2-tmpwidth;
              
               if parameter(2)<0
                   for ww=1:w2
                       for hh=1+abs(parameter(2)):h2+abs(parameter(2))
                           if ww>=tmpwidth+1
                               if blendingmethod ==0
                                panoimage(hh,ww,:) = 0.5*panoimage(hh,ww,:)+0.5*pasteimage(hh-abs(parameter(2)),ww,:);
                               elseif blendingmethod ==2
                                blending_ratio= (ww-tmpwidth)/summother;
                                panoimage(hh,ww,:) = blending_ratio*panoimage(hh,ww,:)+(1-blending_ratio)*pasteimage(hh-abs(parameter(2)),ww,:);   
                               end
                           else
                            panoimage(hh,ww,:) = pasteimage(hh-abs(parameter(2)),ww,:);
                           end
                       end
                   end
               else
                   for ww=1:w2
                       for hh=1:h2
                           if ww>=tmpwidth+1
                            if blendingmethod ==0   
                                panoimage(hh,ww,:) = 0.5*panoimage(hh,ww,:)+0.5*pasteimage(hh,ww,:);
                            elseif blendingmethod ==2
                                blending_ratio= (ww-tmpwidth)/summother;
                                panoimage(hh,ww,:) =  blending_ratio*panoimage(hh,ww,:)+ (1-blending_ratio)*pasteimage(hh,ww,:);
                            end
                           else
                            panoimage(hh,ww,:) = pasteimage(hh,ww,:);   
                           end
                       end
                   end
               end
           end
           
            %重疊邊界區套上gaussian filter 讓邊界變不明顯
%             hsize = GKernelSize(0.5);
%             gfilter =fspecial('gaussian',hsize,0.5);
%             if tmpwidth-4-hsize >1 && w2+5+hsize < neww
%                 tmpimg=imfilter(panoimage(:,tmpwidth-4-hsize:tmpwidth+5+hsize,:),gfilter);
%                 panoimage(:,tmpwidth-4:tmpwidth+5,:) =tmpimg(:,hsize:hsize+9,:);
%                 tmpimg =imfilter(panoimage(:,w2-4-hsize:w2+5+hsize,:),gfilter);  
%                 panoimage(:,w2-4:w2+5,:) = tmpimg(:,hsize:hsize+9,:);
%             end
            
            
%             imshow(uint8(panoimage));
           
            nowsearchidx =stitichidx;
%              imshow(uint8(panoimage))
%             pause
            lefttail=stitichidx;
        else  %一直往左走卻遇到阻礙 結束!            
            disp('end!!');
            repeatflag=1;
        end
        
        
    end    
    
    
    if righttail==lefttail || righttail == 1 || lefttail==1;
        repeatflag=1;
    end    
end
toc


hsize = GKernelSize(0.5);
gfilter =fspecial('gaussian',hsize,0.5);
panoimage = imfilter(panoimage,gfilter);
%imshow(uint8(panoimage))
%pause
imwrite(uint8(panoimage),'./results/panoimage.bmp');


%bundle adjustment

if bundleadjustflag==1
    height = size(panoimage,1);
    bestflag =0;
    best_black_region=0;

    for h=1:100
        up_cut_blackregion=0;
        down_cut_blackregion=0;
        up_cut_region=0;
        down_cut_region=0;
        for w=1:size(panoimage,2)
            hratio = (w/size(panoimage,2))*h;
            for hh=1:hratio
                up_cut_region = up_cut_region+1;
                if sum (panoimage(hh,w,:)) <50
                    up_cut_blackregion =  up_cut_blackregion +1;
                end
            end
            for hhbar=1:(h-hratio)
                down_cut_region = down_cut_region+1;
                if sum (panoimage(height-hh,w,:)) <50
                    down_cut_blackregion =  down_cut_blackregion +1;
                end
            end
        end
        if (up_cut_blackregion + down_cut_blackregion)/ (up_cut_region + down_cut_region) >0.5
            if best_black_region < (up_cut_blackregion + down_cut_blackregion);
                best_region = up_cut_region + down_cut_region;
                bset_black_region = up_cut_blackregion + down_cut_blackregion;
                besth=h;
                bestflag=1;
            end
        end

        up_cut_blackregion=0;
        down_cut_blackregion=0;
        up_cut_region=0;
        down_cut_region=0;
        for w=1:size(panoimage,2)
            hratio = (w/size(panoimage,2))*h;
            for hh=1:hratio
                up_cut_region = up_cut_region+1;
                if sum (panoimage(height-hh,w,:)) <50
                    up_cut_blackregion =  up_cut_blackregion +1;
                end
            end
            for hhbar=1:(h-hratio)
                down_cut_region = down_cut_region+1;
                if sum (panoimage(hh,w,:)) <50
                    down_cut_blackregion =  down_cut_blackregion +1;
                end
            end
        end

        if (up_cut_blackregion + down_cut_blackregion)/ (up_cut_region + down_cut_region) >0.5
            if best_black_region > (up_cut_blackregion + down_cut_blackregion);
                best_region = up_cut_region + down_cut_region;
                bset_black_region = up_cut_blackregion + down_cut_blackregion;
                besth=h;
                bestflag=2;
            end
        end

    end

    besth;
    bestflag;

    bapanoimage =zeros(size(panoimage,1)-besth,size(panoimage,2),3);
    if bestflag ==1
        for h=1:size(panoimage,1)-besth
            for w=1:size(panoimage,2)
                hratio = (w/size(panoimage,2))*besth;
                bapanoimage(h,w,:)=panoimage(round(h+hratio),w,:);
            end
        end

    elseif bestflag ==2
        for h=1:size(panoimage,1)-besth
            for w=1:size(panoimage,2)
                hratio = (1- (w/size(panoimage,2)))*besth;
                bapanoimage(h,w,:)=panoimage(round(h+hratio),w,:);
            end
        end
    end
    imshow(uint8(bapanoimage));
end
imwrite(uint8(bapanoimage),'./results/bapanoimage.bmp');

if clearflag==1
    clear
end