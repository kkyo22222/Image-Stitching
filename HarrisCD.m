function HarrisCD( image ,window_size,k )
%HARRISCD Summary of this function goes here
%   Detailed explanation goes here

height =size(image,1)
width =size(image,2)

biimage = rgb2gray(image);
%imshow(biimage);

[biimage_GX biimage_GY] = gradient(double(biimage));

%imshow(biimage_GX);
%pause();
%imshow(biimage_GY);
G= fspecial('gaussian',[5 5],2);

SX2=imfilter(biimage_GX.^2,G);

SY2=imfilter(biimage_GY.^2,G);

SXY=imfilter(biimage_GX.*biimage_GY,G);

resultimg=zeros(height,width);
resultimg_thre=zeros(height,width);

%image(°ª,¼e)
for h=1:height
    for w=1:width
        nowM=[SX2(h,w) SXY(h,w);SXY(h,w) SY2(h,w)];
        nowR=det(nowM)-k*trace(nowM);
         
         if nowR>80
             resultimg(h,w) =nowR;
         else
             resultimg(h,w) =0;
         end
        
        
    end
end

%find local maximum
interval=10;

%imshow(resultimg);
%pause();

for h=1:height
    
    for w=1:width
        flag=0;
        nowValue = resultimg(h,w);
        if nowValue ==0
            flag=1;
        end
        for i = -interval:interval
            if flag==1
                    break;
            end
            for j = -interval:interval
                if flag==1
                    break;
                end
                
                if resultimg( max(min(h+i,height),1) ,max(min(w+j,width),1) ) > nowValue
                    flag=1;
                end
                
            end
        end
        
        if flag==0
            resultimg_thre(h,w)=1;
        else
            resultimg_thre(h,w)=0;
        end
        
    end
end

imshow(image);
hold on
for h=1:height
    for w=1:width
        if resultimg_thre(h,w)==1
            plot(w,h,'r.','Markersize',7);
        end
    end
end


end

