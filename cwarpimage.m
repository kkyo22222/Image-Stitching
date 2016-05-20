function [imgresult cutx] = cwarpimage( testimg,f,s )
%CWARPIMAGE Summary of this function goes here
%   Detailed explanation goes here
imgresult =zeros(size(testimg,1),size(testimg,2),3);

centraly = round((1+size(testimg,1))/2);
centralx = round((1+size(testimg,2))/2);

changeornot =zeros(size(testimg,1),size(testimg,2));
for i=1:size(testimg,1) %h    
    for j=1:size(testimg,2) %w
        
       [xx yy]= cwarp(j,i,centralx,centraly,f,s,1);
       
       if xx<1 || xx > size(testimg,2) || yy <1 || yy>size(testimg,1)
       else
        imgresult(i,j,:) = testimg(yy,xx,:);
        changeornot(i,j) =1;
       end
    end
end


flag =0;
for j=1:size(testimg,2)  %w
    if flag ==1
        break
    end
    for i=1:size(testimg,1) %h
        %照到第一排有非0的
        if changeornot(i,j)==1
            cutx=j;
            flag=1;
            break
        end


    end
end
oriw = size(imgresult,2);
orih = size(imgresult,1);
imgresult = imgresult(1:orih,cutx:oriw-cutx,:);

end

