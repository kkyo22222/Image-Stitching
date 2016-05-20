function  showkey( append_direction, image,figuren )
%SHOWKEY Summary of this function goes here
%   Detailed explanation goes here


    if ~exist('figuren')
       figuren = 1;
    end

    figure(figuren);
    
imshow(image)

    num_bins = 36;
    hist_step = 2*pi/num_bins;

    hold on    
    for i=1:size(append_direction,1)
        %line([ฐ_??งค?กA???งค?],[ฐ_??งค?กA???งค?])กA
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
   

end

