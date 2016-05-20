function [xbar ybar] = cwarp( x,y,centralx,centraly,f,s,direct)
%CWARP Summary of this function goes here
%   Detailed explanation goes here

    if ~exist('direct')
       direct = 1;
    end
    
    if direct ==1
        xbar1= f*tan( (x-centralx) /s)    ;
        ybar1= sqrt((x-centralx)^2+f^2)/s*(y-centraly);
        xbar = round(xbar1)+centralx;
        ybar = round(ybar1)+centraly;
    elseif direct ==0;
        xbar1 = s*atan2((x-centralx),f);
        ybar1 = s*(y-centraly)/sqrt((x-centralx)^2+f^2);
        xbar = round(xbar1)+centralx;
        ybar = round(ybar1)+centraly;
    end

end

