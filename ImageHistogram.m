%% Author: Colbert Sesanker, 2/8/2011
%% This Matlab funtion returns a Histogram plot of 
%% a grey scale image with 256 grey levels

function [x y] = simple_hist(image)

 im = imread(image);
 L = 256;
 y = zeros(1,L);             %% Stores # of pixels with the kth grey level
 x = 0:255;                  %% Grey Scale Levels
 [width, height] = size(im);
 for i=1: width
     for j = 1: height
            for k = 1:L
         if im(i,j) == x(k)
             y(k) = y(k)+1;
         end
                    
             end    
     end
 end
 
bar(x,y);

