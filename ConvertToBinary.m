%%ConvertToBinary.m
% Converts a full matrix of elements (which values are limited from 0 to 1)
% in a binary matrix, comparing every cell value with a user given
% threshold, and print a grayscale colormap of the binary matrix.
% Threshold value must be contained between 0 and 1.
% USED IN topmain.m
function [v]=ConvertToBinary(x,th)
global nelx nely

v = zeros(nely,nelx);
v(x>=th) = 1;

% Outliers filter
for ii = 2:nely-2
    for jj = 2:nelx-2
        if v(ii,jj)==0 && v(ii-1,jj)==1 && v(ii,jj+1)==1 && v(ii+1,jj)==1 && v(ii,jj-1)==1
            v(ii,jj)=1;
        end
    end
end 

figure
colormap(gray); imagesc(-v); axis equal; axis tight; axis off;pause(1e-6);
