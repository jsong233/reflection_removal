function [y] = isparsePrior(x)
[h,w,ch] = size(x); y = zeros(h,w);
f1 = [-1,1]'; f2 = [-1;1]'; f3 = [0 1 0; 1 -4 1; 0 1 0]';
y = y + imfilter(x(:,:,1),f1);
y = y + imfilter(x(:,:,2),f2);
y = y + imfilter(x(:,:,3),f3);
end

