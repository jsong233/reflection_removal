function [y] = sparsePrior(x)
[h,w] = size(x); y = zeros(h,w,3);
f1 = [-1,1]; f2 = [-1;1]; f3 = [0 1 0; 1 -4 1; 0 1 0];
y(:,:,1) = imfilter(x,f1);
y(:,:,2) = imfilter(x,f2);
y(:,:,3) = imfilter(x,f3);
end

