function [eI] = Sobel(I)
hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
eI = sqrt(Ix.^2+Iy.^2);
end

