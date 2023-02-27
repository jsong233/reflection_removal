function [dx, dy, c] =  estKernel(I) 
%% estimate the shift [dx dy] and the attenuation c
dim = ndims(I);
if dim == 3
    I = rgb2gray(I);
end
% Laplacian = [0 -1 0; -1 4 -1; 0 -1 0];
% resp = imfilter(I, Laplacian); % extract edge

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gI = sqrt(Ix.^2+Iy.^2);
resp = gI;
resp(abs(resp)>0.2) = 0;      % eliminate the strong edge that are likely
                           % belong to transmission layer)

auto_corr = xcorr2(abs(resp), abs(resp));
bdry = min(size(I))/2;         % eliminate the autocorrelation in boundary
                               % which is not accuracy
auto_corr = auto_corr(1+bdry:end-bdry, 1+bdry:end-bdry);

max_1 = ordfilt2(auto_corr, 25, true(5));

center = ceil(size(auto_corr)/2);
auto_corr(center(1) - 6 : center(1) + 6, center(2) - 6 : center(2) + 6) = 0; % remove the autocorrelation in the neighbour of center for the edge have width 

a = (auto_corr == max_1);
r = 5;
mask = zeros(2*r+1);
for i = 1:2*r+1
    for j = 1:2*r+1
        if (i-r-1)^2+(j-r-1)^2 <= r*r
            mask(i,j) = 1;
        end
    end
end
fa = imfilter(double(a),mask);
faa = fa.*a;
candidates = find(faa == 1);
candidates_val = auto_corr(candidates);

cur_max = 0;
dy = 0;
dx = 0;
offset = (size(auto_corr)+1)/2;
for i = 1 : length(candidates)
    if (candidates_val(i) >= cur_max)  
        [dy, dx] = ind2sub(size(auto_corr), candidates(i)); 
        dy = dy - offset(1);
        dx = dx - offset(2);
        cur_max = candidates_val(i); 
    end
end
dx = ceil(dx);
dy = ceil(dy);
c = estAttenuation(I, dx, dy);
