function [dx, dy, c] =  kernel(I) 
%% estimate the shift [dx dy] and the attenuation c
dim = ndims(I);
if dim == 3
    I = rgb2gray(I);
end
% Laplacian = [0 -1 0; -1 4 -1; 0 -1 0]; 
% resp = imfilter(I, Laplacian);

hy = fspecial('sobel');
hx = hy';
Iy = imfilter(double(I), hy, 'replicate');
Ix = imfilter(double(I), hx, 'replicate');
gI = sqrt(Ix.^2+Iy.^2);

cen = 3;
[~,ce] = kmeans(gI(:),cen);
[cv,~] = sort(ce,'descend');
%gI(gI > cv(2)) = 0;
resp = gI;
%resp(abs(resp)>0.2) = 0; 
                                               
auto_corr = xcorr2(abs(resp), abs(resp)); % 对edge map求自相关矩阵
                                          % 若移动越界，则其余部分补为0
                                          % c(i,j):原图与原图下移i右移j的卷积
                                          
bdry = ceil(min(size(I))/2);         % eliminate the autocorrelation in boundary
auto_corr = auto_corr(1+bdry:end-bdry, 1+bdry:end-bdry);

max_1 = ordfilt2(auto_corr, 25, true(5));

center = size(I) - bdry;
auto_corr(center(1), center(2)) = 0;

a = (auto_corr == max_1);
candidates = find(a == 1);
candidates_val = auto_corr(candidates);

offset = (size(auto_corr)+1)/2;
cur_max = 0; dx = 0; dy = 0;

for i = 1 : length(candidates)        
    [cur_y, cur_x] = ind2sub(size(auto_corr), candidates(i)); 
    cur_y = cur_y - offset(1); cur_x = cur_x - offset(2);
    cur_c = attenuation(I,cur_x,cur_y);
    if (candidates_val(i) >= cur_max) && (cur_c < 0.9)
        dy = cur_y; dx = cur_x; c = cur_c;
        cur_max = candidates_val(i); % current maximum
    end
end
end
