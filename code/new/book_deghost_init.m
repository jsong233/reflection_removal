% Yan's method
% initialization for all the XX_deghost.m
% output initT_XX and initR_XX

% we measure the ground-truth [dx,dy,c] in Photoshop for each input

clear;
close all;

addpath('swt')
addpath('kernel')
addpath('figures')
test_case = 'book.jpg'; dx = -4; dy = -10; c = 0.65;
% test_case = 'shadowSyn.png'; dx = -20; dy = -5; c = 0.6; 
% test_case = 'build.png'; [dx,dy,c] = estKernel(I);
% test_case = 'car.jpg'; dx = 10; dy = 6; c = 0.6;
% test_case = 'Flower.jpg'; [dx,dy,c] = estKernel(I);
% test_case = 'Lake.jpg'; dx = -32; dy = -8; c = 0.8;

I = im2double((imread(test_case))); 
figure;imshow(I); 
[h,w,ch] = size(I);

configs.dx = dx;
configs.dy = dy;
configs.c = c;

k = genKernel(dx,dy,c);
figure;imshow(k);
configs.h = h;
configs.w = w;
configs.ch = ch;
configs.num_px = h*w;

% ------------------------- first iter ------------------------------------
initT = I; initR = zeros(size(I));
T = initT; R = initR;

figure,imshow(T);
figure,imshow(R); 

opts.mu = 1e5; 
opts.delta = 1;
opts.nIter = 20;
opts.nIterCG = 40;

lambda = 0.01;
beta = 0.01;
wavelet_order = 1;
wavelet_level =  1;
level_response = 1+wavelet_level*((wavelet_order+2)^2-1);

P1 = @(x) x(1:h, :);
P2 = @(x) x(h+1:end, :);
A = @(x)  P1(x) + imfilter(P2(x),k,'symmetric','conv');
AT = @(x) [P1(x); imfilter(x,rot90(k,2),'symmetric','conv')];
Ar = @(x) imfilter(x,k,'symmetric','conv');

weight = 0.5 * ones(h,w);
figure,imshow(weight)

Z = (lambda+beta)*repmat(weight,1,1,level_response);
Z1 = (lambda+beta)*repmat((1-weight),1,1,level_response); Z(:,:,1) = 0;Z1(:,:,1) = 0;

W = @(x) cat(3,Z.*swt2(P1(x), wavelet_order, wavelet_level),Z1.*swt2(P2(x), wavelet_order, wavelet_level));
WT = @(x) [iswt2(Z.*x(:,:,1:level_response), wavelet_order, wavelet_level);iswt2(Z1.*x(:,:,1+level_response:end), wavelet_order, wavelet_level)];

T_r = zeros(size(T));
R_r = zeros(size(R));
for i = 1:ch
    g = I(:,:,i); 
    opts.u0 = [T(:,:,i);R(:,:,i)]; 
    [r, rec] = SplitBregAnalysis(g, A, AT, W, WT, opts);   
    figure,imshow(P1(r));title('T') 
    figure,imshow(P2(r));title('R') 
    T_r(:,:,i) = P1(r);   
    R_r(:,:,i) = P2(r);
end
figure,imshow(I);title('I')
figure,imshow(R_r);title('R') 
figure,imshow(T_r);title('T')

maxI = max(max(max(I)));
minI = min(min(min(I)));
maxT = max(max(max(T_r)));
minT = min(min(min(T_r)));

T_r = minI + (T_r-minT)*(maxI-minI)/(maxT-minT);
figure,imshow(T_r);

curpath = pwd;
savepath = [curpath,'\results','\',mfilename];
if ~exist(savepath)
    mkdir(savepath)
end

Tpath = [savepath,'\','initT_',test_case];
Rpath = [savepath,'\','initR_',test_case];
imwrite(T_r,Tpath);
imwrite(R_r,Rpath);

% -----------------------------------------------------------------------

initITER = 5;

for j = 1:initITER
    T = T_r;
    R = R_r;
    T_r = zeros(size(T));
    R_r = zeros(size(R));
    for i = 1:ch
        g = I(:,:,i); 
        opts.u0 = [T(:,:,i);R(:,:,i)]; 
        [r, rec] = SplitBregAnalysis(g, A, AT, W, WT, opts);   
        figure,imshow(P1(r));title('T') 
        figure,imshow(P2(r));title('R') 
        T_r(:,:,i) = P1(r);   
        R_r(:,:,i) = P2(r);
    end
    figure,imshow(I);title('I')
    figure,imshow(R_r);title('R') 
    figure,imshow(T_r);title('T')

    maxI = max(max(max(I)));
    minI = min(min(min(I)));
    maxT = max(max(max(T_r)));
    minT = min(min(min(T_r)));

    T_r = minI + (T_r-minT)*(maxI-minI)/(maxT-minT);
    figure,imshow(T_r);

    curpath = pwd;
    savepath = [curpath,'\results','\',mfilename];
    if ~exist(savepath)
        mkdir(savepath)
    end

    Tpath = [savepath,'\','initT_',test_case];
    Rpath = [savepath,'\','initR_',test_case];
    imwrite(T_r,Tpath);
    imwrite(R_r,Rpath);
end
