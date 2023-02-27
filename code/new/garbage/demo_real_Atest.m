clear;
close all;
% clc;
%
%  min 0.5*||I-T-R*k||_2^2 + lambda*||WT||_1 + beta*||WR||_1
%  T,R
addpath('swt')
test_case = 'A9RF22E.png';

I = im2double((imread(test_case)));

[dx,dy,c] = estKernel(I); % [20,1,0.9251]
save([test_case(1:end-4),'.mat'],'dx','dy','c');
% load([test_case(1:end-4),'.mat']);
configs.dx = dx;
configs.dy = dy;
configs.c = c;
figure,imshow(I); % figure 1
[h,w,ch] = size(I);

k = genKernel(dx,dy,c);

figure,imshow(k) % figure 2
configs.h = h;
configs.w = w;
configs.ch = ch;
configs.num_px = h*w;

initT = ['initT_',test_case]; 
initR = ['initR_',test_case];

T = im2double(imread(initT));
R = im2double(imread(initR));
figure,imshow(T); % figure 3
figure,imshow(R); % figure 4

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

weight = 0.5*ones(h,w);
figure,imshow(weight) % figure 5

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
    figure,imshow(P1(r));title('T') % figure 6,8,10
    figure,imshow(P2(r));title('R') % figure 7,9,11
    T_r(:,:,i) = P1(r);   
    R_r(:,:,i) = P2(r);
end
figure,imshow(I);title('I') % figure 12
figure,imshow(R_r);title('R') % figure 13
figure,imshow(T_r);title('T') % figure 14

a = T_r; b = R_r;

T_r = a; 
maxI = max(max(max(I)));
minI = min(min(min(I)));
maxT = max(max(max(T_r)));
minT = min(min(min(T_r)));

T_r = minI + (T_r-minT)*(maxI-minI)/(maxT-minT);
figure,imshow(T_r) % figure 15

curpath = pwd;
savepath = [curpath,'\results','\',mfilename];
if ~exist(savepath)
    mkdir(savepath)
end

Tpath = [savepath,'\','T_Atest_',test_case];
Rpath = [savepath,'\','R_Atest_',test_case];
weightPath =  [savepath,'\','weight_Atest_',test_case(1:end-4),'.mat'];
imwrite(T_r,Tpath);
imwrite(R_r,Rpath);
save(weightPath,'weight','lambda','beta')

d = T_r - T;
norm(d(:),'fro')
