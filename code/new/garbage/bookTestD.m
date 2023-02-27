% 寻找最佳不透明度混合模式
clear;
close all;
% clc;
%
%  min 0.5*||I-T-R*k||_2^2 + lambda*||WT||_1 + beta*||WR||_1
%  T,R
addpath('swt')
addpath('experiments')
test_case = 'book.jpg';

I = im2double((imread(test_case)));
figure,imshow(I); % figure 1
[h,w,ch] = size(I);

% [dx,dy,c] = kernel(I);
% [dx,dy,c] = estKernel(I); %[23 4 0.649]
dx = -5; dy = -17; c = 0.65;

save([test_case(1:end-4),'.mat'],'dx','dy','c');
% load([test_case(1:end-4),'.mat']);
configs.dx = dx;
configs.dy = dy;
configs.c = c;

k = genKernel(dx,dy,c);
figure,imshow(k) % figure 2
configs.h = h;
configs.w = w;
configs.ch = ch;
configs.num_px = h*w;

initT = ['initT_',test_case];
initR = ['initR_',test_case];

T0 = im2double(imread(initT));
R0 = im2double(imread(initR));
figure,imshow(T0);title('T0');
figure,imshow(R0);title('R0');

opts.mu = 1e5;
opts.delta = 1;
opts.nIter = 20;
opts.nIterCG = 40;

lambda = 0.01;
beta = 0.01;
wavelet_order = 1;
wavelet_level =  1;
level_response = 1+wavelet_level*((wavelet_order+2)^2-1);
Laplacian = [0 1 0;1 -4 1;0 1 0];


P1 = @(x) x(1:h, :);
P2 = @(x) x(h+1:end, :);

weight = genWeight(I,dx,dy);
figure,imshow(weight);

Z = (lambda+beta)*repmat(weight,1,1,level_response);
Z1 = (lambda+beta)*repmat((1-weight),1,1,level_response); Z(:,:,1) = 0;Z1(:,:,1) = 0;

W = @(x) cat(3,Z.*swt2(P1(x), wavelet_order, wavelet_level),Z1.*swt2(P2(x), wavelet_order, wavelet_level));
WT = @(x) [iswt2(Z.*x(:,:,1:level_response), wavelet_order, wavelet_level);iswt2(Z1.*x(:,:,1+level_response:end), wavelet_order, wavelet_level)];

opacity = 0.1:0.1:0.9;
T = zeros(size(T0));
R = zeros(size(R0));
for i = 1:ch
    g = I(:,:,i);
    opts.u0 = [T(:,:,i);R(:,:,i)];
    for j = 1:length(opacity)
        d = opacity(j);
        A = @(x)  d * P1(x) + (1-d) * imfilter(P2(x),k,'symmetric','conv');
        AT = @(x) [d * P1(x); (1-d) * imfilter(x,rot90(k,2),'symmetric','conv')];
        u = SplitBregAnalysis(g, A, AT, W, WT, opts);
        r = norm(imfilter(P1(u),Laplacian,'symmetric','conv'),'fro') + norm(imfilter(P2(u),Laplacian,'symmetric','conv'),'fro');
        if j == 1
            tmp = r;
        end
        if r <= tmp
            tmp = r; D = d;
            T(:,:,i) = P1(u); R(:,:,i) = P2(u);
            figure,imshow(P1(u));title('T');
            figure,imshow(P2(u));title('R');
        end
    end
end
figure,imshow(I);title('I')
figure,imshow(R);title('R')
figure,imshow(T);title('T')

a = T; b = R;

T = a;
maxI = max(max(max(I)));
minI = min(min(min(I)));
maxT = max(max(max(T)));
minT = min(min(min(T)));

T = minI + (T-minT)*(maxI-minI)/(maxT-minT);
figure,imshow(T) % figure 15

curpath = pwd;
savepath = [curpath,'\results','\',mfilename];
if ~exist(savepath)
    mkdir(savepath)
end

Tpath = [savepath,'\','T_d_',test_case];
Rpath = [savepath,'\','R_d_',test_case];
weightPath =  [savepath,'\','weight_d_',test_case(1:end-4),'.mat'];
imwrite(T_r,Tpath);
imwrite(R_r,Rpath);
save(weightPath,'weight','lambda','beta')

d = T - T0;
norm(d(:),'fro')
