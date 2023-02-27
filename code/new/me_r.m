% our method
% reflection removal
% output T_r_XX and R_r_XX

clear;
close all;

addpath('swt')
addpath('kernel')
addpath('results')
addpath('figures')
addpath(['results','\',mfilename,'_init']);

test_case = 'me.jpg'; dx = 9; dy = -10; c = 0.8;
% test_case = 'book.jpg'; dx = -4; dy = -10; c = 0.65;
% if we change to the following test_case 'XX.jpg.png'
% we should copy this file and change the filename to 'XX_r.m'

% test_case = 'shadowSyn.png'; dx = -20; dy = -5; c = 0.6; 
% test_case = 'build.png'; [dx,dy,c] = estKernel(I);
% test_case = 'car.jpg'; dx = 10; dy = 6; c = 0.6;
% test_case = 'Flower.jpg'; [dx,dy,c] = estKernel(I);
% test_case = 'Lake.jpg'; dx = -32; dy = -8; c = 0.8;

I = im2double((imread(test_case)));
figure,imshow(I); 

configs.dx = dx;
configs.dy = dy;
configs.c = c;
[h,w,ch] = size(I);

k = genKernel(dx,dy,c);
figure,imshow(k) 
configs.h = h;
configs.w = w;
configs.ch = ch;
configs.num_px = h*w;

initT = ['initT_r_',test_case]; 
initR = ['initR_r_',test_case];

T0 = im2double(imread(initT));
R0 = im2double(imread(initR));
figure;imshow(T0); 
figure;imshow(R0); 

opts1.mu = 1e5; opts2.mu = 1e5;
opts1.delta = 1; opts2.delta = 1; 
opts1.nIter = 20; opts2.nIter = 20;
opts1.nIterCG = 40; opts2.nIterCG = 40;
MAX_ITER = 10;

d = 0.3; % d = 0.5;
lambda = 0.01;
beta = 0.01;
gamma = 0.001/(h*w); % control how much details in T we want
wavelet_order = 1;
wavelet_level =  1;
level_response = 1+wavelet_level*((wavelet_order+2)^2-1);

A = @(x)  (1-d) * imfilter(x,k,'symmetric','conv');
AT = @(x) (1-d) * imfilter(x,rot90(k,2),'symmetric','conv');
B = @(x) d * x;

weight = genWeight(I,dx,dy);
figure,imshow(weight);

Z = beta * repmat((1-weight),1,1,level_response); 
Z1 = lambda * repmat(weight,1,1,level_response);
Z(:,:,1) = 0; Z1(:,:,1) = 0;

eps = 1e-3;
T_r = zeros(size(T0));
R_r = zeros(size(T0));
for i = 1:ch
    Y = I(:,:,i);
    for j = 1:MAX_ITER
        if j == 1
            opts1.u0 = R0(:,:,i); T = T0(:,:,i);
        else
            opts1.u0 = R;
        end
        g = Y - B(T);
        W = @(x) cat(3,Z.*swt2(x, wavelet_order, wavelet_level),gamma*imageLike(T,x)*ones(h,w));
        WT = @(x) iswt2(Z.*x(:,:,1:level_response), wavelet_order, wavelet_level) ...
            + (gamma*h*w) * x(:,:,1+level_response:end).*T/norm(T,'fro')/norm(g/(1-d),'fro');
        r = SplitBregAnalysis(g, A, AT, W, WT, opts1); R = r;
    
        l = Y - A(R);opts2.u0 = T; 
        M = @(x) cat(3,Z1.*swt2(x, wavelet_order, wavelet_level),gamma*imageLike(R,x)*ones(h,w));
        MT = @(x) iswt2(Z1.*x(:,:,1:level_response), wavelet_order, wavelet_level) ...
            + (gamma*h*w) * x(:,:,1+level_response:end).*R/norm(R,'fro')/norm(l/d,'fro');
        t = SplitBregAnalysis(l, B, B, M, MT, opts2);
       
        figure; subplot(121);imshow(t);title('T'); subplot(122);imshow(r);title('R');
        if norm(T - t,2)/norm(T,2) < eps
            break;
        else
            T = t;
        end
    end
    T_r(:,:,i) = T; R_r(:,:,i) = R;
end

figure;imshow(I);title('I');
figure;imshow(T_r);title('T');
figure;imshow(R_r);title('R');
 
maxI = max(max(max(I)));
minI = min(min(min(I)));
maxT = max(max(max(T_r)));
minT = min(min(min(T_r)));

T = minI + (T_r-minT)*(maxI-minI)/(maxT-minT);
figure;imshow(T);title('T_r'); 

curpath = pwd;
savepath = [curpath,'\results','\',mfilename];
if ~exist(savepath)
    mkdir(savepath)
end

Tpath = [savepath,'\','T_r_',test_case];
Rpath = [savepath,'\','R_r_',test_case];

imwrite(T,Tpath);
imwrite(R_r,Rpath);
