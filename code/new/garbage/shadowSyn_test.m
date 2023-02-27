clear;
close all;

% 0.4 0.6 default test
% same initialization as shadowSyn_W
addpath('swt')
addpath('experiments')
test_case = 'shadowSyn.png';

I = im2double((imread(test_case)));
figure,imshow(I); % figure 1

%[dx1,dy1,c1] = estKernel(I); %[24 0 0.8946]
%[dx2,dy2,c2] = kernel(I); %[24 0 0.8513]
dx = -20; dy = -5; c = 0.6; 
configs.dx = dx;
configs.dy = dy;
configs.c = c;
[h,w,ch] = size(I);

k = genKernel(dx,dy,c);
figure,imshow(k) % figure 2
configs.h = h;
configs.w = w;
configs.ch = ch;
configs.num_px = h*w;

initT = ['initT_W_',test_case]; 
initR = ['initR_W_',test_case];

T0 = im2double(imread(initT));
R0 = im2double(imread(initR));
figure,imshow(T0); % figure 3
figure,imshow(R0); % figure 4

opts.mu = 1e5; opts2.mu = 1e5;
opts.delta = 1; opts2.delta = 1;
opts.nIter = 20; opts2.nIter = 20;
opts.nIterCG = 40; opts2.nIterCG = 40;

lambda = 0.01;
beta = 0.01;
gamma = 100;
wavelet_order = 1;
wavelet_level =  1;
level_response = 1+wavelet_level*((wavelet_order+2)^2-1);
Laplacian = [0 1 0;1 -4 1;0 1 0];


P1 = @(x) x(1:h, :);
P2 = @(x) x(h+1:end, :);

weight = genWeight(I,dx,dy);
figure,imshow(weight) % figure 5

A = @(x)  0.4 * P1(x) + 0.6 * imfilter(P2(x),k,'symmetric','conv');
AT = @(x) [0.4 * x; 0.6 * imfilter(x,rot90(k,2),'symmetric','conv')];

Z = (lambda+beta)*repmat(weight,1,1,level_response);
Z1 = (lambda+beta)*repmat((1-weight),1,1,level_response); Z(:,:,1) = 0;Z1(:,:,1) = 0;

M = @(x) cat(3,Z.*swt2(P1(x), wavelet_order, wavelet_level),Z1.*swt2(P2(x), wavelet_order, wavelet_level));
MT = @(x) [iswt2(Z.*x(:,:,1:level_response), wavelet_order, wavelet_level);iswt2(Z1.*x(:,:,1+level_response:end), wavelet_order, wavelet_level)];
D = @(x) gamma * [imfilter(x,[-1,1],'symmetric','conv');imfilter(x,[-1;1],'symmetric','conv')];
DT = @(x) gamma * (imfilter(P1(x),[1,-1],'symmetric','conv')+imfilter(P2(x),[1;-1],'symmetric','conv'));

eps = 1e-03; MAX_ITER = 20;
T_u = zeros(size(T0)); T = zeros(size(T0));
R_u = zeros(size(R0)); R = zeros(size(T0));

for i = 1:ch
    g = I(:,:,i);
    opts.u0 = [T0(:,:,i);R0(:,:,i)];
    u = SplitBregAnalysis(g, A, AT, M, MT, opts);   
    T_u(:,:,i) = P1(u); R_u(:,:,i) = P2(u);
    figure,imshow(P1(u));title('T');
    figure,imshow(P2(u));title('R');
    T(:,:,i) = P1(u);
    R(:,:,i) = P2(u);
end

figure,imshow(I);title('I');
figure,imshow(T);title('T');
figure,imshow(R);title('R');
 
maxI = max(max(max(I)));
minI = min(min(min(I)));
maxT = max(max(max(T)));
minT = min(min(min(T)));

T_r = minI + (T-minT)*(maxI-minI)/(maxT-minT);
figure,imshow(T_r);title('T_r'); 

curpath = pwd;
savepath = [curpath,'\results','\',mfilename];
if ~exist(savepath)
    mkdir(savepath)
end

Tpath = [savepath,'\','T_w_',test_case];
Rpath = [savepath,'\','R_w_',test_case];

imwrite(T_r,Tpath);
imwrite(R,Rpath);

d = T_r - T0;
norm(d(:),'fro')