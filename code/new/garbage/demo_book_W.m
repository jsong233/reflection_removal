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

%[dx,dy,c] = estKernel(I);
dx = -5; dy = -17; c = 0.65;
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

T0 = 2*im2double(imread(initT));
R0 = 2*im2double(imread(initR));
figure,imshow(T0); % figure 3
figure,imshow(R0); % figure 4

opts.mu = 1e5; opts2.mu = 1e5;
opts.delta = 1; opts2.delta = 1;
opts.nIter = 20; opts2.nIter = 20;
opts.nIterCG = 40; opts2.nIterCG = 40;

lambda = 0.01;
beta = 0.01;
gamma = 10;
wavelet_order = 1;
wavelet_level =  1;
level_response = 1+wavelet_level*((wavelet_order+2)^2-1);
Laplacian = [0 1 0;1 -4 1;0 1 0];


P1 = @(x) x(1:h, :);
P2 = @(x) x(h+1:end, :);

weight = genWeight(I,dx,dy);
figure,imshow(weight) % figure 5

Z = (lambda+beta)*repmat(weight,1,1,level_response);
Z1 = (lambda+beta)*repmat((1-weight),1,1,level_response); Z(:,:,1) = 0;Z1(:,:,1) = 0;

M = @(x) cat(3,Z.*swt2(P1(x), wavelet_order, wavelet_level),Z1.*swt2(P2(x), wavelet_order, wavelet_level));
MT = @(x) [iswt2(Z.*x(:,:,1:level_response), wavelet_order, wavelet_level);iswt2(Z1.*x(:,:,1+level_response:end), wavelet_order, wavelet_level)];
D = @(x) gamma * [imfilter(x,[-1,1],'symmetric','conv');imfilter(x,[-1;1],'symmetric','conv')];
DT = @(x) gamma * (imfilter(P1(x),[1,-1],'symmetric','conv')+imfilter(P2(x),[1;-1],'symmetric','conv'));

eps = 1e-05; MAX_ITER = 15;
W = ones(h,w,3);
T_u = zeros(size(T0)); T = zeros(size(T0));
R_u = zeros(size(R0)); R = zeros(size(T0));

for i = 1:ch
    g = I(:,:,i);
    opts.u0 = [T0(:,:,i);R0(:,:,i)]; opts2.u0 = ones(h,w);
    for j = 1:MAX_ITER
        A = @(x)  W(:,:,i) .* P1(x) + (2-W(:,:,i)).*imfilter(P2(x),k,'symmetric','conv');
        AT = @(x) [W(:,:,i) .* x; (2-W(:,:,i)) .* imfilter(x,rot90(k,2),'symmetric','conv')];
        u = SplitBregAnalysis(g, A, AT, M, MT, opts);   
        T_u(:,:,i) = P1(u); R_u(:,:,i) = P2(u);
        figure,imshow(P1(u));title('T');
        figure,imshow(P2(u));title('R');
        figure,imshow(W(:,:,i));title('W');
        l = g - imfilter(R_u(:,:,i),k,'symmetric','conv');
        B = @(x) (T_u(:,:,i)-imfilter(R_u(:,:,i),k,'symmetric','conv')).*x;
        v = SplitBregAnalysis2(l, B, B, D, DT, opts2);
        r = norm(W(:,:,i) - v,2)/norm(W(:,:,i),2);
        if r < eps
            W(:,:,i) = v; break
        else
            W(:,:,i) = v;
        end
    end
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

Tpath = [savepath,'\','T_W_',test_case];
Rpath = [savepath,'\','R_W_',test_case];
Wpath = [savepath,'\','W_',test_case];
weightPath =  [savepath,'\','weight_W_',test_case(1:end-4),'.mat'];

imwrite(T_r,Tpath);
imwrite(R,Rpath);
imwrite(W,Wpath);
save(weightPath,'weight','lambda','beta','gamma')

d = T_r - T0;
norm(d(:),'fro')

%test
W1 = W(:,:,1);
W2 = W(:,:,2);
W3 = W(:,:,3);