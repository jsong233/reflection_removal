% 寻找最佳不透明度混合模式
clear;
close all;

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
figure;imshow(T0); % figure 3
figure;imshow(R0); % figure 4

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
D = zeros(1,3);
T = zeros(size(T0));
R = zeros(size(T0));
for i = 1:ch
    g = I(:,:,i);
    opts.u0 = [T0(:,:,i);R0(:,:,i)];
    for j = 1:length(opacity)
        d = opacity(j);
        A = @(x)  d * P1(x) + (1-d) * imfilter(P2(x),k,'symmetric','conv');
        AT = @(x) [d * P1(x); (1-d) * imfilter(x,rot90(k,2),'symmetric','conv')];
        u = SplitBregAnalysis(g, A, AT, W, WT, opts);
        eT = Sobel(P1(u)); eR = Sobel(P2(u));
        r = imageLike(eT,eR);
        if j == 1
            tmp = r;
        end
        if r <= tmp
            tmp = r; D(i) = d;
            T(:,:,i) = P1(u); R(:,:,i) = P2(u);
        end
        figure;subplot(121);imshow(P1(u));xlabel('T');
        subplot(122);imshow(P2(u));xlabel('R');
        suptitle(['i = ',num2str(i),' d = ',num2str(d),' r = ',num2str(r)]);
    end
end

figure;imshow(I);title('I');
figure;imshow(T);title('T');
figure;imshow(R);title('R');
 
maxI = max(max(max(I)));
minI = min(min(min(I)));
maxT = max(max(max(T)));
minT = min(min(min(T)));

T_r = minI + (T-minT)*(maxI-minI)/(maxT-minT);
figure;imshow(T_r);title('T_r'); 

curpath = pwd;
savepath = [curpath,'\results','\',mfilename];
if ~exist(savepath)
    mkdir(savepath)
end

Tpath = [savepath,'\','T_D_',test_case];
Rpath = [savepath,'\','R_D_',test_case];

imwrite(T_r,Tpath);
imwrite(R,Rpath);
