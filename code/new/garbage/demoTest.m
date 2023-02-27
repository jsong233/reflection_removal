% kernel和estKernel都未去除强边
addpath('experiments');
addpath('results\demo_real_deghost');
I = imread('demo.jpg');
%I = im2double(imread('A9RF22E.png')); I = rgb2gray(I);

[dx1,dy1,c1] = estKernel(I); % [0,0,NAN]
[dx2,dy2,c2] = kernel(I);  % [71,0,0.8107] 
% 用Kernel取出来的还是背景的pattern 因为demo中T没有角点 所以估算的c都不大