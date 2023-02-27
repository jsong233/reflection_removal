clear;
close all;

addpath('experiments')

hy = fspecial('sobel');
hx = hy';

T = im2double(imread('synT5.png'));
R = im2double(imread('synR4.png'));
figure; subplot(121); imshow(T); subplot(122); imshow(R);

dx = -20; dy = -5; c = 0.6;
k = genKernel(dx,dy,c);
d = 0.3;

Y = d*T + (1-d)*imfilter(R,k,'symmetric','conv');
figure; imshow(Y);
imwrite(Y,'shadowSyn.png');
