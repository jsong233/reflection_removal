% compare estKernel.m and kernel.m on synthetic image

% use kernel2 and estKernel2 cause for artificial images
% we do not eliminate strong edges

addpath('figures');
addpath('kernel');

T = imread('cTest_T.jpg'); 
R = imread('cTest_R.jpg'); 
% pattern from T: -22 -22
k = [10 -10 0.8]; dx = k(1); dy = k(2); c = k(3);
Y = mix(T,R,k); Y = rgb2gray(Y);
imwrite(Y, 'cTestdemo.jpg');
I = imread('cTestdemo.jpg');
imshow(I);

[dx1, dy1, c1] =  estKernel2(I); % [21,22,0.6883]
[dx2, dy2, c2] =  kernel2(I); % [10,-10,0.8671]
