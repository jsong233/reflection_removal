addpath('figures');
T = imread('T.jpg'); T = rgb2gray(T);
R = imread('R.jpg'); R = rgb2gray(R);
k = [-5 -5 0.8]; dx = k(1); dy = k(2); c = k(3);
Y = mix(T,R,k);
imwrite(Y, 'demo.jpg');

I = imread('demo.jpg');
imshow(I);