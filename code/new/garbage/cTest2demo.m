% kernel和estKernel都未去除强边
addpath('experiments');

T = imread('cTest2_T.jpg'); T = rgb2gray(T);
R = imread('cTest_R.jpg'); R = rgb2gray(R);
k = [-20 -20 0.8]; dx = k(1); dy = k(2); c = k(3);
Y = T + R + c * imtranslate(R,[dx,dy]);
imwrite(Y, 'cTest2demo.jpg');
I = imread('cTest2demo.jpg');
cns = corner(I);
imshow(I); hold on;
plot(cns(:,1),cns(:,2),'r*');

[dx2, dy2, c2] =  kernel(I); % [5,2,0.8366];
[c,score,w,attn] = attenuation(I,-20,-20); % 0.4355