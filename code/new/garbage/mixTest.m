addpath('figures');
I = imread('A9RF22E.png');
T = imread('T_A9RF22E.png');
R = imread('R_A9RF22E.png');
figure(1); 
subplot(2,3,1); imshow(I);
subplot(2,3,2); imshow(T);
subplot(2,3,3); imshow(R);

k = [20 1 0.6];
Y1 = mix(T,R,k);
Y2 = Y1/2;
subplot(2,3,4); imshow(rgb2gray(I));
subplot(2,3,5); imshow(Y1);
subplot(2,3,6); imshow(Y2);
