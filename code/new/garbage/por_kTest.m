clear;
close all;
% clc;
%
%  min 0.5*||I-T-R*k||_2^2 + lambda*||WT||_1 + beta*||WR||_1
%  T,R
addpath('swt')
addpath('experiments')
test_case = 'portrait.png';

I = im2double((imread(test_case)));

[dx1,dy1,c1] = estKernel(I); %[12 5 0.8526]
[dx2,dy2,c2] = kernel(I); %[40 0 0.4332]
