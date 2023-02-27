% kernel和estKernel里面都去除了强边
addpath('swt')
test_case = 'A9RF22E.png';

I = im2double((imread(test_case)));

[dx1,dy1,c1] = estKernel(I); % [20,1,0.9251] 
[dx2,dy2,c2] = kernel(I); 
% 去除cv(1)以上：[29,2,0.7333]  attenuation(I,20,1) = 0.7147
% 去除0.2以上：kernel(I) = [20,1,0.7147]
% 去除cv(2)以上：kernel(I) = [20,1,0.7147]
