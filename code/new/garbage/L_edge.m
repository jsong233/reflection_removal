function [LeY] = L_edge(Y)
Laplacian = [0 -1 0; -1 4 -1; 0 -1 0];
LeY = imfilter(Y,Laplacian,'symmetric','conv');
end

