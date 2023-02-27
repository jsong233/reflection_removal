function k = genKernel(dx,dy,c)
%% genKernel(dx,dy,c) generate the kernel of shift (dx dy) and 
%  attenuation ratio c

center = [abs(dy),abs(dx)]+1;

k = zeros(2*center(1)-1,2*center(2)-1);
k(center(1),center(2)) = 1;
k(center(1)+dy,center(2)+dx) = c;
