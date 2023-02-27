function [Y] = mix(T,R,k)
dx = k(1); dy = k(2); c = k(3);
Y = T + R + c * imtranslate(R,[dx,dy]);
end

