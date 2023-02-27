function [r] = imageLike(I1,I2)
I1 = double(I1(:)); I2 = double(I2(:));
nI1 = (sum(I1.^2))^0.5; nI2 = (sum(I2.^2))^0.5;
r = sum(I1.*I2)/nI1/nI2;
end

