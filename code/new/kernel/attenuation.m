function [c,score,w,attn] = attenuation(I, dx, dy)
% estimate the attenuation c of double reflection with shift [dx dy]
if size(I,3)~=1
    I = rgb2gray(I);
end
cns = corner(I);
hw = ceil(min(abs(dx),abs(dy))/2);
ave = fspecial('average',[3,3]);
psize = 2*hw + 1; center = (psize^2+1)/2; eps = 0.01;
ind = vec(1:psize^2); beta = exp(-abs(ind-center)/psize^2);
% m = [];
score = zeros(size(cns,1),1); 
attn = score*0;
w = score*0;
for j = 1 : size(cns,1)
    p1 = get_patch(I, cns(j,1), cns(j,2), hw);
    p2 = get_patch(I, cns(j,1) + dx, cns(j,2) + dy, hw);
    p1_ave = imfilter(p1,ave); p2_ave = imfilter(p2,ave);
    if ~isempty(p1) && ~isempty(p2)
        p1 = double(p1(:)); p2 = double(p2(:));
        p1_ave = double(p1_ave(:)); p2_ave = double(p2_ave(:));
        np1 = (sum(p1.^2))^0.5; np2 = (sum(p2.^2))^0.5;
        %p1 = p1 - mean(p1); p2 = p2 - mean(p2);
        if (np1 == 0) || (np2 == 0)
            score(j) = 0;
        else
            score(j) = sum(p1.*p2)/np1/np2;
        end
        attn(j) = sum(beta.*((p2_ave+eps)./(p1_ave+eps))) / sum(beta);
        if (attn(j) < 1) && ( attn(j) > 0)
            w(j) = exp(10*(score(j)-1)); % w(j) = exp(-score(j)/(2*0.2^2));
        end
    end
end
c = sum(w.*attn)/sum(w);
 

function p = get_patch(I, x, y, hw)
if (x>hw)&&(x<size(I,2)-hw)&&(y>hw)&&(y<size(I,1)-hw)
    p=I(y-hw:y+hw, x-hw:x+hw);
else
    p=[];
end


