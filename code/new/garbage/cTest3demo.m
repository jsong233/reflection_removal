R1 = imread('cTest3_R.jpg'); R2 = imread('R.jpg');
Y = R1 + 0.8 * imtranslate(R2,[-20,-20]);
imwrite(Y,'cTest3demo.jpg');

Y = imread('cTest3demo.jpg');
Y = rgb2gray(Y);
cns = corner(Y); dx = -20; dy = -20; % 实际 c = 0.5

max_text = {};
for i = 1:length(cns)
    max_text = [max_text,num2str(i)];
end

figure(1);
imshow(Y); hold on; plot(cns(:,1),cns(:,2),'r*');
hold on; 
text(cns(:,1)+0.05,cns(:,2)+0.05,max_text,'Color','yellow');

[c1,score1,w1,attn1] = estAttenuation(Y,dx,dy); % 0.0679 偏小
[c2,score2,w2,attn2] = attenuation(Y,dx,dy); % 0.8081较为准确


j= 2;hw = ceil(min(abs(dx),abs(dy))/2);
ave = fspecial('average',[3,3]);
psize = 2*hw + 1; center = (psize^2+1)/2; eps = 0.01;
ind = vec(1:psize^2); beta = exp(-abs(ind-center)/psize^2);
p1 = get_patch(Y, cns(j,1), cns(j,2), hw);
p2 = get_patch(Y, cns(j,1) + dx, cns(j,2) + dy, hw);
p1_ave = imfilter(p1,ave); p2_ave = imfilter(p2,ave);
figure(2); subplot(2,2,1); imshow(p1); subplot(2,2,2); imshow(p2);
subplot(2,2,3); imshow(p1_ave); subplot(2,2,4); imshow(p2_ave);

p1 = double(p1(:)); p2 = double(p2(:));
p1_ave = double(p1_ave(:)); p2_ave = double(p2_ave(:));
np1 = (sum(p1.^2))^0.5; np2 = (sum(p2.^2))^0.5;
sum(beta.*((p2_ave+eps)./(p1_ave+eps))) / sum(beta)