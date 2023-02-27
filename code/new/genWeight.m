function weight = genWeight(I,dx,dy)
    [h,w,ch] = size(I);
    if ch == 3
        I = rgb2gray(I);
    end
    weight = 0.5*ones(h,w);
    
    gI = double(edge(I,'Canny'));
    gau = fspecial('gaussian',3,0.5);
    gI = imfilter(gI,gau);
    gI = imfilter(gI,gau);
    gI = imfilter(gI,gau);  
%     gI = imfilter(gI,gau);  

    patSize = 8;
    center = floor((1+patSize)/2);
    abs_dy = abs(dy);
    abs_dx = abs(dx);
    padI = padarray(gI,[2*abs_dy+center,2*abs_dx+center],'both');
    for m = 2*abs_dy+center + 1: 2*abs_dy+center + h
        for n = 2*abs_dx+center + 1: 2*abs_dx+center + w
           Igp = getPatchP(padI,m,n,[patSize,patSize]);
           Igp_shift_r1 = getPatchP(padI,m+dy,n+dx,[patSize,patSize]);
           Igp_shift_r2 = getPatchP(padI,m+2*dy,n+2*dx,[patSize,patSize]);
           Igp_shift_l1 = getPatchP(padI,m-dy,n-dx,[patSize,patSize]);
           Igp_shift_l2 = getPatchP(padI,m-2*dy,n-2*dx,[patSize,patSize]);
           r1 = 0; r2 = 0; l1 = 0; l2 = 0;
           if padI(m,n) > 5e-3
               if sum(sum(Igp_shift_r1)) > 5e-3
                   r1 =  (sum(sum(Igp.*Igp_shift_r1)))/(norm(Igp,'fro'))/(norm(Igp_shift_r1,'fro'));  
               end
               if sum(sum(Igp_shift_r2)) > 5e-3
                   r2 =  (sum(sum(Igp.*Igp_shift_r2)))/(norm(Igp,'fro'))/(norm(Igp_shift_r2,'fro'));  
               end
               if sum(sum(Igp_shift_l1)) > 5e-3
                   l1 =  (sum(sum(Igp.*Igp_shift_l1)))/(norm(Igp,'fro'))/(norm(Igp_shift_l1,'fro'));  
               end
               if sum(sum(Igp_shift_l2)) > 5e-3
                   l2 =  (sum(sum(Igp.*Igp_shift_l2)))/(norm(Igp,'fro'))/(norm(Igp_shift_l2,'fro'));  
               end

               v = sort([r1,r2,l1,l2],'descend');
               weight(m - 2*abs_dy - center,n - 2*abs_dx - center) = abs(v(1) - v(2));  

           end
        end
    end

    weight(weight<0.5) = 0.5;
    weight(weight>0.5) = 1;

    hy = fspecial('sobel');
    hx = hy';
    Iy = imfilter(double(I), hy, 'replicate');
    Ix = imfilter(double(I), hx, 'replicate');
    ggI = sqrt(Ix.^2+Iy.^2);
    ggI = imfilter(ggI,gau);
    ggI = imfilter(ggI,gau);

    cen = 3;
    [v,ce] = kmeans(ggI(:),cen);
    [cv,idx] = sort(ce,'descend');

    ggI(ggI<cv(2)) = 0;
    weight(ggI~=0 & weight~=1) = 0;
