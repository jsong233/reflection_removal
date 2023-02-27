%inverse stationary wavelet transform for 2D image
%s---signal
%order---spline order: 1 or 3
%level---decomposition level
function rs=iswt2(ds,order,level);

%get the dimension
[height, width, arrayno] = size(ds);

%get the mask
allmask = splmask(order,1,'U');
[maskno, masklen] = size(allmask);

%space allocation
%|  |    |   |...| ......|   |    |...|
%app(le)  det(le)             det(1)
%ie, first array is approx part, followed by detail part on level le,
%then detail part on level (le-1) all the way up to detail on level 1
for le=level:-1:1
    %approx part on current level
    if le==level
        att=reshape(ds(:,:,1),height,width);
    else
        att=rs;
    end
    %mask with padding zeros
    allmask = splmask(order,le,'U');
    allmask = allmask(:,end:-1:1); %flip the masks
    [maskno, masklen] = size(allmask);
    %reconstruction one level up
    %storage of wavelets on current level starts from windex+2
    windex = (level-le)*(maskno^2-1);
    rs = zeros(size(att));
    for fi=1:maskno
        %first reconstruct on Y
        tmpy = zeros(size(rs'));
        for fj=1:maskno
            if(fi==1)&&(fj==1)
                tmpim = att;
            else
                tmpim = reshape(ds(:,:,windex+maskno*(fi-1)+fj),height,width);
            end
            if mod(fj,2)==0
                es = MirrorExtY(tmpim',(masklen-1)/2,2);
            else
                es = MirrorExtY(tmpim',(masklen-1)/2,1);
            end
            tmpy = tmpy+conv2(es,allmask(fj,:)','valid');            
        end
        %then reconstruct on X
        if mod(fi,2)==0
            es = MirrorExtY(tmpy',(masklen-1)/2,2);
        else
            es = MirrorExtY(tmpy',(masklen-1)/2,1);
        end
        rs = rs+conv2(es,allmask(fi,:)','valid');
    end
    %rs=rs/4;
end
%================================================%