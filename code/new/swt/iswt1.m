%inverse stationary wavelet transform for 1D signnal
%s---signal
%order---spline order: 1 or 3
%level---decomposition level
function rs=iswt1(ds,order,level);

%convert from column-wise storing to row-wise storing
ds = ds'; 

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
        att=ds(1,:);
    else
        att=rs;
    end
    %mask with padding zeros
    allmask = splmask(order,le,'U');
    allmask = allmask(:,end:-1:1);
    [maskno, masklen] = size(allmask);
    %reconstruction one level up
    rs = zeros(1,length(att));
    es = MirrorExtX(att,(masklen-1)/2,1);
    rs = convn(es,allmask(1,:),'valid');
    windex = (level-le)*(maskno-1);
    for i=2:maskno
        %extend the coeffs
        coeffs = ds(windex+i,:);
        if mod(i,2)==0
            es = MirrorExtX(coeffs,(masklen-1)/2,2);%odd mask 
        else
            es = MirrorExtX(coeffs,(masklen-1)/2,1);%even mask
        end            
        %inverse transform
        rs = rs+convn(es,allmask(i,:),'valid');
    end
    %rs = rs/2;
end

%column-wise storing
rs = rs';
%================================================%