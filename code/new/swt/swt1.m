%stationary wavelet transform for 1D signal
%s---signal
%order---spline order: 0, 1 or 3
%level---decomposition level
function ds=swt1(s,order,level)

s = s(:); %column vector
slen=length(s);

%get the mask
allmask = splmask(order,1,'U');
[maskno, masklen] = size(allmask);

%space allocation
%|  |   |   |...| ......|   |    |...|
%app    det(le)           det(1)
%ie, first array is approx part, followed by detail part on level le,
%then detail part on level (le-1) all the way up to detail on level 1
wavecfs = zeros(1+level*(maskno-1), slen);
wavecfs(1,:) = s;
for le=1:level
    %mask with padding zeros
    allmask = splmask(order,le,'U');
    [maskno, masklen] = size(allmask);
    %extend signal
    es = MirrorExtX(wavecfs(1,:),(masklen-1)/2,1);
    %decomposition part
    wavecfs(1,:) = convn(es,allmask(1,:),'valid');
    windex = (level-le)*(maskno-1);
    for i=2:maskno
        wavecfs(windex+i,:) = convn(es,allmask(i,:),'valid');
    end
end

ds=wavecfs'; %columnwise storing
%==================================%