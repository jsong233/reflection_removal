%stationary frame transform for 2D image
%im---image
%order---spline order: 1 or 3
%level---decomposition level
function ds=swt2(im,order,level)

[height, width] = size(im);

%get the mask
allmask = splmask(order,1,'U');
[maskno, masklen] = size(allmask);

%space allocation
%|  |   |   |...| ......|   |    |...|
%app    det(le)           det(1)
%ie, first array is approx part, followed by detail part on level le,
%then detail part on level (le-1) all the way up to detail on level 1
wavecfs = zeros(height, width, 1+level*(maskno^2-1));
wavecfs(:,:,1) = im;
for le=1:level
    %mask with padding zeros
    allmask = splmask(order,le,'U');
    [maskno, masklen] = size(allmask);
    %filtering with filter on X, then with filter on Y
    app = reshape(wavecfs(:,:,1),height,width);
    windex = (level-le)*(maskno^2-1);
    for fi=1:maskno
        es = MirrorExtY(app,(masklen-1)/2,1);
%         allmask(fi,:)'
        tmpx = conv2(es,allmask(fi,:)','valid'); %on X
        %taking transpose is for filtering on Y
        es = MirrorExtY(tmpx',(masklen-1)/2,1); 
        for fj=1:maskno            
            tmpxy = conv2(es,allmask(fj,:)','valid'); %on Y
%             allmask(fj,:)'
            if(fi==1)&&(fj==1)
                wavecfs(:,:,1)=tmpxy'; %transpose again
            else
                wavecfs(:,:,windex+maskno*(fi-1)+fj)=tmpxy';                
            end
        end
    end
end

ds=wavecfs;
%=========%