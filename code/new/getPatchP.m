function p = getPatchP(I, y, x, patSize)
    if size(patSize,2) == 1
        patSize = [patSize patSize];
    end
    [h, w] = size(I);
    indy = mod(patSize(1),2)~=0;
    indx = mod(patSize(2),2)~=0;
    hw = floor((1+patSize)/2);

    p=I(y-hw(1)+1:y+hw(1)-1*indy,x-hw(2)+1:x+hw(2)-1*indx);

end