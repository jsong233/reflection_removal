function p = get_patch(I, x, y, hw)
if (x>hw)&&(x<size(I,2)-hw)&&(y>hw)&&(y<size(I,1)-hw)
    p=I(y-hw:y+hw, x-hw:x+hw);
else
    p=[];
end