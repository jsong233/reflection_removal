%hard or soft  thresholding
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ds1 = wavethresh(ds,t,string)

wave = ds(:,:,2:end);

if string=='hard' %hard thresholding    

    ds1 = wave .* (abs(wave) > t);

elseif string=='soft' %soft thresholding
    
    res = (abs(wave) - t);
    res = (res + abs(res))/2;
    ds1 = sign(wave).*res;

end

%keep approximation part
ds1 = cat(3, ds(:,:,1), ds1);
%==========================%