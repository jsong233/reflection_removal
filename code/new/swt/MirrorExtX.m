%mirror extension for a signal or image on X axis
function MX=MirrorExtX(M,sn,lr);

[height, width] = size(M);

if height == 1 %1D signal  
    if(lr==1)%symmetric extension
    MX=[M(sn:-1:1) M M(end:-1:end+1-sn)];
    end 
    if(lr==2)%asymmetric extension
    MX=[-M(sn:-1:1) M -M(end:-1:end+1-sn)];
    end    
else %2D image
    if(lr==1)%symmetric extension
    MX=[M(:,sn:-1:1) M M(:,end:-1:end+1-sn)];
    end    
    if(lr==2)%asymmetric extension
    MX=[-M(:,sn:-1:1) M -M(:,end:-1:end+1-sn)];
    end  
end