function x = zx(y)
%Find zero crossings (zx)

zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %find zx indices 
x = zci(y);  

end