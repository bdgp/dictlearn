function [up,down] = cutImg(d,position,direction,cutWidth,width,height)
I = reshape(d,height,width);
mask = zeros(height,width);
up = [];
down = [];
if strmatch(direction,'horizontal')
    trans = 1:(-1/(2*cutWidth)):0;
    trans = trans(1:(2*cutWidth))';
    lineMask1 = [ones(position-cutWidth,1);trans;zeros(height- ...
                                                      position - cutWidth,1)];
    mask = repmat(lineMask1,1,width);
    up = I.*mask;
    
    lineMask2 = lineMask1(end:-1:1);
    mask = repmat(lineMask2,1,width);
    down = I.*mask;
    up = up(:);
    down = down(:);
end