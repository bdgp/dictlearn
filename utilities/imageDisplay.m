function imageDisplay(pc,width,height,new,r)
    
    if ~exist('r')
        r = [0, max(pc)];
    end

    if new == 1
        figure('Position',[1436 684 794 366]);
    end;
    %axis ij;
    
    imagesc(reshape(pc,height,width));
    caxis(r);
    axis off;
    
    if new == 1
        colorbar;
    end
end