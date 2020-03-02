function R=aggregate(AP,blk)
    
    ndoms = length(blk)-1;
    n = length(AP);
    R = spalloc(ndoms,n,n);
    for i=1:ndoms
        R(i,blk(i):blk(i+1)-1) = 1;%/length(blk(i):blk(i+1)-1);
    end
    
end