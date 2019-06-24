function [blockvar, sizearrblock] = blockave(quantity)

lenquant = length(quantity);
blockmax = floor(0.4*lenquant);
blockvar = zeros(floor(0.2*blockmax),1);
sizearrblock  = zeros(floor(0.2*blockmax),1);
cnt = 0;
for blocksize = 1:blockmax
    
    % Split entire length into "nchunks" blocks and assign the rest to
    % finrem if the blocksize is not a perfect divisor of lenquant
    cnt = cnt + 1;
    nchunks = floor(lenquant/blocksize);
    blockmean = zeros(nchunks,1);    
    
    % Compute the mean of each chunk
    for blcnt = 1:nchunks
        init = 1 + (blcnt-1)*blocksize;
        fin  = init + blocksize - 1;
        blockmean(blcnt,1) = mean(quantity(init:fin));
    end
 
    % Compute the variance in chunks (block variance) and save block
    % variance along with the number of points per chunk.
    blockvar(cnt,1)     = std(blockmean);
    sizearrblock(cnt,1) = blocksize;
    clear blockmean;
    
end

if cnt < 10
    fprintf('Warning: Small Number of Block Points\n');
end

