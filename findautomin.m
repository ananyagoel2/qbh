function shift = findautomin(s, i, start, stop, corrthresh)

shift = start;
daprev = 0;
da = 0;
sprev = 0;

tot = 1;

if(max( s(i:i+1200,1) ) > 0.1)
    
    while( (~(daprev < 0 & da > 0) | ~( sprev < corrthresh)) & shift < stop)
        sprev = tot;
        tot = 0;
        daprev = da;
        
        for j = 0 : 599
            tot = tot + abs( s(i + j, 1) - s(i +j + shift, 1));
        end
        
        tot = tot / 600;
        shift = shift + 1;
        da = tot - sprev;
        
    end
    if(shift == stop)
        shift = -1;
    else
        shift = shift - 2;
    end
    
else
    shift = -1;
end