function [WATS_start]=Aria_WATS(UT_NACS_1, UT_WATS)
 
L=length(UT_WATS);
% Looking for aria where WATS & NACS data agreement and interpolated wind data from WATS

if UT_WATS(1)<UT_NACS_1 % NACS start work later then WATS
    for i=1:L-1
        if (UT_WATS(i)<=UT_NACS_1) && (UT_NACS_1<=UT_WATS(i+1))
            WATS_start=i+1;
        end
    end
    
elseif UT_WATS(1)>=UT_NACS_1 % NACS start work earler then WATS
        WATS_start=1;

end


end


