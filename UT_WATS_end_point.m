function [WATS_end]=UT_WATS_end_point(UT_NACS_end, UT_WATS)
 
L=length(UT_WATS);
% Looking for end point where WATS & NACS data agreement or close
    for i=1:L-1
            if (UT_WATS(i)<=UT_NACS_end) && (UT_NACS_end<=UT_WATS(i+1))
                WATS_end=i;
            elseif (UT_WATS(end)==UT_NACS_end) || (UT_WATS(end)<UT_NACS_end)
                WATS_end=L;
            end
    end
end
