function [breake_points]=BreakePoints(UT, sampling)

	breake_points=[length(UT)-1,length(UT),UT(end)-UT(end-1)];
	% Looking for break points and put in massive [UT(k-1) UT(k) BreakeSeconds]
for i=1:length(UT)-1            
    diff_UT=UT(i+1)-UT(i); % normally diff_UT=1, ...
    samp_diff(i,:)=[i i+1 diff_UT]; % we are looking for unnormal events
end

for j=1:length(samp_diff)
    if samp_diff(j,3)~=sampling
       breake_points=[breake_points; samp_diff(j,:)];
    end
end
end