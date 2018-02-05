function [VerticalWind_data]=Vertical_wind_WATS(K_wats)

L_WATS=length(K_wats(:,3));

%% Vertical wind data selected from all WATS data set
w=K_wats(:,3)./100; % key for separated Vz and Vy
j_wats=0;
% for condition of key w = 5 and 6, than it is Vz data sets
        for k=1:L_WATS
            if w(k,:)>4
                Kvert=K_wats(k,:);
                j_wats=j_wats+1;
            end
                if j_wats~=0
            VertWind(j_wats,:)=[Kvert];
                end
        end
        
%% Vertical WIND data set
        VerticalWind_data=VertWind; % means of Vertical Wind
end
