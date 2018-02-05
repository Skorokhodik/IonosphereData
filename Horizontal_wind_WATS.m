function [HorizontalWind]=Horizontal_wind_WATS(K_wats)

L=length(K_wats(:,3));

%% Vertical wind data selected from all WATS data set
w=round(K_wats(:,3)./100);  % separate key
j=0;
% for condition w = 4 and more, than it is Vz data sets
        for k=1:L
            if w(k,:)<5
                Khor=K_wats(k,:);
                j=j+1;
            end
                if j~=0
            HorWind(j,:)=[Khor];
                end
        end
        
%% Vertical WIND

        HorizontalWind=HorWind;
end
