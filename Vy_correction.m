function [Vy_corr]=Vy_correction(Vy, Lat_Vy)
                                    %Lat_Vy = Latitude_Vy
k=0;
        
for i=2:length(Lat_Vy)-1
    if (Lat_Vy(i)-Lat_Vy(i-1))*(Lat_Vy(i+1)-Lat_Vy(i))<0
        k=[i;k];
    end
end
    if length(k)>2
        Vy_corr=[Vy(1:k(2));(-1)*Vy(k(2)+1:k(1));Vy(k(1)+1:end)];
    else
        Vy_corr=[Vy(1:k(1));(-1)*Vy(k(1)+1:end)];
    end
end