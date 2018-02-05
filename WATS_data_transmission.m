function [K_wats]=WATS_data_transmission(WATS)

l_WATS=length(WATS(:,1));   % determinate data length
N_wats=permute(WATS,[2 1]); % turn around massive

for i_wats=2:2:l_WATS
    n1=N_wats(:,i_wats-1);  % select each second element 
    n2=N_wats(:,i_wats);
      n_wats=[n1;n2];
        k_wats(:,i_wats/2)=[n_wats];
end
 
K_wats=k_wats';% turn back
end


                
                

