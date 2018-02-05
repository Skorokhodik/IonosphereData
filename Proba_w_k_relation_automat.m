%% Looking for minimum in spectra and nearby min in all parameters
[noZero_min_O]=MinimumInSpectra(NFFT, FFT_GW_O);
    [noZero_min_N2]=MinimumInSpectra(NFFT, FFT_GW_N2);
        [noZero_min_dz]=MinimumInSpectra(NFFT, FFT_dz);
            [noZero_min_Vz]=MinimumInSpectra(NFFT, FFT_GW_Vz);
                [noZero_min_Vy]=MinimumInSpectra(NFFT, FFT_GW_Vy);
                    [noZero_min_dp]=MinimumInSpectra(NFFT, FFT_GW_dp);
                        [noZero_min_dT]=MinimumInSpectra(NFFT, FFT_GW_T_Vz);
%                                 SpectaMin=importdata('E:\Sciense\DISER\work in Matlab\programs\METHODICS_2016\1982333T152820\1982333T152820_SpecMin.txt');
                    
                   for i=1:length(noZero_min_N2)
                    minimum_spectra_all(i,:)=[noZero_min_O(i) noZero_min_N2(i) noZero_min_dz(i) noZero_min_Vz(i) noZero_min_Vy(i) noZero_min_dp(i) noZero_min_dT(i)];
                   end
    
    for i=1:length(noZero_min_O)
        for j=1:length(noZero_min_N2)
if (noZero_min_O(i)-noZero_min_N2(j))<50
	cross_min_O_N2(j,:)=[noZero_min_O(i) noZero_min_N2(j)];
end
        end
    end
               
    for i=1:length(noZero_min_O)
        for j=1:length(noZero_min_dz)
if (noZero_min_O(i,1)-noZero_min_dz(j))<40
	cross_min_O_dz(j,:)=[noZero_min_O(i) noZero_min_dz(j)];
end
        end
    end
        
    for i=1:length(noZero_min_O)
        for j=1:length(noZero_min_Vz)
if (noZero_min_O(i,1)-noZero_min_Vz(j))<40
	cross_min_O_Vz(j,:)=[noZero_min_O(i) noZero_min_Vz(j)];
end
        end
    end

    for i=1:length(noZero_min_O)
        for j=1:length(noZero_min_Vy)
if (noZero_min_O(i,1)-noZero_min_Vy(j))<40
	cross_min_O_Vy(j,:)=[noZero_min_O(i) noZero_min_Vy(j)];
end
        end
    end
    
    for i=1:length(noZero_min_O)
        for j=1:length(noZero_min_dp)
if (noZero_min_O(i,1)-noZero_min_dp(j))<40
	cross_min_O_dp(j,:)=[noZero_min_O(i) noZero_min_dp(j)];
end
        end
    end
    
    for i=1:length(cross_min_O_N2)
        for j=1:length(cross_min_O_dz)   
            for k=1:length(cross_min_O_Vz)  
if (cross_min_O_N2(i,1)==cross_min_O_dz(j,1))&&(cross_min_O_N2(i,1)==cross_min_O_Vz(k,1))
    cross_min_O_N2_dz_Vz(i,:)=[cross_min_O_N2(i,1) cross_min_O_N2(i,2) cross_min_O_dz(j,2) cross_min_O_Vz(k,2)];
end
            end
        end
    end
    
    
for i=1:length(cross_min_O_N2_dz_Vz)
	for j=1:length(cross_min_O_Vy) 
        for k=1:length(cross_min_O_dp)
if (cross_min_O_N2_dz_Vz(i,1)==cross_min_O_Vy(j,1))&&(cross_min_O_N2_dz_Vz(i,1)==cross_min_O_dp(k,1))
    cross_min(i,:)=[cross_min_O_N2_dz_Vz(i,1) cross_min_O_N2_dz_Vz(i,2) cross_min_O_N2_dz_Vz(i,3) cross_min_O_N2_dz_Vz(i,4) cross_min_O_Vy(j,2) cross_min_O_dp(k,2)];
end
        end
    end
end

noZero=find(cross_min(:,1));
for i=1:length(k)
    minimum_spectra=cross_min(noZero,:);
end

for i=1:length(minimum_spectra)
minimum_spectra_and_km(i,:)=[minimum_spectra(i,:) 2^16*7.8/minimum_spectra(i,1)];
end