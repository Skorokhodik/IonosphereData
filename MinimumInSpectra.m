function [noZero_min]=MinimumInSpectra(NFFT, FFT)


for i=round(NFFT*7.8/7000):round(NFFT*7.8/100)
   if (abs(FFT(i))-abs(FFT(i-1))<0)&&(abs(FFT(i+1))-abs(FFT(i))>0)
      minimum_poins_spectra(i,:)=[i-1 abs(FFT(i-1))];
   end
end

noZero_min=find(minimum_poins_spectra)-1;

end