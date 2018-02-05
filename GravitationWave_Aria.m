function [FFT_GW_Aria, iFFT_GW_Aria_norm]=GravitationWave_Aria(FFT, Trend, L, a, b)

NFFT=2^16;

for j=[1:a-1,(b+1):NFFT-(b+1),NFFT-(a-1):NFFT]
      FFT_GW_Aria(j)=0;
end
for j=[a:b,NFFT-b:NFFT-a]
    FFT_GW_Aria(j)=FFT(j);
end

       
FFTOwMA1=permute(FFT_GW_Aria,[2 1]);
iFFT_GW_Aria=ifft(FFTOwMA1);
iFFT_GW_Aria_norm=iFFT_GW_Aria(1:L)./Trend(1:L);

end