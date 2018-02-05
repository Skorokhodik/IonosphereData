function [Trend, wave, FFT, GravWave]=GravitationWave(Element)

        L=length(Element);
O = [Element; zeros(2^16-L,1)];
        %% trend by moving average
Trend=[smooth(Element,701,'moving');zeros(2^16-L,1)]; % calculate trend for data set
        
wave=O-Trend; % trend removing from data set
        
%% fft by wave, what without trend
            NFFT=2^16;                              
        FFT=fft(wave, NFFT);
 
        %% GW creation just for DE 2 NACS 365T042320

        a=round(NFFT*7.8/2500); % in spectr with L=2500 km
        b=round(NFFT*7.8/100); % in spectr with L=250 km
        
        for n=[1:a-1,(b+1):NFFT-(b+1),NFFT-(a-1):NFFT]
            FFTsAria(n)=0;
        end
        for n=[a:b,NFFT-b:NFFT-a]
            FFTsAria(n)=FFT(n);
        end
        
            FFTp=permute(FFTsAria,[2 1]);
            iFFT=ifft(FFTp);
            GravWave=iFFT(1:L)./Trend(1:L); % norm wave by height
            
end