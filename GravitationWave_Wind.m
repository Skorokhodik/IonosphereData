function [Trend, wave, FFT, FFTp, GravWave]=GravitationWave_Wind(Wind)

 L=length(Wind);
 V=[Wind'; zeros(2^16-L,1)];
        %% trend by moving average
        Trend=[smooth(Wind,701,'moving');zeros(2^16-L,1)]; % calculate trend for data set
        
        wave=V-Trend; % trend removing from data set
        
%% fft by wave, what without trend
        NFFT=2^16;                              
        FFT=fft(wave,NFFT); % f=linspace(0,NFFT-1,NFFT);
 
        %% GW creation just for DE 2 NACS 365T042320

        startPoint=round(NFFT*7.8/2500); % in spectr with L=2500 km
        finishPoint=round(NFFT*7.8/100); % in spectr with L=250 km
        
        for n=[1:startPoint-1,(finishPoint+1):NFFT-(finishPoint+1),NFFT-(startPoint-1):NFFT]
            FFTsAria(n)=0;
        end
        for n=[startPoint:finishPoint,NFFT-finishPoint:NFFT-startPoint]
            FFTsAria(n)=FFT(n);
        end
        
            FFTp=permute(FFTsAria,[2 1]);
            iFFT=ifft(FFTp);
            GravWave=iFFT; % norm wave by height


end