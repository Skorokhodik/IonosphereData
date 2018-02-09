clear
%% create mathematic function for method proba

x=1:2024;

 figure
y1=18*10^8*exp(10^-6*(-(x-1000).^2)); % trend
        hold on    
    subplot(431); plot(x,y1,'g','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
        hold on
y2=(sin(93*pi*x./length(x))+0.5*sin(97*pi*x./length(x))+0.5*sin(109*pi*x./length(x))+0.3*sin(117*pi*x./length(x))); % Lx=260-290-320-340 km
%         hold on
%      subplot(431); plot(x,y2,'k','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
%         hold on    
y3=3*(0.5*sin(47*pi*x./length(x))+sin(51*pi*x./length(x))+0.7*sin(63*pi*x./length(x))+0.3*sin(70*pi*x./length(x))+0.4*sin(77*pi*x./length(x))); % Lx=460-520-650-750 km
%         hold on
%     subplot(431); plot(x,y3,'b','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
%         hold on
y4=7*(sin(20*pi*x./length(x))+0.5*sin(31*pi*x./length(x))); % Lx=1050,1570 km
%         hold on
%     subplot(431); plot(x,y4,'c','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
%         hold on    

%% GW without trend   
y234=8*10^6*exp(10^-5*(-(x-1000).^2)).*(y2+y3+y4);
noTrend=[y234';zeros(2^16-length(y2),1)]; % no trend
    subplot(433); plot(noTrend,'g','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
    hold on
%% spectrum by GW without trend
        sGW=length(noTrend);
        NFFTGW=2^nextpow2(sGW);                              
        FFTGW=fft(noTrend,NFFTGW);
        fGW=linspace(0,NFFTGW-1,NFFTGW);
        	subplot(432); plot(fGW,abs(FFTGW(1:NFFTGW)),'g','LineWidth',2); set(gca,'XLim',[0 2600]); grid on   
            set(gca,'YLim',[0 5*10^10]);
%% recipe of Filter
        for m=[1:(200-1),(2600+1):NFFTGW-(2600+1),NFFTGW-(200-1):NFFTGW]
            FFTGW1(m)=0;
        end
        for m=[200:2600,NFFTGW-2600:NFFTGW-200]
            FFTGW1(m)=FFTGW(m);
        end
            hold on
            subplot(432), plot(fGW,abs(FFTGW1(1:NFFTGW)),'m','LineWidth',1); grid on 
            hold on
                  
            FFTGW2=permute(FFTGW1,[2 1]);
            iFFTGW=ifft(FFTGW2);
            %norm_GW=iFFTGW./trendMA;
            hold on
            subplot(433), plot(1:NFFTGW,iFFTGW(1:NFFTGW),'r','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
            hold on
            
%% GW with trend and noises
y5=[(y1+y234+10^7*randn(size(length(x))))';zeros(2^16-length(y1),1)]; % mix with trend
    subplot(434); plot(y5,'g','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
    
%% spectrum by GW with trend
        s_realGWm=length(y5);
        NFFT_realGWm=2^nextpow2(s_realGWm);                              
        FFT_realGWm=fft(y5,NFFT_realGWm);
        f_realGWm=linspace(0,NFFT_realGWm-1,NFFT_realGWm);
            subplot(435); plot(f_realGWm,abs(FFT_realGWm(1:NFFT_realGWm)),'g','LineWidth',2); set(gca,'XLim',[0 2600]); grid on   
            set(gca,'YLim',[0 5*10^10]);
%% recipe of Filter for GW with trend
        for i=[1:(200-1),(2600+1):NFFT_realGWm-(2600+1),NFFT_realGWm-(200-1):NFFT_realGWm]
            FFTGWt(i)=0;
        end
        for i=[200:2600,NFFT_realGWm-2600:NFFT_realGWm-200]
            FFTGWt(i)=FFT_realGWm(i);
        end
                hold on
            subplot(435), plot(f_realGWm,abs(FFTGWt(1:NFFT_realGWm)),'m','LineWidth',1); grid on 
                hold on
                  
            FFTGWtrend=permute(FFTGWt,[2 1]);
            iFFT_realGWm=ifft(FFTGWtrend);
                hold on
            subplot(436), plot(1:NFFT_realGWm,iFFT_realGWm(1:NFFT_realGWm),'r','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
                hold on
            

%% our METHOD Lave=701 points
trend=[smooth(y5(1:length(x)),701,'moving');zeros(2^16-length(x),1)];
    
GW=y5-trend;
	subplot(437); plot(GW,'g','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
        s=length(GW);
        NFFT=2^nextpow2(s);                              
        FFT=fft(GW,NFFT);
        f=linspace(0,NFFT-1,NFFT);
            subplot(438); plot(f,abs(FFT(1:NFFT)),'g','LineWidth',2); set(gca,'XLim',[0 2600]); grid on   
            set(gca,'YLim',[0 5*10^10]);
%% recipe of Filter for GW=signal-trend
        for j=[1:(200-1),(2600+1):NFFT-(2600+1),NFFT-(200-1):NFFT]
            FFTGWmat(j)=0;
        end
        for j=[200:2600,NFFT-2600:NFFT-200]
            FFTGWmat(j)=FFT(j);
        end
            hold on
            subplot(438), plot(f,abs(FFTGWmat(1:NFFT)),'m','LineWidth',1); grid on 
            hold on
                  
            FFTGWmethod=permute(FFTGWmat,[2 1]);
            iFFT=ifft(FFTGWmethod);
            %norm_GW=iFFT./trend;
            hold on
            subplot(439), plot(1:NFFT,iFFT(1:NFFT),'r','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
            hold on                

%% CHECK OUT THE DIFFERENCEs
figure
        norm_noTrend=noTrend(1:length(x))/trend(1:length(x));
        norm_iFFT=iFFT(1:length(x))/trend(1:length(x));
        
    r_diff=(norm_iFFT-norm_noTrend)/norm_noTrend;
plot(1:length(x),r_diff,'r','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on

% %% polinom
% %% calculate trend by polinom (Polinom by Lagrange 1 power)
% xl=[1:round(length(x)/5):length(x) length(x)];
% yl=y1(xl)+y2(xl)+y3(xl)+y4(xl);
% 
% xx=linspace(1,length(x),length(x));
% yy=lagrange(xl,yl,xx);
% 
% hold on
% %subplot(4,3,10), plot(xx,yy,'b','LineWidth',1); set(gca,'XLim',[0 length(x)]); grid on
% hold on
% 
% GWpol=[y5(1:length(x))-yy';zeros(2^16-length(x),1)];
% subplot(4,3,10), plot(GWpol,'b','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
% hold on
% 
%                 
%                 
% %% Alla's METHODICS
% % trend71=[smooth(y5(1:length(x)),71,'moving');zeros(2^16-length(x),1)];
% %     
% % GW71=y5-trend71;
% % 	subplot(4,3,10); plot(GW71,'g','LineWidth',2); set(gca,'XLim',[0 length(x)]); grid on
%         s71=length(GWpol);
%         NFFT71=2^nextpow2(s71);                              
%         FFT71=fft(GWpol,NFFT71);
%         f71=linspace(0,NFFT71-1,NFFT71);
%             subplot(4,3,11); plot(f71,abs(FFT71(1:NFFT71)),'g','LineWidth',2); set(gca,'XLim',[0 2600],'YLim',[0 6*10^7]); grid on   
%             
% %% recipe of Filter for GW=signal-trend71
%         for j=[1:(200-1),(2600+1):NFFT71-(2600+1),NFFT71-(200-1):NFFT71]
%             FFTGW71(j)=0;
%         end
%         for j=[200:2600,NFFT_realGWm-2600:NFFT_realGWm-200]
%             FFTGW71(j)=FFT71(j);
%         end
%             hold on
%             subplot(4,3,11), plot(f71,abs(FFTGW71(1:NFFT71)),'m','LineWidth',1); grid on 
%             hold on
%                   
%             FFTGWm71=permute(FFTGW71,[2 1]);
%             iFFT71=ifft(FFTGWm71);
%             %norm_GW=iFFT./trend;
%             hold on
%             subplot(4,3,12), plot(1:NFFT71,iFFT71(1:NFFT71),'r','LineWidth',2); set(gca,'XLim',[0 length(x)],'YLim',[-4*10^5 4*10^5]); grid on
%             hold on

