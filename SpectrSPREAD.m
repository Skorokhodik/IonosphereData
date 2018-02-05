clear

x=linspace(0,10,1000);
y=sin(2*8*pi*x/10);

    subplot(2,2,1), plot(x,y,'b','LineWidth',2); grid on
        set(gca,'XLim',[0 max(x)]);

    % fft by 
    y_long = [y';zeros(2^16-length(y),1)];

        NFFT = 2^16;                              
    FFTO = fft(y_long,NFFT);

    subplot(2,2,2), plot(1:NFFT,abs(FFTO),'b','LineWidth',2); grid on
        set(gca,'XLim',[0 30]);
        
        
%% Spectr Spread

y1=exp(-5*10^-9*(x-2^15).^2).*sin(2*8*pi*x/2^16);

    subplot(2,2,3), plot(x,y1,'b','LineWidth',2); grid on
        set(gca,'XLim',[0 max(x)]);

    % fft by 
    y1_long = [y1';zeros(2^16-length(y1),1)];

    FFTO1 = fft(y1_long,NFFT);

    subplot(2,2,4), plot(1:NFFT,abs(FFTO1),'b','LineWidth',2); grid on
        set(gca,'XLim',[0 30]);