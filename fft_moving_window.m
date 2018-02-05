clear
%% filtr FATKULLIN
    % 72 signal point = 561 km
    for i=1:72
        y(i)=1;
    end
        for i=length(y)+1:1642;
            y(i)=0;
        end

    subplot(2,2,1), plot(y,'b','LineWidth',2); grid on
    set(gca,'XLim',[0 length(y)]);

    % fft by 
    y_long = [y';zeros(2^16-length(y),1)];

        NFFT = 2^16;                              
    FFTO = fft(y_long,NFFT)/72;

    subplot(2,2,2), plot(1:NFFT,1-FFTO(1:NFFT),'b','LineWidth',2); grid on
        set(gca,'XLim',[0 2600]);


    % 11 signal points = 85 km
    for j=1:11
        y1(j)=1;
    end
        for j=length(y1)+1:1642
            y1(j)=0;
        end

    subplot(2,2,3), plot(y1,'r','LineWidth',2); grid on
    set(gca,'XLim',[0 length(y1)]);

    % fft by 
    y1_long = [y1';zeros(2^16-length(y1),1)];
        FFTO1 = fft(y1_long,NFFT)/11;

    subplot(2,2,4), plot(1:NFFT,(1-FFTO(1:NFFT)).*FFTO1(1:NFFT),'r','LineWidth',2); grid on
    set(gca,'XLim',[0 2600]);




%% FILTR LIZUNOV-SKOROKHOD
    % %% 72 signal point = 561 km
    % for i=1:701
    %     y2(i)=1;
    % end
    %     for i=length(y2)+1:1642
    %         y2(i)=0;
    %     end
    % 
    % subplot(2,2,1), plot(y2,'b','LineWidth',2); grid on
    % set(gca,'XLim',[0 length(y2)]);
    % 
    % %% fft by 
    % y2_long = [y2';zeros(2^16-length(y2),1)];
    % 
    %     NFFT = 2^16;                              
    % FFT2 = fft(y2_long,NFFT)/701;
    % 
    % subplot(2,2,2), plot(1:NFFT,1-FFT2(1:NFFT),'b','LineWidth',2); grid on
    %     set(gca,'XLim',[0 2600]);
