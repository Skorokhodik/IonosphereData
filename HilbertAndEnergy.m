%% 1
figure
        
subplot(211), plot(1:L_Ox,iFFT_GW_O_1(1:L_Ox),'r','LineWidth',1); grid on
        hold on
              plot(1:L_Ox,abs(hilbert(iFFT_GW_O_1(1:L_Ox))),'g','LineWidth',2);
                title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);
                set(gca,'XLim',[0 L_Ox]);
subplot(212), plot(1:L_Ox,abs(hilbert(iFFT_GW_O_1(1:L_Ox))).^2,'k','LineWidth',2); grid on
                set(gca,'XLim',[0 L_Ox]);

figure

subplot(211), imagesc(abs(hilbert(iFFT_GW_O_1(1:L_Ox)))); colorbar;
	title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

subplot(212), imagesc(abs(hilbert(iFFT_GW_O_1(1:L_Ox))).^2); colorbar;


%% 2
figure
        
subplot(211), plot(1:L_Ox,iFFT_GW_O_2(1:L_Ox),'r','LineWidth',1); grid on
        hold on
              plot(1:L_Ox,abs(hilbert(iFFT_GW_O_2(1:L_Ox))),'g','LineWidth',2);
                title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);
                set(gca,'XLim',[0 L_Ox]);
subplot(212), plot(1:L_Ox,abs(hilbert(iFFT_GW_O_2(1:L_Ox))).^2,'k','LineWidth',2); grid on
                set(gca,'XLim',[0 L_Ox]);

figure

subplot(211), imagesc(abs(hilbert(iFFT_GW_O_2(1:L_Ox)))); colorbar;
	title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

subplot(212), imagesc(abs(hilbert(iFFT_GW_O_2(1:L_Ox))).^2); colorbar;


%% 3
figure
        
subplot(211), plot(1:L_Ox,iFFT_GW_O_3(1:L_Ox),'r','LineWidth',1); grid on
        hold on
              plot(1:L_Ox,abs(hilbert(iFFT_GW_O_3(1:L_Ox))),'g','LineWidth',2);
                title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);
                set(gca,'XLim',[0 L_Ox]);
subplot(212), plot(1:L_Ox,abs(hilbert(iFFT_GW_O_3(1:L_Ox))).^2,'k','LineWidth',2); grid on
                set(gca,'XLim',[0 L_Ox]);

figure

subplot(211), imagesc(abs(hilbert(iFFT_GW_O_3(1:L_Ox)))); colorbar;
	title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

subplot(212), imagesc(abs(hilbert(iFFT_GW_O_3(1:L_Ox))).^2); colorbar;

