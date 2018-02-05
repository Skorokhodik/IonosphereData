%% Wave in diff parameters and Hilbert
figure
% Wave
subplot(311), plot(UT_NACS_1sec./3600,iFFT_GW_O_1,'Color',[1,0.5,0.1],'LineWidth',1); grid on
 hold on
subplot(311), plot(UT_NACS_1sec./3600,iFFT_GW_N2_1,'Color',[0.2,1,0.5],'LineWidth',1);
        hold on
    subplot(311), plot(UT_NACS_1sec./3600,HILBERT_1,'Color',[0 0 0],'LineWidth',1);
legend('[O]','[N2]','Hilbert');
xlabel('UT, hours','FontName','Raleway Light','fontsize',10);
ylabel('Relative consentration','FontName','Raleway Light','fontsize',10)
title(['Datafile ' dayOrbit  ', Orbit ' num2str(Orbit(1)) ', UT start ' num2str(UT_NACS(1)/3600) 'hour'],...
    'fontsize',10);
set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600],...
    'YLim',[min(real(iFFT_GW_N2_1)) max(real(iFFT_GW_N2_1))]);

        relation_dz_Vz=rms(abs(iFFT_GW_dz_1))/rms(abs(iFFT_GW_Vz_1));
	subplot(312), plot(UT_NACS_1sec./3600,iFFT_GW_dz_3(1:L_Ox),'Color',[0.2,0.8,0.5],'LineWidth',1); grid on
    hold on
    subplot(312), plot(UT_WATS_Vz_1sec./3600,iFFT_GW_Vz_3(1:length(UT_WATS_Vz_1sec)).*relation_dz_Vz,...
        'Color',[0,1,0],'LineWidth',1);
    legend('dz',['V_z' '*' num2str(round(relation_dz_Vz))]);
    xlabel('UT, hours','FontName','Raleway Light','fontsize',10);
    ylabel('dz, meter','FontName','Raleway Light','fontsize',10)
    set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);

        	relation_Vy_Vz=rms(abs(iFFT_GW_Vy_1))/rms(abs(iFFT_GW_Vz_1));
            relation_Vy_dp=rms(abs(iFFT_GW_Vy_1))/rms(abs(iFFT_GW_dp_1));
        subplot(313), plot(UT_WATS_Vy_1sec(1:L_Vy)./3600,iFFT_GW_Vy_1(1:L_Vy),...
            'Color',[0.8,0.2,1],'LineWidth',1); grid on
        hold on
        subplot(313), plot(UT_WATS_Vz_1sec./3600,iFFT_GW_Vz_1(1:length(UT_WATS_Vz_1sec)).*relation_Vy_Vz,...
            'Color',[0,1,0],'LineWidth',1);
        hold on
        subplot(313), plot(UT_WATS_Vz_1sec./3600,iFFT_GW_dp_1(1:length(dp_p)).*relation_Vy_dp,...
            'Color',[0,0.7,0.9],'LineWidth',1);
        legend('Vy',['V_z' '*' num2str(round(relation_Vy_Vz))],['dp' '*' num2str(round(relation_Vy_dp))]);
        xlabel('UT, hours','FontName','Raleway Light','fontsize',10);
        ylabel('V_y, m/s','FontName','Raleway Light','fontsize',10)
        set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);



%% in space with points
figure
% Wave
plot(1:L_Ox,iFFT_GW_O_2,'Color',[1,0.5,0.1],'LineWidth',1); grid on
 hold on
plot(1:L_Ox,iFFT_GW_N2_2,'Color',[0.2,1,0.5],'LineWidth',1);
        hold on
plot(1:L_Ox,HILBERT_2,'Color',[0 0 0],'LineWidth',1);
legend('[O]','[N2]','Hilbert');
set(gca,'XLim',[0 L_Ox]);


%% noise level in space

for i=1:L_Ox
    s=rms(FFT_noise)/rms(FFT_GW_O);
    SpaceLineNoise(i)=s*rms(GravWave_Oxigen);
end
% hold on, subplot(212), plot(1:L_Ox, SpaceLineNoise,'k','LineWidth',1)


 % max&min HILBERT
                                %HILBERT=abs(hilbert(iFFT_GW_dO_O_Aria(1:L_Ox)));
	for h=2:L_Ox-1
        if ((HILBERT_2(h)-HILBERT_2(h-1))>0)&&((HILBERT_2(h+1)-HILBERT_2(h))<0)
            HILBERT_2_maximum(h)=h;
            
        end
    end
                                        
        for h=2:L_Ox-1
            if ((HILBERT_2(h)-HILBERT_2(h-1))<0)&&((HILBERT_2(h+1)-HILBERT_2(h))>0)
                HILBERT_2_minimum(h)=h;
            end
        end
HILBERT_2_max=find(HILBERT_2_maximum);
HILBERT_2_min=find(HILBERT_2_minimum);


%% Maximum Amplitude and in what Latitudional aria is there
a=[10 524 810 1506 L_Ox];
            Lat_aria=[0 0];
    for i=1:length(a)-1
maxA(i)=max(abs(iFFT_GW_O_2(a(i):a(i+1))));
Lat_aria(i,:)=[Lat(a(i)) Lat(a(i+1))];

    end

%% SPECTRUM OF each LAGW
figure
% Wave
subplot(211), plot(1:L_Ox,iFFT_GW_O_2,'r','LineWidth',1); grid on
set(gca,'XLim',[0 L_Ox]);
        hold on

        colors=['m' 'c' 'b' 'g' 'k' 'y'];
for i=1:length(a)-1
    hold on
subplot(211), plot(a(i):a(i+1),HILBERT_2(a(i):a(i+1)),'Color',colors(i),'LineWidth',2);
end

subplot(211), title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);
    
% Spectrum
    GW_scale_zero=[iFFT_GW_O_2; zeros(2^16-length(iFFT_GW_O_2),1)];
        FFT_GW_scale=fft(GW_scale_zero,NFFT);
        subplot(212), plot(0:NFFT-1, abs(FFT_GW_scale),'r','LineWidth',1); grid on
        set(gca,'XLim',[0 3000]);

        FFT_LWP_summ=0;
for i=1:length(a)-1

    LWP_zero=[iFFT_GW_O_2(a(i):a(i+1)); zeros(2^16-length(iFFT_GW_O_2(a(i):a(i+1))),1)];
        FFT_LWP=fft(LWP_zero,NFFT);

        FFT_LWP_summ=FFT_LWP_summ+abs(FFT_LWP);

	hold on
subplot(212), plot(0:NFFT-1, abs(FFT_LWP),'Color',colors(i),'LineWidth',2);
  	hold on
    
end
  	hold on
subplot(212), plot(0:NFFT-1, FFT_LWP_summ,'r','LineWidth',2);

  