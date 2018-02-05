%% dz and Wind

%% Oxigen in GW aria
L=length(O);
a=480;
b=1196;

[FFT_GWO_Aria, iFFT_GWO_Aria_norm]=GravitationWave_Aria(FFT_GW_O, Trend_Oxigen, L, a, b);
% figure
%     subplot(211), plot(1:L, iFFT_GWO_Aria_norm,'r','LineWidth',2);  grid on
%         set(gca,'XLim',[0 L]);

%% Nitrogen in GW aria
a=510;
b=1250;

[FFT_GWN2_Aria, iFFT_GWN2_Aria_norm]=GravitationWave_Aria(FFT_GW_N2, Trend_Nitrogen, L, a, b);
%     hold on
%     subplot(211), plot(1:L, iFFT_GWN2_Aria_norm,'g','LineWidth',2);  grid on

%% dz
[dz, FFT_dz]=Vertical_displacement_dz(Temperature_VertWind(WATS_Vy_start:WATS_Vy_end), O, N2, L);

[Trend_dz, FFT_GW_dz, GravWave_dz]=GravitationWave(dz);
% subplot(212), plot(1:2^16, abs(FFT_GW_dz),'k','LineWidth',2); set(gca,'XLim',[0 2600]); grid on

 a=665;
 b=1767;

[FFT_GWdz_Aria, iFFT_GWdz_Aria_norm]=GravitationWave_Aria(FFT_GW_dz, Trend_dz, L, a, b);


%% Vz
Vz=VerticalWind(WATS_Vy_start:WATS_Vy_end);
	Temperature=Temperature_VertWind(WATS_Vy_start:WATS_Vy_end);
        L_Vz=UT_WATS_VertWind(WATS_Vy_end)-UT_WATS_VertWind(WATS_Vy_start); % how many seconds lost
        Length_Vz=length(Vz);
[Vz_interpolated]=Naiquist_theorem(Vz, L_Vz, Length_Vz);

     UT_WATS=(UT_WATS_VertWind(WATS_Vy_start):1:UT_WATS_VertWind(WATS_Vy_end))';
     
Vz=(Vz_interpolated(7:1647))';     

[Trend_Vz, FFT_GW_Vz, GravWave_Vz]=GravitationWave(Vz);
% subplot(211), plot(1:2^16, abs(FFT_GW_Vz),'r','LineWidth',2);

 a=575;
 b=1771;

[FFT_GWVz_Aria, iFFT_GWVz_Aria_norm]=GravitationWave_Aria(FFT_GW_Vz, Trend_Vz, L, a, b);

 %% Vy

 a_Vy=471;
 b_Vy=1425;
 L_Vy=length(Vy);
 
 NFFT=2^16;

for j=[1:a_Vy-1,(b_Vy+1):NFFT-(a_Vy+1),NFFT-(a_Vy-1):NFFT]
      FFT_GW_Aria(j)=0;
end
for j=[a_Vy:b_Vy,NFFT-b_Vy:NFFT-a_Vy]
    FFT_GWVy_Aria(j)=FFT_GW_Vy(j);
end

       
FFTOwMA1=permute(FFT_GWVy_Aria,[2 1]);
iFFT_GWVy_Aria=ifft(FFTOwMA1);
iFFT_GWVy_Aria_norm=iFFT_GWVy_Aria(1:L_Vy);

% [FFT_GWVy_Aria, iFFT_GWVy_Aria_norm]=GravitationWave_Aria(FFT_GW_Vy, Trend_Vy, L_Vy, a_Vy, b_Vy);


%% FIGUREs
figure
 subplot(211), plot(1:L,iFFT_GWO_Aria_norm,'r','LineWidth',2); grid on
 hold on
  subplot(211), plot(1:L,iFFT_GWN2_Aria_norm,'g','LineWidth',2);
    set(gca,'XLim',[0 L]);    
    
 subplot(212), plot(1:L,iFFT_GWdz_Aria_norm,'k','LineWidth',2); grid on
 hold on
  subplot(212), plot(1:L,iFFT_GWVz_Aria_norm,'b','LineWidth',2); set(gca,'XLim',[0 L]);
 hold on
  subplot(212), plot(iFFT_GWVy_Aria_norm,'r','LineWidth',2);
    

   
