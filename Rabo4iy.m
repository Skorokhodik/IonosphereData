%% Trend checking as middle line and polinom
clear
    way_NACS='E:\Sciense\DISER\DATA files\Data file of AGW\DE2\neutral_gas_nacs\n_T_1s_ascii\'; % way on disk data are there
    day='1983031';
    dayOrbit=[day 'T141140'];
    datafile=[dayOrbit '_0_DE2_NACS_1S_V01.ASC']; % datafile name
    dataway=[way_NACS datafile];
%% DOWNLOAD datafile    
NACS=importdata(dataway);

%% Separate parameters to different dataset
UT_NACS1=round((NACS.data(:,1))./1000);	% from msec to sec
  UT_NACS=UT_NACS1(448:end);
    Alt=NACS.data(:,13);                % km
    Lat=NACS.data(:,14);                % deg 
    Long=NACS.data(:,15);               % deg
    LST=NACS.data(:,16);                % hr
Oxigen1=NACS.data(:,2);                  % cm-3
  Oxigen=Oxigen1(448:end);
    Oxigen_error=NACS.data(:,3);
Nitrogen1=NACS.data(:,4);                % cm-3
  Nitrogen=Nitrogen1(448:end);
Helium=NACS.data(:,6);                  % cm-3
Argentum=NACS.data(:,8);                % cm-3
    Orbit=NACS.data(:,12);

%% looking for breaking point (BP) in NACS data
    sampling_NACS=1;
    
[breake_points_NACS]=BreakePoints(UT_NACS, sampling_NACS);

%% NACS_O
if length(breake_points_NACS(:,1))>1
    O=Oxigen(1:breake_points_NACS(2,1)); % why is it impotant?.. if I need not full data set, 
    N2=Nitrogen(1:breake_points_NACS(2,1));
    UT_NACS_1sec=UT_NACS(1:breake_points_NACS(2,1)); % all UT data in NACS with 1 sec from start to breake point
else
	O=Oxigen(1:breake_points_NACS(1,2));  
    N2=Nitrogen(1:breake_points_NACS(1,2));
    UT_NACS_1sec=UT_NACS(1:breake_points_NACS(1,2));
end
L_Ox=length(O); % length of data set after BP ... so length of data have changed

        % function create trend and detect GW
        [Trend_O, Wave_O, FFT_GW_O, GravWave_Oxigen]=GravitationWave(O);
        
%% NACS_N2
        % function create trend and detect GW
        [Trend_N2, Wave_N2, FFT_GW_N2, GravWave_Nitrogen]=GravitationWave(N2);
        

figure % data from SC, trend and GW in Oxigen and Nitrogen
    
        % Concentration and trend for O and N2        
    subplot(221), plot(UT_NACS_1sec./3600,O,'r','LineWidth',2); grid on
        hold on
    subplot(221), plot(UT_NACS_1sec./3600,Trend_O(1:L_Ox),'m','LineWidth',2);
        hold on
    subplot(221), plot(UT_NACS_1sec./3600,N2,'g','LineWidth',2); grid on
        hold on
    subplot(221), plot(UT_NACS_1sec./3600,Trend_N2(1:L_Ox),'c','LineWidth',2);        
            set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
            xlabel('Oxigen variations - red line, Oxigen trend - magenta line, Nitrogen - green line, Nitrogen trend - cyant line','fontsize',12);
            title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

        % O and N2 spectrum
    subplot(222), plot(0:2^16-1, abs(FFT_GW_O),'r','LineWidth',2); grid on
        hold on
    subplot(222), plot(0:2^16-1, abs(FFT_GW_N2),'g','LineWidth',2);
            set(gca,'XLim',[0 2600]);
            xlabel('Oxigen (red line), N2 (green line)','fontsize',12);  
            
        % iFFT wave in GW area in O and N2 data
    subplot(212), plot(UT_NACS_1sec./3600, GravWave_Oxigen(1:L_Ox),'r','LineWidth',2); grid on; 
        hold on
    subplot(212), plot(UT_NACS_1sec./3600, GravWave_Nitrogen(1:L_Ox),'g','LineWidth',2);    
            set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
            xlabel('Oxigen (red line) and Nitrogen (green line) gravitation wave area (from spectrum)','fontsize',12);
   
                
%% WATS
        way_WATS='E:\Sciense\DISER\DATA files\Data file of AGW\DE2\neutral_gas_wats\';
        WATS_file_type='_de2_wats_2s_v01.asc';

WATS=textread([way_WATS day WATS_file_type]); % whole way to the datafile

[K_wats]=WATS_data_transmission(WATS);

        % Date_WATS=K_wats(:,1); 
                    % [yyddd], 82365
        % UT_WATS=K_wats(:,2); 
                    % [ms]
        % Mode_Slot_Outin_WATS=K_wats(:,3); 
                    % tree parameters Mode = 3,4 (horisontal wind) 5,6 (vertical)
                    % 1,2,3,4 steps during each 8 sec
                    % 1 - buffle going out, 0 - in
        % Mass_WATS=K_wats(:,4); 
                    % 28 or 32
        % Density_WATS=K_wats(:,5);
                    % [cm-3] Density neutrals with Mass upper
        % Tn_WATS=K_wats(:,6);
                    % [K] Neutral temperature
        % Tn_corr_WATS=K_wats(:,7);
                    % [K] after WATSCOR correction
        % Wind_WATS=K_wats(:,8);
                    % [m/s] in spacecraft spacecraft coordinates
        % =K_wats(:,9); C1
        % =K_wats(:,10); C2
        % =K_wats(:,11); T1 time & T2 time
        % Wind_geo_WATS=K_wats(:,12); 
                    % [m/s]  in corotating Earth frame (positive in the eastward/upward direction)
        % Wind_geo_corr_WATS=K_wats(:,13);
                    % after WATSCOR correction*
        % Orbit_WATS=K_wats(:,14);
        % Altitude_WATS=K_wats(:,15);
        % Latitude_WATS=K_wats(:,16);
        % Longitude_WATS=K_wats(:,17);
        % LST_WATS=K_wats(:,18);
        % ...


%% Vertical wind (Vz)
[WATS_Vz]=Vertical_wind_WATS(K_wats);

    VerticalWind=WATS_Vz(:,12); % [m/s]
    UT_WATS_VertWind=WATS_Vz(:,2); % [ms]
        UT_WATS_Vz=round(UT_WATS_VertWind./1000); % [sec]
    Temperature_Vz=WATS_Vz(:,6); 

% Harmonization of data NACS and WATS by UT
[Vz_start]=Aria_WATS(UT_NACS_1sec(1), UT_WATS_Vz);
[Vz_end]=UT_WATS_end_point(UT_NACS_1sec(end), UT_WATS_Vz);

%Looking for breake poins there
sampling_WATS_Vz=8;
[breake_points_Vz]=BreakePoints(UT_WATS_Vz(Vz_start:Vz_end), sampling_WATS_Vz);

if (length(breake_points_Vz(:,1))>1) && (breake_points_Vz(2,1)~=breake_points_Vz(1,1))
    	Vz_end=Vz_start+breake_points_Vz(2,1)-1;
elseif (length(breake_points_Vz(:,1))==1) || (breake_points_Vz(2,1)==breake_points_Vz(1,1))
    	Vz_end=Vz_end;
end

    Vz=VerticalWind(Vz_start:Vz_end);
	T_Vz=Temperature_Vz(Vz_start:Vz_end);
        L_Vz=UT_WATS_Vz(Vz_end)-UT_WATS_Vz(Vz_start); % how many seconds lost ... so how long should be data set
        Length_Vz=length(Vz); % carent data set with sampling = 8 second
[Vz_interpolated]=Naiquist_theorem(Vz, L_Vz, Length_Vz);
[T_Vz_interpolated]=Naiquist_theorem(T_Vz, L_Vz, Length_Vz);

     UT_WATS_Vz_1sec=(UT_WATS_Vz(Vz_start):1:UT_WATS_Vz(Vz_end)-1)'; % each one second from the startVz to endVz
     
[Trend_Vz, wave_Vz, FFT_Vz, FFT_GW_Vz, GravWave_Vz]=GravitationWave_Wind(Vz_interpolated);
    L_Vz=length(Vz_interpolated);
    
[Trend_T_Vz, wave_T_Vz, FFT_dT_Vz, FFT_GW_T_Vz, dT_T]=GravitationWave_Wind(T_Vz_interpolated);

%% Horizontal wind (Vy)
[WATS_Vy]=Horizontal_wind_WATS(K_wats);

    HorizontalWind=WATS_Vy(:,12); % [m/s]
        UT_WATS_HorizWind=WATS_Vy(:,2); % [ms]
    UT_WATS_Vy=round(UT_WATS_HorizWind./1000); % [sec]
    Temperature_Vy=WATS_Vy(:,6);
    Latitude_WATS_Vy=WATS_Vy(:,16);

[Vy_start]=Aria_WATS(UT_NACS_1sec(1), UT_WATS_Vy);
    	[Vy_end]=UT_WATS_end_point(UT_NACS_1sec(end), UT_WATS_Vy);

        UT_WATS_Vy_bp=UT_WATS_Vy(Vy_start:2:Vy_end);
        sampling_WATS_Vy=8;
        [breake_points_Vy]=BreakePoints(UT_WATS_Vy_bp, sampling_WATS_Vy);

if length(breake_points_Vy(:,1))>1
    	Vy_end=Vy_start+breake_points_Vy(2,1)*2-1;
else
    	Vy_end=Vy_end;
end
% Axis of Vy changen when satellite cross poles 
%...so we nead take into account this fact
Vy=HorizontalWind(Vy_start:2:Vy_end);
Lat_Vy=Latitude_WATS_Vy(Vy_start:2:Vy_end);

[Vy_corr]=Vy_correction(Vy, Lat_Vy);
UT_WATS_Vy_8sec=UT_WATS_Vy(Vy_start:2:Vy_end);
   
        T_Vy=Temperature_Vy(Vy_start:Vy_end);
      L_Vy=UT_WATS_Vy(Vy_end)-UT_WATS_Vy(Vy_start); % how many seconds lost
      Length_Vy=length(Vy_corr);
        
[Vy_interpolated]=Naiquist_theorem(Vy_corr, L_Vy, Length_Vy);
     UT_WATS_Vy_1sec=(UT_WATS_Vy(Vy_start):1:UT_WATS_Vy(Vy_end)-1)';
     
[Trend_Vy, wave_Vy, FFT_Vy, FFT_GW_Vy, GravWave_Vy]=GravitationWave_Wind(Vy_interpolated);
     
figure % Vy and Vy_correct coz Vy axe change direction
	subplot(211), plot(UT_WATS_Vy_8sec./3600,Vy,'r','LineWidth',2); grid on
        set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
    subplot(212), plot(UT_WATS_Vy_8sec./3600,Vy_corr,'r','LineWidth',2); grid on
        set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
        xlabel('Vy and Vy_correction coz of Y axe change direction','fontsize',12);
        title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

figure % Wind and Temperature
    % Vz
        subplot(321), plot(UT_WATS_Vz_1sec./3600,Vz_interpolated,'b','LineWidth',1); grid on
        	hold on
        subplot(321), plot(UT_WATS_Vz_1sec./3600,Trend_Vz(1:L_Vz),'c','LineWidth',1);
        	hold on
        subplot(321), plot(UT_WATS_Vz_1sec./3600,wave_Vz(1:L_Vz),'b','LineWidth',2);
                set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
                xlabel('Vz (blue line)','fontsize',12);
                title(['Datafile' '   ' dayOrbit  '   ' 'Noises in data from NACS'],'fontsize',14);
        % FFT Vz-Trend
        subplot(322), plot(1:2^16,abs(FFT_Vz),'b','LineWidth',2); grid on
                set(gca,'XLim',[0 2600]);
    % Vy            
        subplot(323), plot(UT_WATS_Vy_1sec./3600,Vy_interpolated,'r','LineWidth',1); grid on
         	hold on
        subplot(323), plot(UT_WATS_Vy_1sec./3600,Trend_Vy(1:L_Vy),'m','LineWidth',1);
        	hold on
        subplot(323), plot(UT_WATS_Vy_1sec./3600,wave_Vy(1:L_Vy),'r','LineWidth',2); 
                set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
                xlabel('Vy (red line)','fontsize',12);

        % FFT Vy-Trend
        subplot(324), plot(1:2^16,abs(FFT_Vy),'r','LineWidth',2); grid on
                set(gca,'XLim',[0 2600]);
    % T
        subplot(325), plot(UT_WATS_Vz_1sec./3600,T_Vz_interpolated,'m','LineWidth',2); grid on
        	hold on
        subplot(325), plot(UT_WATS_Vz_1sec./3600,Trend_T_Vz(1:length(T_Vz_interpolated)),'m','LineWidth',2);
        	hold on
        subplot(325), plot(UT_WATS_Vz_1sec./3600,wave_T_Vz(1:length(T_Vz_interpolated)),'r','LineWidth',2); 
                set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
                xlabel('Temperature','fontsize',12);
        % FFT T
        subplot(326), plot(1:2^16,abs(FFT_dT_Vz),'m','LineWidth',2); grid on
                set(gca,'XLim',[0 2600]);
                xlabel('SPECTRUM of parameters','fontsize',12);
                title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);
              

%% dz
[dz, FFT_dz]=Vertical_displacement_dz(Temperature_Vz(Vz_start:Vz_end), GravWave_Oxigen, GravWave_Nitrogen, L_Ox);

%[Trend_dz, Wave_dz, FFT_wave_dz, FFT_GW_dz, GravWave_dz]=GravitationWave_Wind(dz');
NFFT=2^16;
    
figure %4  dp/p and dz
    % dz_normalised
	subplot(221), plot(UT_NACS_1sec./3600, dz,'k','LineWidth',2); grid on
        hold on
            set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
            xlabel('dz','fontsize',12); 
	subplot(222), plot(1:NFFT, abs(FFT_dz),'k','LineWidth',2); grid on
            set(gca,'XLim',[0 2600]);
            xlabel('Spectrum dz','fontsize',12); 

            %% dp/p = (dO+dN2)/(trendO+trendN2)+dT/T
    % p=nkT: n - from NACS (each 1 sec); T - from WATS (each 8 sec)
    % for coherence of dO/O, dN2/N2 and dT/T fluctuations
    for i_NACS=1:length(UT_NACS_1sec)
        if UT_WATS_Vz_1sec(1)==UT_NACS_1sec(i_NACS)
            s=i_NACS; % startNACS
        end
        if UT_WATS_Vz_1sec(end)==UT_NACS_1sec(i_NACS)
            e=i_NACS; % endNACS
        end
    end
    
[dp_p, FFT_GW_dp]=PressureVariation(Wave_O(s:e), Wave_N2(s:e), Trend_O(s:e), Trend_N2(s:e), dT_T(1:length(UT_WATS_Vz_1sec)));
    % dp_p    
	subplot(223), plot(UT_NACS_1sec(s:e)./3600, dp_p,'m','LineWidth',2); grid on
            xlabel('dp','fontsize',12); 
            set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
    subplot(224), plot(1:NFFT,abs(FFT_GW_dp),'m','LineWidth',2); grid on
            set(gca,'XLim',[0 2600]);
            xlabel('Spectrum dp','fontsize',12); 
            title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

            
figure % Spectrum of all parameters
    plot(0:2^16-1, abs(FFT_GW_O).*1e-10,'r','LineWidth',2); grid on
        hold on
    plot(0:2^16-1, abs(FFT_GW_N2).*1e-10,'g','LineWidth',2);
        hold on
    plot(1:2^16,abs(FFT_Vz).*1e-3,'b','LineWidth',2); 
        hold on
    plot(1:2^16,abs(FFT_Vy).*1e-4,'m','LineWidth',2);
        hold on
    plot(1:2^16,abs(FFT_GW_T_Vz).*1e-5,'r','LineWidth',1);
        hold on
    plot(1:2^16,abs(FFT_dz).*1e-6,'k','LineWidth',2);
        hold on
    plot(1:2^16,abs(FFT_GW_dp).*1e-4,'c','LineWidth',2);
        set(gca,'XLim',[0 2600],'YLim',[0 4]);
        xlabel('O (red), N2 (green), Vz (blue), Vy (magenta), T (red thik), dz (black), dp (cyant)','fontsize',12);  
        title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

        
%% NOISE calculate in fraquencies aria

                Poit100km=round(2^16*7.8/100);% 100 km - bottom line of GW
            mN=mean(abs(FFT_GW_O(Poit100km:2^16-Poit100km))); 
for i=1:2^16
	Noise(i)=mN;
end
        for i_noise=[1:Poit100km-1, 2^16-Poit100km+1:2^16]
            FFT_noise(i_noise)=mN;
        end
        for i_noise=Poit100km:2^16-Poit100km
            FFT_noise(i_noise)=FFT_GW_O(i_noise);
        end
        
        ifft_noise=ifft(FFT_noise);

        
figure % Noise in spectr - all scales less 100 km
        subplot(311), plot(0:2^16-1,abs(FFT_GW_O(1:2^16)),'r','LineWidth',1); grid on
        	hold on
        subplot(311), plot(0:2^16-1, Noise,'m','LineWidth',2);
                set(gca,'XLim',[0 2600], 'YLim', [0 max(abs(FFT_GW_O(1:2^16)))/10]);
                xlabel('Oxigen (red line), Noise (magenta line)','fontsize',12);
                title(['Datafile' '   ' dayOrbit  '   ' 'Noises in data from NACS'],'fontsize',14);
        
        subplot(312), plot(UT_NACS_1sec./3600,real(ifft_noise(1:L_Ox)),'m','LineWidth',2); grid on
                set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
                xlabel('Noise from Oxigen data','fontsize',12);
                title(['Datafile   ' dayOrbit  '   Orbit Nomber  ' num2str(Orbit(1)) '   UT start   ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);
                
% Histogramm of NOISE
                x=linspace(min(real(ifft_noise(1:L_Ox))), max(real(ifft_noise(1:L_Ox))),100);
                
        subplot(313), hist(real(ifft_noise(1:L_Ox)), x); grid on
                 xlabel('Noise from Oxigen data','fontsize',12);
