%%  H=kT/mg     wg=Cs*Kg;
	%constants
        R=8.32*1000; % gramm*m^2/mol*K*s^2
        m_O=16; %gramm/mol
        m_N2=28;%gramm/mol
        g=9.8; % m/s2;
H_O=R*mean(T_Vz_interpolated)/(m_O*g); %[m]
H_N2=R*mean(T_Vz_interpolated)/(m_N2*g); %[m]
        gamma=5/3; % gamma
wg=((1-1/gamma)*g/H_O)^(1/2); % 1/s
	Tg=2*pi/(wg*60); % min
Kg=1000/(2*H_O); % [1/m]
	Lg=2*pi/Kg; % [km]
            Cs2=gamma*R*mean(T_Vz_interpolated)/m_O; % [m^2/s^2]
Cs=sqrt(Cs2); % [m/s]

%                         SpM_O=SpectaMin(:,1)';
%                         SpM_N2=SpectaMin(:,2)';
%                         SpM_dz=SpectaMin(:,3)';
%                         SpM_Vz=SpectaMin(:,4)';
%                         SpM_Vy=SpectaMin(:,5)';
%                         SpM_dp=SpectaMin(:,6)';
%                         SpM_T=SpectaMin(:,7)';

%% O  N2  dz,[m]  Vz,[m/s]  Vy,[m/s]  dp/p
% calc=[1 2 5 13];
% for calcul=1:4
%     a_O(calcul)=SpM_O(calc(calcul));   
%                         a_O_km(calcul)=2^16*7.8/a_O(calcul);
%     a_N2(calcul)=SpM_N2(calc(calcul));    a_dz(calcul)=SpM_dz(calc(calcul));  
%     a_Vz(calcul)=SpM_Vz(calc(calcul));   a_Vy(calcul)=SpM_Vy(calc(calcul));    
%     a_dp(calcul)=SpM_dp(calc(calcul));    a_T(calcul)=SpM_T(calc(calcul)); 
% end

a_O=[223 489 767 1409];   
for i=1:4,     a_O_km(i)=2^16*7.8/a_O(i); end
a_N2=[204 491 770 1418];   a_dz=[202 490 775 1497];  
a_Vz=[204 448 690 1412];   a_Vy=[226 490 721 1347];   
a_dp=[280 462 761 1515];   a_T=[280 463 761 1518];

figure % Waves of different arias and thier spectrum

for i=1:3   
    % O & N2
[FFT_GWO_Aria, iFFT_GW_dO_O_Aria]=GravitationWave_Aria(FFT_GW_O, Trend_O, L_Ox, a_O(i), a_O(i+1));
                    % create separate aria and save data 
                    nameAria_O=['iFFT_GW_O_' num2str(i)];
                    eval([nameAria_O '= iFFT_GW_dO_O_Aria']); clear var_O;
                    var_O=eval(['iFFT_GW_O_' num2str(i)]);

                    hold on
            subplot(2,3,i), plot(UT_NACS_1sec./3600,var_O(1:L_Ox),'r','LineWidth',1); grid on
                
                middle_O(i)=(a_O(i)+a_O(i+1))/2;
                Lx(i)=2^16*7.8/middle_O(i); % km
                    l2=a_O_km(i+1);
                    l1=a_O_km(i);
                        ERROR_Kx(i)=((a_O(i+1)-a_O(i))/(2*middle_O(i)))*100; % error of Kx in [%]

                    xlabel(['L_x=' num2str(round(Lx(i))) ' km, O (red line), N_2 (green)'],'fontsize',12);
                    set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
                    
                            % Noises level   
                            relationOandNoise_FFT=rms(FFT_noise)/rms(FFT_GWO_Aria);
                                        clear Noises_Level;
                                    for j=1:length(UT_NACS_1sec)
                                        Noises_Level(j)=relationOandNoise_FFT*rms(iFFT_GW_dO_O_Aria);
                                    end
                                                     nameAria_N=['Noises_Level_' num2str(i)];
                                                     eval([nameAria_N '= Noises_Level']); clear var_N;
                                                     var_N=eval(['Noises_Level_' num2str(i)]);
                                                       
                                    
                                 HILBERT=abs(hilbert(iFFT_GW_dO_O_Aria(1:L_Ox)));
                                                     nameAria_H=['HILBERT_' num2str(i)];
                                                     eval([nameAria_H '= HILBERT']); clear var_H;
                                                     var_H=eval(['HILBERT_' num2str(i)]);

  
                                    
[FFT_GWN2_Aria, iFFT_GW_dN2_N2_Aria]=GravitationWave_Aria(FFT_GW_N2, Trend_N2, L_Ox, a_N2(i), a_N2(i+1));
                hold on
            subplot(2,3,i), plot(UT_NACS_1sec./3600,iFFT_GW_dN2_N2_Aria(1:L_Ox),'g','LineWidth',1);
                hold on, plot(UT_NACS_1sec./3600,HILBERT,'m','LineWidth',1);
%                 hold on, plot(UT_NACS_1sec./3600,var_N,'k','LineWidth',1);
%                 hold on, plot(UT_NACS_1sec./3600,-var_N,'k','LineWidth',1);

                    % create separate aria and save data 
                    nameAria_N2=['iFFT_GW_N2_' num2str(i)];
                    eval([nameAria_N2 '= iFFT_GW_dN2_N2_Aria']); clear var_N2;
                    var_N2=eval(['iFFT_GW_N2_' num2str(i)]);

                
                relation_N2_O(i)=rms(abs(iFFT_GW_dN2_N2_Aria))/rms(abs(var_O)); % relative amplitude of N2 and O
    
    % dz & Vz
[FFT_dz_GW_Aria, iFFT_dz_Aria]=GravitationWave_Aria_Wind(FFT_dz, a_dz(i), a_dz(i+1)); % [m]
%                 hold on
%             subplot(3,3,3+i), plot(UT_NACS_1sec./3600,iFFT_dz_Aria(1:L_Ox),'k','LineWidth',1); grid on
% 
                    % create separate aria and save data 
                    nameAria_dz=['iFFT_GW_dz_' num2str(i)];
                    eval([nameAria_dz '= iFFT_dz_Aria']); clear var_dz;
                    var_dz=eval(['iFFT_GW_dz_' num2str(i)]);

            
[FFT_Vz_GW_Aria, iFFT_Vz_GW_Aria]=GravitationWave_Aria_Wind(FFT_Vz, a_Vz(i), a_Vz(i+1)); % [m/s]
%                 hold on
%                         relation_dz_Vz=rms(abs(iFFT_dz_Aria))/rms(abs(iFFT_Vz_GW_Aria));
%             subplot(3,3,3+i), plot(UT_WATS_Vz_1sec./3600,iFFT_Vz_GW_Aria(1:length(UT_WATS_Vz_1sec)).*relation_dz_Vz,'b','LineWidth',1);
%                     xlabel('V_z (blue), dz (black)','fontsize',10); 
%                     set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
% 
                    % create separate aria and save data 
                    nameAria_Vz=['iFFT_GW_Vz_' num2str(i)];
                    eval([nameAria_Vz '= iFFT_Vz_GW_Aria']); clear var_Vz;
                    var_Vz=eval(['iFFT_GW_Vz_' num2str(i)]);

                    
% Vy & dp/p
[FFT_Vy_GW_Aria, iFFT_Vy_GW_Aria]=GravitationWave_Aria_Wind(FFT_Vy, a_Vy(i), a_Vy(i+1)); % [m/s]
                hold on
            subplot(2,3,3+i), plot(UT_WATS_Vy_1sec(1:L_Vy)./3600,iFFT_Vy_GW_Aria(1:L_Vy),'m','LineWidth',1); grid on
                hold on
                        relation_Vy_Vz=rms(abs(iFFT_Vy_GW_Aria))/rms(abs(iFFT_Vz_GW_Aria));
            subplot(2,3,3+i), plot(UT_WATS_Vz_1sec./3600,iFFT_Vz_GW_Aria(1:length(UT_WATS_Vz_1sec)).*relation_Vy_Vz,'b','LineWidth',1);

                    % create separate aria and save data 
                    nameAria_Vy=['iFFT_GW_Vy_' num2str(i)];
                    eval([nameAria_Vy '= iFFT_Vy_GW_Aria']); clear var_Vy;
                    var_Vy=eval(['iFFT_GW_Vy_' num2str(i)]);

            
[FFT_dp_GW_Aria,iFFT_dp_GW_Aria]=GravitationWave_Aria_Wind(FFT_GW_dp, a_dp(i), a_dp(i+1));
                hold on
                        relation_Vy_dp=rms(abs(iFFT_Vy_GW_Aria))/rms(abs(iFFT_dp_GW_Aria));
            subplot(2,3,3+i), plot(UT_WATS_Vz_1sec./3600,iFFT_dp_GW_Aria(1:length(dp_p)).*relation_Vy_dp,'c','LineWidth',1);
                    xlabel('V_y (magenta), dp_p (cyant)','fontsize',12); 
                    set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);

                    % create separate aria and save data 
                    nameAria_dp=['iFFT_GW_dp_' num2str(i)];
                    eval([nameAria_dp '= iFFT_dp_GW_Aria']); clear var_dp;
                    var_dp=eval(['iFFT_GW_dp_' num2str(i)]);

                    
[FFT_dT_GW_Aria,iFFT_GW_dT_T_Aria]=GravitationWave_Aria(FFT_GW_T_Vz, Trend_T_Vz, L_Vz, a_T(i), a_T(i+1));
                hold on

                    % create separate aria and save data 
                    nameAria_T=['iFFT_GW_T_' num2str(i)];
                    eval([nameAria_T '= iFFT_GW_dT_T_Aria']); clear var_T;
                    var_T=eval(['iFFT_GW_T_' num2str(i)]);

                
%% w=dVz(k)/dz(k)
    w(i)=rms(abs(FFT_Vz_GW_Aria))/rms(abs(FFT_dz_GW_Aria)); % 1/s
                T_minut(i)=2*pi/(w(i)*60); % min
        W0(i)=w(i)/wg;

%% Kx
            
            % Kx=2*pi/Lx;
    Kx(i)=2*pi/Lx(i); % 1/km
    	qx(i)=Kx(i)/Kg;
                        
                        dKx(i)=(1/l2-1/l1)*pi;
                        
                        dqx(i)=dKx(i)/Kg;


%% Ky=(y*w*Vy)/(Cs2*dp)
                FFT_dO_O=fft([var_O; zeros(NFFT-L_Ox,1)], NFFT);
                FFT_dN2_N2=fft([iFFT_GW_dN2_N2_Aria; zeros(NFFT-L_Ox,1)], NFFT);
                FFT_dT_T=fft([iFFT_GW_dT_T_Aria; zeros(NFFT-L_Vz,1)], NFFT);
    Ky(i)=(gamma/Cs2)*w(i)*1000*(rms(abs(iFFT_Vy_GW_Aria(1:L_Vy)))/(rms(abs(var_O))+rms(abs(iFFT_GW_dT_T_Aria)))); % 1/km
        qy(i)=Ky(i)/Kg;
        
            % check out
            check_K(i)=Ky(i)/Kx(i);
            check_V(i)=rms(abs(iFFT_Vy_GW_Aria))/rms(abs(iFFT_Vz_GW_Aria));

%% Ly=2*pi/Ky
    Ly(i)=2*pi/Ky(i); % km

%% Kh2=Kx^2+Ky^2
        Kh2=Kx(i)^2+Ky(i)^2; % horisontal wave nomber, 1/km2
	Kh(i)=sqrt(Kh2);
        qh(i)=Kh(i)/Kg;

%% Lh=2*pi/Kh
    Lh(i)=2*pi/Kh(i);

%% Kz2=(wg^2/w^2-1)*Kh2-Kg2
        Kz2(i)=(wg^2/w(i)^2-1)*Kh2-Kg^2;
    Kz(i)=sqrt(Kz2(i));
        qz(i)=Kz(i)/Kg;

%% Lz=2*pi/Kz
    Lz(i)=2*pi/Kz(i);
                            
%% K=sqrt(Kh2+Kz2)
    K(i)=sqrt(Kh2+Kz2(i));
        q(i)=K(i)/Kg;
end        
             title(['Datafile ' dayOrbit  ', Orbit ' num2str(Orbit(1)) ', UT start ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',10);
             
    filename=[dayOrbit,'_GW_WP'];
saveas(gcf, fullfile(fpath,filename),'jpeg');
    saveas(gcf, fullfile(fpath,filename),'pdf');

        
figure % DISPERTION ARIA
% create mathematic function
        q_math=0.01:0.001:10;
        W=q_math./sqrt(1+q_math.^2); % upper border condition
            plot(q_math, W,'--k','LineWidth',2); grid on
                hold on
        Wdiss=0.1*q_math.^(2/3); % dissipative bottom border condition
                hold on
            plot(q_math, Wdiss,'--k','LineWidth',2);
                hold on
                
        ws=sqrt(1+q_math.^2); %normalaised accoustic
                hold on
            plot(q_math, ws,'--k','LineWidth',2);
            
        cs=q_math; % line of sound velosity
                hold on
            plot(q_math, cs,'--k','LineWidth',2);

        cg=q_math.*0.98; % line of sound velosity
                hold on
            plot(q_math, cg,'--k','LineWidth',2);

% plot point of spectrum arias
    for i=1:length(qh)
            hold on
        plot(qh(i),W0(i),'or','LineWidth',2); hold on
    end
xlabel(['H=' num2str(round(H_O/1000)) ' km, w_g=' num2str(wg) ' 1/s, T_g=' num2str(round(Tg)) ' min, Lh_1=' num2str(round(Lg)) ' km'],'fontsize',10);
title(['Datafile ' dayOrbit  ',  Orbit Nomber ' num2str(Orbit(1)) ',  UT start ' num2str(UT_NACS(1)/3600) ' hour'],'fontsize',10);


% Errors in spectrum
for i=1:length(qh)
        d1=qh(i)-dqx(i);
        d2=qh(i)+dqx(i);
line([qh(i) d1],[W0(i) W0(i)],'Marker','.','LineStyle','-'); hold on
line([qh(i) d2],[W0(i) W0(i)],'Marker','.','LineStyle','-'); hold on        
        if qh(i)<1
          dw0(i)=dqx(i);
            dw1=W0(i)-dw0(i);
            dw2=W0(i)+dw0(i);
    line([qh(i) qh(i)],[W0(i) dw1],'Marker','.','LineStyle','-'); hold on
    line([qh(i) qh(i)],[W0(i) dw2],'Marker','.','LineStyle','-'); hold on

        elseif qh(i)>1
          Vg(i)=qz(i)/q(i)^2; % km/sec
                
            if real(qz(i))==0
              Vg(i)=1/(1+qx(i)^2)^(3/2); % km/sec
            end
            
          dw0(i)=Vg(i)*dqx(i);
            dw1=W0(i)-dw0(i);
            dw2=W0(i)+dw0(i);
    line([qh(i) qh(i)],[W0(i) dw1],'Marker','.','LineStyle','-'); hold on
    line([qh(i) qh(i)],[W0(i) dw2],'Marker','.','LineStyle','-'); hold on
        end
            Vgrup(i)=(qh(i)*qz(i)/(1+q(i)^2)^(3/2))*Cs; % m/s
end

set(gca,'XLim',[0 6],'YLim',[0 3]);

    filename=[dayOrbit,'_dispers_curve'];
saveas(gcf, fullfile(fpath,filename),'jpeg');
    saveas(gcf, fullfile(fpath,filename),'pdf');


for i=1:3 
    V_faz_goriz(i)=(w(i)*1000)/Kh(i); % m/sec
end

KandQ=[Kx',dqx',Ky',qh',Kz',qz',q'];
LandW=[Lx',w',W0',T_minut',dw0',round(V_faz_goriz)',Vgrup'];
check=[relation_N2_O',check_K',check_V'];