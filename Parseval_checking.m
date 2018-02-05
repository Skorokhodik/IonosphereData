%% Checking out of Parseval`s theorem

% O
for i=1:3
    middle(i)=a(i)+(a(i+1)-a(i))/2;
            Lx(i)=2^16*7.8/middle(i);
                l2(i)=2^16*7.8/a(i+1);
                l1(i)=2^16*7.8/a(i);
            
            % Kx=2*pi/Lx;
                Kx(i)=2*pi/Lx(i); % 1/km
                    qx(i)=Kx(i)/Kg;
                        dqx(i)=(1/l2(i)-1/l1(i))*2*pi/Kg;

[FFT_GWO_Aria, iFFT_GWO_Aria_norm]=GravitationWave_Aria(FFT_GW_O, Trend_O, L_Ox, a(i), a(i+1));
            hold on
        subplot(3,3,i), plot(UT_NACS_1sec./3600,iFFT_GWO_Aria_norm(1:L_Ox),'r','LineWidth',2); grid on
                    xlabel(['L_x=' num2str(round(Lx(i))) ' km, O (red line), N_2 (green)'],'fontsize',12);
                    set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
end
            title(['Datafile ' dayOrbit  ', Orbit ' num2str(Orbit(1)) ', UT start ' num2str(UT_NACS(1)/3600) 'hour'],'fontsize',14);

% N2
            for i=1:3
    [FFT_GWN2_Aria, iFFT_GWN2_Aria_norm]=GravitationWave_Aria(FFT_GW_N2, Trend_N2, L_Ox, a(i), a(i+1));
            hold on
        subplot(3,3,i), plot(UT_NACS_1sec./3600,iFFT_GWN2_Aria_norm(1:L_Ox),'g','LineWidth',2);

            end

% else
for i=1:3            
    % dz & Vz
[FFT_dz_GW_Aria, iFFT_dz_Aria]=GravitationWave_Aria_Wind(FFT_dz, a_dz(i), a_dz(i+1)); % [m]
                hold on
            subplot(3,3,3+i), plot(UT_NACS_1sec./3600,iFFT_dz_Aria(1:L_Ox).*1e-3,'k','LineWidth',2); grid on
[FFT_Vz_GW_Aria, iFFT_Vz_GW_Aria]=GravitationWave_Aria_Wind(FFT_Vz, a_Vz(i), a_Vz(i+1)); % [m/s]
                hold on
            subplot(3,3,3+i), plot(UT_WATS_Vz_1sec./3600,iFFT_Vz_GW_Aria(1:L_Vz),'b','LineWidth',2);
                    xlabel('V_z (blue), dz (black)','fontsize',12); 
                    set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);
	% Vy & dp/p
[FFT_Vy_GW_Aria, iFFT_Vy_GW_Aria]=GravitationWave_Aria_Wind(FFT_Vy, a_Vy(i), a_Vy(i+1)); % [m/s]
                hold on
            subplot(3,3,6+i), plot(UT_WATS_Vy_1sec./3600,iFFT_Vy_GW_Aria(1:L_Vy)./10,'m','LineWidth',2); grid on
                hold on
            subplot(3,3,6+i), plot(UT_WATS_Vz_1sec./3600,iFFT_Vz_GW_Aria(1:L_Vz)./10,'b','LineWidth',2);
[FFT_dp_GW_Aria,iFFT_dp_GW_Aria]=GravitationWave_Aria_Wind(FFT_GW_dp, a_dp(i), a_dp(i+1));
                hold on
%             subplot(3,3,6+i), plot(1:length(dp_p),iFFT_dp_GW_Aria(1:length(dp_p)),'c','LineWidth',2);
                    xlabel('V_y (magenta), V_z (blue)','fontsize',12); 
                    set(gca,'XLim',[UT_NACS_1sec(1)/3600 UT_NACS_1sec(end)/3600]);

%% w=d|Vz(k)|/d|z(k)|
            Vz_Parseval_space=sqrt(1/L_Vz*sum((abs(iFFT_Vz_GW_Aria(1:L_Vz))).^2));
            dz_Parseval_space=sqrt(1/L_Ox*sum((abs(iFFT_dz_Aria(1:L_Ox))).^2));
            
w_Pars_space(i)=Vz_Parseval_space/dz_Parseval_space; % 1/s
    
    
                T_minut_Pars_space(i)=2*pi/(w_Pars_space(i)*60); % min
        W0_Pars_space(i)=w_Pars_space(i)/wg;

%% Ky=(y*w*Vy)/(Cs2*dp)
    Ky(i)=(gamma*w_Pars_space(i)*mean(abs(FFT_GW_Vy(a_Vy(i):a_Vy(i+1)))))/(Cs2*mean(abs(FFT_GW_dp(a_dp(i):a_dp(i+1))))); % 1/km
        qy(i)=Ky(i)/Kg;

%% Ly=2*pi/Ky
    Ly(i)=2*pi/Ky(i); % km

%% Kh2=Kx^2+Ky^2
        Kh2=Kx(i)^2+Ky(i)^2; % horisontal wave nomber, 1/km2
	Kh(i)=sqrt(Kh2);
        qh(i)=Kh(i)/Kg;

%% Lh=2*pi/Kh
    Lh(i)=2*pi/Kh(i);

%% Kz2=(wg^2/w^2-1)*Kh2-Kg2
        Kz2(i)=(wg^2/w_Pars_space(i)^2-1)*Kh2-Kg^2;
    Kz(i)=sqrt(Kz2(i));
        qz(i)=Kz(i)/Kg;

%% Lz=2*pi/Kz
    Lz(i)=2*pi/Kz(i);
                            
%% K=sqrt(Kh2+Kz2)
    K(i)=sqrt(Kh2+Kz2(i));
        q(i)=K(i)/Kg;
end        
        
figure % DISPERTION ARIA
    % create mathematic function
        q_math=0.01:0.001:4;
        W=q_math./sqrt(1+q_math.^2); % upper border condition
            plot(q_math, W,'--k','LineWidth',2); grid on
                hold on
        Wdiss=0.1*q_math.^(2/3); % dissipative bottom border condition
                hold on
            plot(q_math, Wdiss,'--k','LineWidth',2);
                hold on
                
	% plot point of spectrum arias
    for i=1:3
            hold on
        plot(qh(i),W0_Pars_space(i),'or','LineWidth',2); hold on
    end
xlabel(['H=' num2str(round(H_O/1000)) ' km, w_g=' num2str(wg) ' 1/s, T_g=' num2str(round(Tg)) ' min, Lh_1=' num2str(round(Lg)) ' km'],'fontsize',12);
set(gca,'XLim',[0 4]);
title(['Datafile ' dayOrbit  ',  Orbit Nomber ' num2str(Orbit(1)) ',  UT start ' num2str(UT_NACS(1)/3600) ' hour'],'fontsize',14);
