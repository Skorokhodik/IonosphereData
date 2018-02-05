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

            % Whole GW spectrum
%             WholeGW_2_zeros=[iFFT_GW_O_1; zeros(2^16-length(iFFT_GW_O_1), 1)];
%             FFT_GW_1=fft(WholeGW_1_zeros,NFFT);
% a - space interval from the file maximum.m
aria1=a(1);
aria2=a(2);

% O_scale aria and fft
Oxigen_aria=iFFT_GW_O_2(aria1:aria2);
Ox_aria_zero=[Oxigen_aria ; zeros(2^16-length(Oxigen_aria), 1)];

FFT_Ox_aria=fft(Ox_aria_zero, NFFT);
            
        subplot(211), plot(aria1:aria2,Oxigen_aria,'r','LineWidth',1); grid on
            set(gca,'XLim',[0 L_Ox]);
        subplot(212), plot(0:NFFT-1, abs(FFT_Ox_aria),'r','LineWidth',2); grid on
        set(gca,'XLim',[0 3000]);
        % line of rms in spectrum
            for i=1:3000
                line_rms_height=rms(abs(FFT_Ox_aria(1:3000)));
                line_rms(i)=line_rms_height;
            end
            hold on, subplot(212), plot(0:length(line_rms)-1, line_rms./2,'k','LineWidth',1);
            legend('fft dO','line of rms(dO)/2');
            title(['Datafile ' dayOrbit  ', Orbit ' num2str(Orbit(1)) ', UT start ' num2str(UT_NACS(1)/3600) 'hour'],...
                    'fontsize',10);


           aria1_spectrum=966;
           aria2_spectrum=1204;
                    
                middle_O_aria=(aria1_spectrum+aria2_spectrum)/2;
                Lx_aria=2^16*7.8/middle_O_aria; % km
                    l2_aria=2^16*7.8/aria2_spectrum;
                    l1_aria=2^16*7.8/aria1_spectrum;
                        
% dz aria and fft
dz_aria=iFFT_GW_dz_2(aria1:aria2);
dz_aria_zero=[dz_aria ; zeros(2^16-length(dz_aria), 1)];

FFT_dz_aria=fft(dz_aria_zero, NFFT);

% Vz aria and fft
Vz_aria=iFFT_GW_Vz_2(aria1:aria2);
Vz_aria_zero=[Vz_aria ; zeros(2^16-length(Vz_aria), 1)];

FFT_Vz_aria=fft(Vz_aria_zero, NFFT);


%% w=dVz(k)/dz(k)
    w_aria=rms(abs(FFT_Vz_aria))/rms(abs(FFT_dz_aria)); % 1/s
                T_aria_minut=2*pi/(w_aria*60); % min
        W0_aria=w_aria/wg;

%% Kx
            
            % Kx=2*pi/Lx;
    Kx_aria=2*pi/Lx_aria; % 1/km
    	qx_aria=Kx_aria/Kg;
                        
                        dKx_aria=(1/l2_aria-1/l1_aria)*pi;
                        
                        dqx_aria=dKx_aria/Kg;


%% Ky=(y*w*Vy)/(Cs2*dp)
    Ky_aria=(gamma/Cs2)*w_aria*1000*(rms(abs(iFFT_GW_Vy_2(aria1:aria2)))/(rms(abs(iFFT_GW_O_2(aria1:aria2)))...
        +rms(abs(iFFT_GW_T_2(aria1:end))))); % 1/km
        qy_aria=Ky_aria/Kg;

%% Ly=2*pi/Ky
    Ly_aria=2*pi/Ky_aria; % km

%% Kh2=Kx^2+Ky^2
        Kh2_aria=Kx_aria^2+Ky_aria^2; % horisontal wave nomber, 1/km2
	Kh_aria=sqrt(Kh2_aria);
        qh_aria=Kh_aria/Kg;

%% Lh=2*pi/Kh
    Lh_aria=2*pi/Kh_aria;

%% Kz2=(wg^2/w^2-1)*Kh2-Kg2
        Kz2_aria=(wg^2/w_aria^2-1)*Kh2_aria-Kg^2;
    Kz_aria=sqrt(Kz2_aria);
        qz_aria=Kz_aria/Kg;

%% Lz=2*pi/Kz
    Lz_aria=2*pi/Kz_aria;
                            
%% K=sqrt(Kh2+Kz2)
    K_aria=sqrt(Kh2_aria+Kz2_aria);
        q_aria=K_aria/Kg;
        
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

% plot point of spectrum aria
    hold on
plot(qh_aria,W0_aria,'or','LineWidth',2); hold on
    xlabel(['H=' num2str(round(H_O/1000)) ' km, w_g=' num2str(wg) ' 1/s, T_g=' num2str(round(Tg))...
        ' min, Lh_1=' num2str(round(Lg)) ' km'],'fontsize',10);
    title(['Datafile ' dayOrbit  ',  Orbit Nomber ' num2str(Orbit(1)) ',  UT start ' num2str(UT_NACS(1)/3600)...
        ' hour'],'fontsize',10);


% Errors in spectrum
        d1=qh_aria-dqx_aria;
        d2=qh_aria+dqx_aria;
line([qh_aria d1],[W0_aria W0_aria],'Marker','.','LineStyle','-'); hold on
line([qh_aria d2],[W0_aria W0_aria],'Marker','.','LineStyle','-'); hold on        
        if qh_aria<1
          dw0_aria=dqx_aria;
            dw1=W0_aria-dw0_aria;
            dw2=W0_aria+dw0_aria;
    line([qh_aria qh_aria],[W0_aria dw1],'Marker','.','LineStyle','-'); hold on
    line([qh_aria qh_aria],[W0_aria dw2],'Marker','.','LineStyle','-'); hold on

        elseif qh_aria>1
          Vg_aria=qz_aria/q_aria^2; % km/sec
                
            if real(qz_aria)==0
              Vg_aria=1/(1+qx_aria^2)^(3/2); % km/sec
            end
            
          dw0_aria=Vg_aria*dqx_aria;
            dw1=W0_aria-dw0_aria;
            dw2=W0_aria+dw0_aria;
    line([qh_aria qh_aria],[W0_aria dw1],'Marker','.','LineStyle','-'); hold on
    line([qh_aria qh_aria],[W0_aria dw2],'Marker','.','LineStyle','-'); hold on
        end
        
Vgrup_aria=(qh_aria*qz_aria/(1+q_aria^2)^(3/2))*Cs; % m/s

set(gca,'XLim',[0 6],'YLim',[0 3]);

 V_faz_goriz_aria=(w_aria*1000)/Kh_aria;