%% w determine by quadratic deviation

for i=1:3
        l=length(a_Vz(i):a_Vz(i+1));
        FFT_Vz_kvadratOtklon=abs(FFT_Vz(a_Vz(i):a_Vz(i+1))).^2;
    Vz_dev=sqrt(1/l*sum(FFT_Vz_kvadratOtklon));
    
        l=length(a_dz(i):a_dz(i+1));
        FFT_dz_kvadratOtklon=abs(FFT_dz(a_dz(i):a_dz(i+1))).^2;
    dz_dev=sqrt(1/l*sum(FFT_dz_kvadratOtklon));
    
    w_dev=Vz_dev/dz_dev;
    
    w_norm_by_dev(i)=w_dev/wg;
    
    %% Ky=(y*w*Vy)/(Cs2*dp)
    Ky(i)=(gamma*w_dev*mean(abs(FFT_GW_Vy(a_Vy(i):a_Vy(i+1)))))/(Cs2*mean(abs(FFT_GW_dp(a_dp(i):a_dp(i+1))))); % 1/km
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
        Kz2(i)=(wg^2/w(i)^2-1)*Kh2-Kg^2;
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
        plot(qh(i),w_norm_by_dev(i),'or','LineWidth',2); hold on
    end
xlabel(['H=' num2str(round(H_O/1000)) ' km, w_g=' num2str(wg) ' 1/s, T_g=' num2str(round(Tg)) ' min, Lh_1=' num2str(round(Lg)) ' km'],'fontsize',12);
set(gca,'XLim',[0 4]);
title(['Datafile ' dayOrbit  ',  Orbit Nomber ' num2str(Orbit(1)) ',  UT start ' num2str(UT_NACS(1)/3600) ' hour'],'fontsize',14);

