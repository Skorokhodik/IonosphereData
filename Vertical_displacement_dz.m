function [dz, FFT_dz]=Vertical_displacement_dz(Temperature, dO_O, dN2_N2, L)

T=mean(Temperature); % Temperature
%% dz
	 %constants
        R=8.32*1000; % gramm*m^2/mol*K*c^2
        m_O=16; %gramm/mol
        m_N2=28;%gramm/mol
        g=9.8; % m/s2;
        
        H_O=R*T/(m_O*g); %[m]
        H_N2=R*T/(m_N2*g); %[m]
        
dz_all=[(H_O*H_N2/(H_N2-H_O))*(dO_O-dN2_N2); zeros(2^16-L,1)]; %[m]
dz=dz_all(1:L);

        NFFT=2^16;
        FFT_dz=fft(dz_all,NFFT);



end