function [dp_p,FFT_dp_p]=PressureVariation(dO, dN2, O, N2, dT_T)

%kB=1.38e-29; % kg*km2/c2*K 
dp_p=(dO+dN2)./(O+N2)+dT_T;

        p_more=[dp_p; zeros(2^16-length(dp_p),1)]; % increase resolution
    
	FFT_dp_p=fft(p_more,2^16); % Fourier analysis
    
end