clear
%% filtr Fedorenko 
f=1:0.05:2600;

        Lave1=500000; % 700 km (1point = 7.8km)
%         Lave2=557*2; % 70 km
 
for i=1:length(f)
    H(i)=1-297*sin(f(i)/297)/f(i); % A4Kh of filtr
end

figure
    subplot(312); plot(f,abs(H),'g','LineWidth',2); 
        set(gca,'XLim',[0 max(f)]); grid on
