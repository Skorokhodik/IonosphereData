clear
    figure
%% create mathematic function for method proba

q=0.01:0.001:4;

W=q./sqrt(1+q.^2);

plot(q,W,'g','LineWidth',2); grid on
    hold on
    

Wdiss=0.1*q.^(2/3);
        hold on
plot(q,Wdiss,'g','LineWidth',2);
        hold on
        
%% POINTs

    W0=[.17 .19 .21 .21 .53 .59 .75 .77 .51 .61 .6 .61];
    qh=[.49 1.34 2.41 3.41 .68 1.17 1.87 2.59 1.01 1.4 1.78 2.78];
        dqx=[.25 .6 .46 .54 .23 .27 .43 .28 .2 .19 .19 .47];
    qz=[2.69 6.72 11.29 15.76 .39 1.24 1.29 1.85 1.36 1.53 2.17 3];
    q=[2.73 6.85 11.55 16.12 .78 1.7 2.28 3.18 1.69 2.07 2.81 3.86];


    
for i=1:length(W0)
    
        hold on
plot(qh(i),W0(i),'or','LineWidth',2);
        hold on
        
        d1=qh(i)-dqx(i);
        d2=qh(i)+dqx(i);
line([qh(i) d1],[W0(i) W0(i)],'Marker','.','LineStyle','-'); hold on
line([qh(i) d2],[W0(i) W0(i)],'Marker','.','LineStyle','-'); hold on        
        if qh(i)<1
          dw(i)=dqx(i);
            dw1=W0(i)-dw(i);
            dw2=W0(i)+dw(i);
    line([qh(i) qh(i)],[W0(i) dw1],'Marker','.','LineStyle','-'); hold on
    line([qh(i) qh(i)],[W0(i) dw2],'Marker','.','LineStyle','-'); hold on
set(gca,'YLim',[0 1.5]);
        elseif qh(i)>1
          Vg=qz(i)/q(i)^2;
          dw(i)=Vg*dqx(i);
            dw1=W0(i)-dw(i);
            dw2=W0(i)+dw(i);
    line([qh(i) qh(i)],[W0(i) dw1],'Marker','.','LineStyle','-'); hold on
    line([qh(i) qh(i)],[W0(i) dw2],'Marker','.','LineStyle','-'); hold on
        end
end

line([0 4],[1 1],'Marker','.','LineStyle','--','LineWidth',1,'Color','k'); hold on
