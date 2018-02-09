 for i=1:length(Ut_WATS)
     for j=1:length(Ut_NACS)
         for k=1:length(Ut_WATS)*length(Ut_NACS)
            UT_diff(k)=Ut_WATS(i)-Ut_NACS(j);
         end
     end 
 end
 
 plot(1:length(Ut_WATS),UT_diff,'k','LineWidth',1); grid on