
%fname = ??? [without '.mat']
%cind = ???

Fres = 100;
tstat = 'chi2';

%-------------------------------------------------------------------------------

fpath = fullfile(getenv('DATADIR'),'accu_test',[fname,'.mat']);
load(fpath);

fprintf('c = %f\n',c(cind));

FRc = linspace(min(min(F1c(:,cind)),min(F2c(:,cind))),max(max(F1c(:,cind)),max(F2c(:,cind))),Fres)';

P1ce = empirical_cdf(FRc,F1c(:,cind));
P1ct  = mvgc_cdf(FRc,Finf(cind),1,2*m,1,1,1,0,tstat);

P2ce = empirical_cdf(FRc,F2c(:,cind));
P2ct  = mvgc_cdf(FRc,Fone(cind),1,m,1,1,1,0,tstat);

subplot(1,2,1);
plot(FRc,[P1ct P1ce P2ct P2ce]);
title('causal');
legend('MVGC(t)','MVGC(e)','GCCA(t)','GCCA(e)');

FRn = linspace(min(min(F1n(:,cind)),min(F2n(:,cind))),max(max(F1n(:,cind)),max(F2n(:,cind))),Fres)'/10;

P1ne = empirical_cdf(FRn,F1n(:,cind));
P1nt  = mvgc_cdf(FRn,0,1,1.5*m,1,1,1,0,tstat);

P2ne = empirical_cdf(FRn,F2n(:,cind));
P2nt  = mvgc_cdf(FRn,0,1,m,1,1,1,0,tstat);

subplot(1,2,2);
plot(FRn,[P1nt P1ne P2nt P2ne]);
title('null');
legend('MVGC(t)','MVGC(e)','GCCA(t)','GCCA(e)');
