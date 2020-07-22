
%fname = ??? [without '.mat']
%xname = ???

alpha = 0.05;
mhtc  = 'NONE';
addpath '../../gpmat/';

%-------------------------------------------------------------------------------

fpath = fullfile(getenv('DATADIR'),'accu_test',[fname,'.mat']);
load(fpath);

eval(['x = ' xname ';']);

F1cm = mean(F1c)';
F1nm = mean(F1n)';
F2cm = mean(F2c)';
F2nm = mean(F2n)';

F1cs = std(F1c)';
F1ns = std(F1n)';
F2cs = std(F2c)';
F2ns = std(F2n)';

for i = 1:I;
    for s = 1:S
        sig1c(s,i)  = significance(mvgc_pval(F1c(s,i),mo1(s,i),m,1,1,1,0,'F'),alpha,mhtc);
        sig1n(s,i)  = significance(mvgc_pval(F1n(s,i),mo1(s,i),m,1,1,1,0,'F'),alpha,mhtc);
        sig2c(s,i)  = significance(mvgc_pval(F2c(s,i),mo2(s,i),m,1,1,1,0,'F'),alpha,mhtc);
        sig2n(s,i)  = significance(mvgc_pval(F2n(s,i),mo2(s,i),m,1,1,1,0,'F'),alpha,mhtc);
    end
end

PI1  =   mean(sig1n)';
PII1 = 1-mean(sig1c)';
PI2  =   mean(sig2n)';
PII2 = 1-mean(sig2c)';

gpwrite(sprintf('%s_alpha_%.2f_%s',fpath,alpha,mhtc),[x Finf Fone F1cm F1nm F2cm F2nm F1cs F1ns F2cs F2ns PI1 PII1 PI2 PII2]);

subplot(3,2,1);
plot(x,[F1cm F2cm Finf Fone]);
title('causal mean');
legend('MVGC','GCCA','inf','one');

subplot(3,2,2);
plot(x,[F1nm F2nm]);
title('null mean');
legend('MVGC','GCCA');

subplot(3,2,3);
plot(x,[F1cs F2cs]);
title('causal std. dev.');
legend('MVGC','GCCA');

subplot(3,2,4);
plot(x,[F1ns F2ns]);
title('null std. dev.');
legend('MVGC','GCCA');

subplot(3,2,5);
plot(x,[PII1 PII2]);
title('Type II error rate');
legend('MVGC','GCCA');

subplot(3,2,6);
plot(x,[PI1 PI2]);
title('Type I error rate');
legend('MVGC','GCCA');
