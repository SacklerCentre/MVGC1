n = 5;
p = 10;

rho = 0.99;

aclags = 50;

x = 1:2;

%----------------------------

[A,oldrho]= var_specrad(randn(n,n,p),rho);

V = eye(n);

ptic('\n*** var_to_autocov... ');
[G,info] = var_to_autocov_test(A,V,aclags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,false); % report results (and bail out on error)

%figure(1); clf;
%plot_autocov(G,[],[],[],false,false);

ptic('*** autocov_to_var... ');
[A1,V1,DV1] = autocov_to_var_test(G);
ptoc;

for k = 1:aclags
    fprintf('k = %d : DV1 = %g\n',k,DV1(k));
end

AA = reshape(cat(3,A,zeros(n,n,aclags-p)),n,n*aclags);
AA1 = reshape(A1,n,n*aclags);
DA = norm(AA-AA1)/norm(AA);
DV = norm(V-V1)/norm(V);

fprintf('\nDA = %g\n',DA);
fprintf('eDV = %g\n\n',DV);

ptic('*** autocov_to_var... ');
[A2,V2,DV2] = autocov_to_var_test(G(x,x,:));
ptoc;

for k = 1:aclags
    fprintf('k = %d : DV2 = %g\n',k,DV2(k));
end

figure(1); clf;
mo = 1:aclags;
dec = exp(-rho*mo);
semilogy(mo,abs([DV1;DV2;dec]));
legend('full','reduced','exp');
