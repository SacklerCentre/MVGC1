
seed = 8272;
fprintf('seed with %d and save state\n',seed);
s = rng_seed(seed);

fprintf('generated ='); disp(randi(5,1,10));
fprintf('generated ='); disp(randi(5,1,10));

fprintf('restore state\n');
rng_restore(s);

fprintf('generated ='); disp(randi(5,1,10));
fprintf('generated ='); disp(randi(5,1,10));
