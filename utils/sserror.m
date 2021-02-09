function err = sserror(rep,y,tol)

if nargin < 3 || isempty(tol), tol = sqrt(eps); end % DARE relative residual tolerance

err = rep < 0;   % DARE failed: show-stopper!

if nargin < 2 || isempty(y)
    if err
        fprintf(2,'DARE ERROR: ');
        switch rep
            case -1, fprintf(2,'eigenvalues on/near unit circle\n');
            case -2, fprintf(2,'no stablising solution\n');
        end
        return
    end
    if rep > tol
        fprintf(2,'DARE WARNING: large relative residual = %e\n',rep);
    end
else
    if err
        fprintf(2,'DARE ERROR for source %d: ',y);
        switch rep
            case -1, fprintf(2,'eigenvalues on/near unit circle\n');
            case -2, fprintf(2,'no stablising solution\n');
        end
        return
    end
    if rep > tol
        fprintf(2,'DARE WARNING for source %d: large relative residual = %e\n',y,rep);
    end
end
