function fit = extract_fit(fit,FitResults,Model,x,y)
%extract_fit(fit 'struct', FitResults 'struct', Model 'string', x 'double', y 'double')
%transform fit results to "human readable" variables 
%
%   fit         : struct containing fitting results
%   FitResults  : fitting results
%   Model       : "Biexp", "Biexp_T1corr" or "Triexp"
%   x,y         : current x and y coordinate

    switch Model
        case {"S0fit"}
            fit.S_0(x,y)     = FitResults{x,y}.a;
            fit.T1(x,y)      = FitResults{x,y}.b;
            fit.T2(x,y)      = FitResults{x,y}.c;
        case {"biexp", "Biexp"}
            fit.f_fast(x,y)  = FitResults{x,y}.a;
            fit.D_slow(x,y)  = FitResults{x,y}.b;
            fit.D_fast(x,y)  = FitResults{x,y}.c;
            fit.S_0(x,y)     = FitResults{x,y}.d;    
        case {"biexp_T1corr","Biexp_T1corr"}
            fit.f_fast(x,y)  = FitResults{x,y}.a;
            fit.D_slow(x,y)  = FitResults{x,y}.b;
            fit.D_fast(x,y)  = FitResults{x,y}.c;
            fit.S_0(x,y)     = FitResults{x,y}.d;
            fit.T1(x,y)      = FitResults{x,y}.e;
        case {"triexp", "Triexp"}
            fit.f_inter(x,y) = FitResults{x,y}.a;
            fit.f_fast(x,y) = FitResults{x,y}.b;
            fit.D_slow(x,y) = FitResults{x,y}.c;
            fit.D_inter(x,y) = FitResults{x,y}.d;
            fit.D_fast(x,y) = FitResults{x,y}.e;
            fit.S_0(x,y) = FitResults{x,y}.f;
    end
end