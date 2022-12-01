function fit = extract_fit(fit,FitResults,Model,x,y)
    switch Model
        case "biexp"
            fit.f_fast(x,y)  = FitResults{x,y}.a;
            fit.D_slow(x,y)  = FitResults{x,y}.b;
            fit.D_fast(x,y)  = FitResults{x,y}.c;
            fit.S_0(x,y)     = FitResults{x,y}.d;
        case "triexp"
            fit.f_inter(x,y) = FitResults{x,y}.a;
            fit.f_fast(x,y) = FitResults{x,y}.b;
            fit.D_slow(x,y) = FitResults{x,y}.c;
            fit.D_inter(x,y) = FitResults{x,y}.d;
            fit.D_fast(x,y) = FitResults{x,y}.e;
            fit.S_0(x,y) = FitResults{x,y}.f;
    end
end