function [fit, FitQuality] = extract_fit_quality(fit, gof, output, sz, Model)
% extract quality paramters for final fit
% sz is the final matrix size

    FitQuality = struct();
    for x = 1:sz(1)
        for y = 1:sz(2) 
            if ~isempty(gof{x,y})
                FitQuality.SSE(x,y) = gof{x,y}.sse;
                FitQuality.Rsq(x,y) = gof{x,y}.rsquare;
                FitQuality.Dfe(x,y) = gof{x,y}.dfe;
                FitQuality.AdjRsq(x,y) = gof{x,y}.adjrsquare;
                FitQuality.RMSE(x,y) = gof{x,y}.rmse;
                FitQuality.Residuals(x,y,:) = output{x,y}.residuals;
                switch Model
                    case {"biexp","Biexp"}
                        fit.f_slow = 1 - fit.f_fast(x,y);
                    case {"triexp","Triexp"}
                        fit.f_slow = 1 - fit.f_interm(x,y) - fit.f_fast(x,y);
                end
            end
        end
    end
end