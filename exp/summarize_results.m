function [fit, FitQuality] = summarize_results(fit, gof, output, sz, Model)
%summarize_results(fit 'struct', gof 'struct', output 'struct', sz 'double')
%'double', Model 'string')
% extract quality paramters for final fit and summarize them in FitQuality
% also calculates missing fraction from fitted ones
%
%   fit     : struct containing fitting results
%   gof     : fitting results
%   output  : fitting output struct
%   sz      : is the final matrix size
%   Model   : "Biexp" or "Triexp"

    FitQuality = struct();
    fit.f_slow = NaN*ones(sz);
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
                        fit.f_slow(x,y) = 1 - fit.f_fast(x,y);
                    case {"triexp","Triexp"}
                        fit.f_slow(x,y) = 1 - fit.f_inter(x,y) - fit.f_fast(x,y);
                end
            end
        end
    end
end