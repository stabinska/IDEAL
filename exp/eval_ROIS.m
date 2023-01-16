function ROIstat = eval_ROIS(ROIs,fit,Model)
%eval_ROIs(ROIs 'cell', fit 'struct', Model 'string')
% evaluate ROIs for calculated parameters
%
%   ROIs    : cell array containing ROIs as matrix
%   fit     : struct containing fitting results
%   Model   : "Biexp" or "Triexp"


    switch Model
        case {"biexp","Biexp"}
            params = ["f_fast";"D_slow";"D_fast";"S_0"];
        case {"triexp","Triexp"}
            params = ["f_inter";"f_fast";"D_slow";...
                "D_inter";"D_fast";"S_0"];
    end
       
    for idx_roi = 1:length(ROIs)
        for idx_IDEAL = 1:length(params)
            ROIstat.(params(idx_IDEAL)).mean(idx_roi) = ...
                mean(fit.(params(idx_IDEAL))(ROIs{idx_roi}==1),"all","omitnan");
            ROIstat.(params(idx_IDEAL)).median(idx_roi) = ...
                median(fit.(params(idx_IDEAL))(ROIs{idx_roi}==1),"all","omitnan");
            ROIstat.(params(idx_IDEAL)).std(idx_roi) = ...
                std(fit.(params(idx_IDEAL))(ROIs{idx_roi}==1),0,"all","omitnan");
            ROIstat.(params(idx_IDEAL)).CV(idx_roi) = ...
                ROIstat.(params(idx_IDEAL)).std(idx_roi)/...
                ROIstat.(params(idx_IDEAL)).mean(idx_roi);
            ROIstat.(params(idx_IDEAL)).iqr(idx_roi) = ...
                iqr(reshape(fit.(params(idx_IDEAL))(ROIs{idx_roi}==1),[],1));
            ROIstat.(params(idx_IDEAL)).q1(idx_roi) = ...
                prctile(reshape(fit.(params(idx_IDEAL))(ROIs{idx_roi}==1),[],1),25);
            ROIstat.(params(idx_IDEAL)).q3(idx_roi) = ...
                prctile(reshape(fit.(params(idx_IDEAL))(ROIs{idx_roi}==1),[],1),1);
        end
    end
end