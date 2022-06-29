function [FitQuality,ROIstat] = IDEALevalBiexp(FitResults,gof,output,DataNii,P,ROIs,ROINii,Data,MaskNii,slice)
%%function [FitQuality,ROIstat] = IDEALevalBiexp(FitResults,gof,output,DataNii,P,ROIs,ROINii,Data,MaskNii,slice)
% organizes Fitting results and plots data
% save workspace at the end 
%
%


% initiate Parameters
f_slow = nan(size(Data(:,:)));
f_fast = nan(size(Data(:,:)));
D_slow = nan(size(Data(:,:)));
D_fast = nan(size(Data(:,:)));
S_0 = nan(size(Data(:,:)));

SSE = nan(size(Data(:,:)));
Rsq = nan(size(Data(:,:)));
Dfe = nan(size(Data(:,:)));
AdjRsq = nan(size(Data(:,:)));
RMSE = nan(size(Data(:,:)));
Residuals = nan(size(Data(:,:), 1), size(Data(:,:), 2), 16);

FitQuality = struct();

% extract Parameters
for x = 1:size(Data, 1)
    for y = 1:size(Data, 2)
        if ~isempty(FitResults{x,y})
            %  fitting Parameters 
            %  d = S_0
            %  a = f_fast
            %  b = D_slow
            %  c = D_fast
            
            f_slow(x, y) = 1 - FitResults{x,y}.a;
            f_fast(x, y) = FitResults{x,y}.a;
            
            D_slow(x, y) = FitResults{x,y}.b;
            D_fast(x, y) = FitResults{x,y}.c;
            S_0(x, y) = FitResults{x,y}.d;
            
            FitQuality.SSE(x, y) = gof{x,y}.sse;
            FitQuality.Rsq(x, y) = gof{x,y}.rsquare;
            FitQuality.Dfe(x, y) = gof{x,y}.dfe;
            FitQuality.AdjRsq(x, y) = gof{x,y}.adjrsquare;
            FitQuality.RMSE(x, y) = gof{x,y}.rmse;
            FitQuality.Residuals(x, y, :) = output{x,y}.residuals;
        end
    end
end

%% Plot the parameter maps
[~,file_name,~] = fileparts(DataNii);
if P.plot
    figure('Visible','on')
    
    subplot(2,3,1)
    imagesc(f_slow);
    caxis(gca,[0 1])
    title('f_{slow}')
    colormap gray;
    axis off;
    
    subplot(2,3,2)
    imagesc(f_fast);
    caxis(gca,[0 1])
    title('f_{fast}')
    colormap gray;
    axis off;
    
    subplot(2,3,3)
    imagesc(D_slow);
    title('D_{slow}')
    colormap gray;
    axis off;
    
    subplot(2,3,4)
    imagesc(D_fast);
    title('D_{fast}')
    colormap gray;
    axis off;
    
    subplot(2,3,5)
    imagesc(S_0);
    title('S_{0,fit}');
    colormap gray;
    axis off;
    
    subplot(2,3,6)
    imagesc(squeeze(Data(:, :, 1)));
    title('S_{0}');
    colormap gray;
    axis off;
    
    if ~exist(P.outputFolder)
        mkdir(P.outputFolder);
    end
    fignm_param = sprintf('%s%sIDEALfit_%s_slice_%d_steps_%d_param.fig',...
        P.outputFolder, filesep, file_name, slice, size(P.Dims_steps,1));
    savefig(gcf, fignm_param);
    close(gcf);
end


ROIstat = eval_rois(ROIs,ROINii,MaskNii,f_slow,f_fast,D_slow,D_fast,S_0);

filenm = sprintf('%s%sIDEALfit_%s_slice_%d_steps_%d.mat',...
    P.outputFolder, filesep, file_name, slice, size(P.Dims_steps, 1));
save(filenm);


end


function ROIstat = eval_rois(ROIs,ROINii,MaskNii,f_slow,f_fast,D_slow,D_fast,S_0)
%% Perform ROI-based analysis
ROIstat.ROIname = cell(1,length(ROIs));
IVIMPars = {'f_slow','f_fast','D_slow','D_fast','S_0'};
    for par = 1 : length(IVIMPars)
        ROIstat.(IVIMPars{par}).mean = zeros(1,length(ROIs));
        ROIstat.(IVIMPars{par}).median = zeros(1,length(ROIs));
        ROIstat.(IVIMPars{par}).std = zeros(1,length(ROIs));
        ROIstat.(IVIMPars{par}).CV = zeros(1,length(ROIs));
        ROIstat.(IVIMPars{par}).iqr = zeros(1,length(ROIs));
        ROIstat.(IVIMPars{par}).q1 = zeros(1,length(ROIs));
        ROIstat.(IVIMPars{par}).q3 = zeros(1,length(ROIs));
    end

    version_num = split(version(),' ');
    version_num = split(version_num{1},'.');
    version_num = join(version_num(1:2),'.');
    version_num = str2num(version_num{1});

    if version_num > 9.4
        for rois = 1 : length(ROIs)
            if rois == 1
                [~,name,~] = fileparts(MaskNii);
            else
                [~,name,~] = fileparts(ROINii{rois-1});
            end
            ROIstat.ROIname{rois} = name;
            for par = 1 : length(IVIMPars)
                eval(['ROIstat.' IVIMPars{par} '.mean(rois) = nanmean(' IVIMPars{par} '(ROIs{rois}==1),''all'');']);
                eval(['ROIstat.' IVIMPars{par} '.median(rois) = nanmedian(' IVIMPars{par} '(ROIs{rois}==1),''all'');']);
                eval(['ROIstat.' IVIMPars{par} '.std(rois) = nanstd(' IVIMPars{par} '(ROIs{rois}==1),0,''all'');']);
                eval(['ROIstat.' IVIMPars{par} '.CV(rois) = ROIstat.' IVIMPars{par} '.std(rois) /ROIstat.' IVIMPars{par} '.mean(rois);']);
                eval(['ROIstat.' IVIMPars{par} '.iqr(rois) = iqr(reshape(' IVIMPars{par} '(ROIs{rois}==1),[],1));']);
                eval(['ROIstat.' IVIMPars{par} '.q1(rois) = prctile(reshape(' IVIMPars{par} '(ROIs{rois}==1),[],1),25);']);
                eval(['ROIstat.' IVIMPars{par} '.q3(rois) = prctile(reshape(' IVIMPars{par} '(ROIs{rois}==1),[],1),1);']);
            end
        end
    elseif version_num <= 9.4
        for rois = 1 : length(ROIs)
            if rois == 1
                [~,name,~] = fileparts(MaskNii);
            else
                [~,name,~] = fileparts(ROINii{rois-1});
            end            
            ROIstat.ROIname{rois} = name;
            for par = 1 : length(IVIMPars)
                eval(['ROIstat.' IVIMPars{par} '.mean(rois) = nanmean(' IVIMPars{par} '(ROIs{rois}==1));']);
                eval(['ROIstat.' IVIMPars{par} '.median(rois) = nanmedian(' IVIMPars{par} '(ROIs{rois}==1));']);
                eval(['ROIstat.' IVIMPars{par} '.std(rois) = nanstd(' IVIMPars{par} '(ROIs{rois}==1),0);']);
                eval(['ROIstat.' IVIMPars{par} '.CV(rois) = ROIstat.' IVIMPars{par} '.std(rois) /ROIstat.' IVIMPars{par} '.mean(rois);']);
                eval(['ROIstat.' IVIMPars{par} '.iqr(rois) = iqr(reshape(' IVIMPars{par} '(ROIs{rois}==1),[],1));']);
                eval(['ROIstat.' IVIMPars{par} '.q1(rois) = prctile(reshape(' IVIMPars{par} '(ROIs{rois}==1),[],1),25);']);
                eval(['ROIstat.' IVIMPars{par} '.q3(rois) = prctile(reshape(' IVIMPars{par} '(ROIs{rois}==1),[],1),1);']);
            end
        end
    end
end