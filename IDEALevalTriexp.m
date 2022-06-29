function [FitQuality,ROIstat] = IDEALevalTriexp(Mask,FitResults,gof,output,DataNii,P,ROIs,ROINii,Data,MaskNii,slice)

%% Extract the parameters maps

f_slow = nan(size(Mask));
f_interm = nan(size(Mask));
f_fast = nan(size(Mask));
D_slow = nan(size(Mask));
D_interm = nan(size(Mask));
D_fast = nan(size(Mask));
S_0 = nan(size(Mask));

SSE = nan(size(Mask));
Rsq = nan(size(Mask));
Dfe = nan(size(Mask));
AdjRsq = nan(size(Mask));
RMSE = nan(size(Mask));
Residuals = nan(size(Mask, 1), size(Mask, 2), 16);

FitQuality = struct();

for kx = 1:size(Mask, 1)
    for ky = 1:size(Mask, 2)
        if ~isempty(FitResults{kx,ky})
            f_slow(kx, ky) = 1 - FitResults{kx,ky}.a - FitResults{kx,ky}.b;
            f_interm(kx, ky) = FitResults{kx,ky}.a;
            f_fast(kx, ky) = FitResults{kx,ky}.b;
            
            D_slow(kx, ky) = FitResults{kx,ky}.c;
            D_interm(kx, ky) = FitResults{kx,ky}.d;
            D_fast(kx, ky) = FitResults{kx,ky}.e;
            S_0(kx, ky) = FitResults{kx,ky}.f;
            
            FitQuality.SSE(kx, ky) = gof{kx,ky}.sse;
            FitQuality.Rsq(kx, ky) = gof{kx,ky}.rsquare;
            FitQuality.Dfe(kx, ky) = gof{kx,ky}.dfe;
            FitQuality.AdjRsq(kx, ky) = gof{kx,ky}.adjrsquare;
            FitQuality.RMSE(kx, ky) = gof{kx,ky}.rmse;
            FitQuality.Residuals(kx, ky, :) = output{kx,ky}.residuals;
        end
    end
end

%% Plot the parameter maps
[~,file_name,~] = fileparts(DataNii);
if P.plot
    figure('Visible','on')
    
    subplot(3,3,1)
    imagesc(f_slow);
    caxis(gca,[0 1])
    title('f_{slow}')
    colormap gray;
    axis off;
    
    subplot(3,3,2)
    imagesc(f_interm);
    caxis(gca,[0 1])
    title('f_{interm}')
    colormap gray;
    axis off;
    
    subplot(3,3,3)
    imagesc(f_fast);
    caxis(gca,[0 1])
    title('f_{fast}')
    colormap gray;
    axis off;
    
    subplot(3,3,4)
    imagesc(D_slow);
    title('D_{slow}')
    colormap gray;
    axis off;
    
    subplot(3,3,5)
    imagesc(D_interm);
    title('D_{interm}')
    colormap gray;
    axis off;
    
    subplot(3,3,6)
    imagesc(D_fast);
    title('D_{fast}')
    colormap gray;
    axis off;
    
    subplot(3,3,7)
    imagesc(S_0);
    title('S_{0,fit}');
    colormap gray;
    axis off;
    
    subplot(3,3,8)
    imagesc(squeeze(Data(:, :, 1)));
    title('S_{0}');
    colormap gray;
    axis off;
    
    if ~exist(P.outputFolder)
        mkdir(P.outputFolder);
    end
    fignm_param = sprintf('%s%sIDEALfit_%s_slice_%d_steps_%d_param.fig',...
        P.outputFolder, filesep, slice ,file_name, size(P.Dims_steps,1));
    savefig(gcf, fignm_param);
    close(gcf);
end
   
ROIstat = eval_rois(ROIs,ROINii,MaskNii,f_slow,f_interm,f_fast,D_slow,D_interm,D_fast,S_0);

filenm = sprintf('%s%sIDEALfit_%s_slice_%d_steps_%d.mat',...
    P.outputFolder, filesep,slice, file_name, size(P.Dims_steps, 1));
save(filenm);
end
function ROIstat = eval_rois(ROIs,ROINii,MaskNii,f_slow,f_interm,f_fast,D_slow,D_interm,D_fast,S_0)
%% Perform ROI-based analysis
ROIstat.ROIname = cell(1,length(ROIs));
IVIMPars = {'f_slow','f_interm','f_fast','D_slow','D_interm','D_fast','S_0'};
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