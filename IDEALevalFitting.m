function [FitQuality,ROIstat] = IDEALevalFitting(FitResults,gof,output,P,Data,DataNii,ROIs,ROINii,MaskNii,slice)


% initiate Parameters
if strcmp(P.Model,'Biexp') || strcmp(P.Model,'Triexp')
    f_slow = nan(size(Data(:,:,1)));
    f_fast = nan(size(Data(:,:,1)));
    D_slow = nan(size(Data(:,:,1)));
    D_fast = nan(size(Data(:,:,1)));
    S_0 = nan(size(Data(:,:,1)));
end
if strcmp(P.Model,'Triexp')
    f_interm = nan(size(Data(:,:,1)));
    D_interm = nan(size(Data(:,:,1)));
end

SSE = nan(size(Data(:,:,1)));
Rsq = nan(size(Data(:,:,1)));
Dfe = nan(size(Data(:,:,1)));
AdjRsq = nan(size(Data(:,:,1)));
RMSE = nan(size(Data(:,:,1)));
Residuals = nan(size(Data(:,:,1), 1), size(Data(:,:,1), 2), 16);


FitQuality = struct();

% extract fitted parameters from fitting-outputs

for x = 1:size(Data, 1)
    for y = 1:size(Data, 2)
        if ~isempty(FitResults{x,y})
            if strcmp(P.Model,'Biexp')
                %  fitting Parameters for Biexp model
                %  d = S_0
                %  a = f_fast
                %  b = D_slow
                %  c = D_fast

                f_slow(x, y) = 1 - FitResults{x,y}.a;
                f_fast(x, y) = FitResults{x,y}.a;

                D_slow(x, y) = FitResults{x,y}.b;
                D_fast(x, y) = FitResults{x,y}.c;
                S_0(x, y) = FitResults{x,y}.d;
            end
            if strcmp(P.Model,'Triexp')
                %  fitting Parameters for Triexp model
                %  f = S_0
                %  a = f_inter
                %  b = f_fast
                %  c = D_slow
                %  d = D_inter
                %  e = D_fast
                f_slow(x, y) = 1 - FitResults{x,y}.a - FitResults{x,y}.b;
                f_interm(x, y) = FitResults{x,y}.a;
                f_fast(x, y) = FitResults{x,y}.b;

                D_slow(x, y) = FitResults{x,y}.c;
                D_interm(x, y) = FitResults{x,y}.d;
                D_fast(x, y) = FitResults{x,y}.e;
                S_0(x, y) = FitResults{x,y}.f;
            end
            
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
    switch P.Model
        case 'Biexp'
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
            fignm_param = sprintf('%s%sIDEALfit_%s_%s_slice_%d_steps_%d_param.fig',...
                P.outputFolder, filesep, file_name, P.Model, slice, size(P.Dims_steps,1));
            savefig(gcf, fignm_param);
            close(gcf);            
            ROIstat = eval_rois_biexp(ROIs,ROINii,MaskNii,f_slow,f_fast,...
                        D_slow,D_fast,S_0);
        
        
        case 'Triexp'
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
            fignm_param = sprintf('%s%sIDEALfit_%s_%s_slice_%d_steps_%d_param.fig',...
                P.outputFolder, filesep, file_name, P.Model, slice, size(P.Dims_steps,1));
            savefig(gcf, fignm_param);
        %     close(gcf);  
            ROIstat = eval_rois_triexp(ROIs,ROINii,MaskNii,f_slow,f_interm,...
                                f_fast,D_slow,D_interm,D_fast,S_0);
    end
end

filenm = sprintf('%s%sIDEALfit_%s_%s_slice_%d_steps_%d.mat',...
    P.outputFolder, filesep, file_name, P.Model, slice, size(P.Dims_steps, 1));
save(filenm);

end

function ROIstat = eval_rois_biexp(ROIs,ROINii,MaskNii,f_slow,f_fast,D_slow,D_fast,S_0)
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

function ROIstat = eval_rois_triexp(ROIs,ROINii,MaskNii,f_slow,f_interm,f_fast,D_slow,D_interm,D_fast,S_0)
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