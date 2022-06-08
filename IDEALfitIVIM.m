function [FitResults,FitQuality,P,ROIstat] = IDEALfitIVIM(DataNii,P,MaskNii,ROINii)

%%  Function IDEALfit
% 
%   Fit a bi-exponential or tri-exponential model to DWI data using the
%   Iterative Downsampling Adaptive Least-squares (IDEAL) approach
%  
%
%   Input:    Data(x_pos, y_pos, z_pos, b-value)   - 4D-DWI data (*.nii.gz)
%             P                                    - structure with experiment and fit parameters
%             MaskNii                              - mask (e.g. both kidneys)  (*.nii.gz)
%             ROINii                               - regions of interest (e.g. cortex, medulla, pelvis) (*.nii.gz) as cell 
%             
%
%   Output:   FitResults          - structure with fit results
%             FitQuality          - structure with goodnes-of-fit measures (SSE, R2, adjR2, stdError, CVs)
%             P                   - updated structure with experiment and fit parameters
%
%
%   Authors: Julia Stabinska (jstabin3@jhmi.edu)
%            Helge Jörn Zöllner (hzoelln2@jhmi.edu)
%            Thomas Andreas Thiel (thomas.thiel@hhu.de)
%
%% Parse input arguments
if nargin < 3
    error('You need to supply at least a 3D-DWI data matrix as *.nii.gz, a struct with settings (see IDEAL script), and a mask volume of the volume of interest as *.nii.gz');
end

%% Load data, mask ,and ROIs
data_nifti = double(rot90(niftiread(DataNii)));
Data = squeeze(data_nifti(:,:,P.slice,:));

Mask_nifti = double(rot90(niftiread(MaskNii)));
Mask = Mask_nifti(:,:,P.slice); 
Mask_nan = Mask;
Mask_nan(Mask == 0) = NaN;

if nargin > 3
    ROIs{1} = double(rot90(niftiread(MaskNii)));
    ROIs{1} = ROIs{1}(:,:,P.slice); 
    ROIs{1}(ROIs{1}==0) =NaN; 
    for rois = 2 : length(ROINii)+1
        ROIs{rois} = double(rot90(niftiread(ROINii{rois-1})));
        ROIs{rois}(ROIs{rois}==0) =NaN;
    end
else
   ROIs{1} = double(rot90(niftiread(MaskNii)));
   ROIs{1} = ROIs{1}(:,:,P.slice); 
   ROIs{1}(ROIs{1}==0) =NaN; 
   sprintf('No ROI supplied. We will do statistics on the VOI mask.')
end

%% Perform IDEAL fitting

clearvars a b c d e f fitresults gof output
Data_masked = Data.*Mask_nan;

tStart = tic;
for res = 1 : size(P.Dims_steps, 1)
    
    fprintf('Downsampling step no: %s\n', num2str(res));
    
        a = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2));
        b = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2));
        c = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2));
        d = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2));
        e = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2));
        f = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2));
        
    if res < size(P.Dims_steps, 1)
        Mask_res = imresize(Mask, [P.Dims_steps(res,2) P.Dims_steps(res,1)],'bilinear');
        Data_res = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2),16);  
        
        for i = 1:numel(P.b_values) 
             Data_res(:,:,i) = imresize(squeeze(Data(:, :, i)),[P.Dims_steps(res,2) P.Dims_steps(res,1)],'bilinear');               
        end
        
        Mask_res = abs(Mask_res);
        Mask_res(Mask_res < 0.025) = 0;

    else
        Data_res = Data_masked;
        Mask_res = Mask;
    end
    
    FitResults = cell(size(Mask_res, 1), size(Mask_res, 2));
    gof = FitResults;
    output = FitResults;
    
    for x = 1 : size(Mask_res, 1)
        for y = 1: size(Mask_res, 2)
                op.Display = 'Off';
                op.Algorithm = 'Trust-Region';
                op.MaxIter = 600;
            if res == 1
                op.Lower =      P.op.Lower;
                op.StartPoint = P.op.StartPoint;
                op.Upper =      P.op.Upper;
            else
%                                 f_inter                  f_fast                   D_slow                   D_inter                  D_fast                   S0
                 op.Lower =      [a_res(x,y)*(1-P.Tol(1))  b_res(x,y)*(1-P.Tol(1))  c_res(x,y)*(1-P.Tol(2))  d_res(x,y)*(1-P.Tol(2))  e_res(x,y)*(1-P.Tol(2))  f_res(x,y)*(1-P.Tol(3))];  
                 op.StartPoint = [a_res(x,y)               b_res(x,y)               c_res(x,y)               d_res(x,y)               e_res(x,y)               f_res(x,y)];
                 op.Upper =      [a_res(x,y)*(1+P.Tol(1))  b_res(x,y)*(1+P.Tol(1))  c_res(x,y)*(1+P.Tol(2))  d_res(x,y)*(1+P.Tol(2))  e_res(x,y)*(1+P.Tol(2))  f_res(x,y)*(1+P.Tol(3))];
 
            end
            if Mask_res(x,y)
               
                if ~isnan(Data_res(x,y,1))
                   [FitResults{x,y}, gof{x,y}, output{x,y}] = TriexpFit(P.b_values, squeeze(Data_res(x,y,:))', op);
                   a(x,y) = FitResults{x,y}.a;
                   b(x,y) = FitResults{x,y}.b;
                   c(x,y) = FitResults{x,y}.c;
                   d(x,y) = FitResults{x,y}.d;
                   e(x,y) = FitResults{x,y}.e;
                   f(x,y) = FitResults{x,y}.f;
                end
              
            end
        end

    end
    
   % Interpolate parameters
   if res < size(P.Dims_steps, 1)
        a_res = imresize(a,[P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        b_res = imresize(b,[P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        c_res = imresize(c,[P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        d_res = imresize(d,[P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        e_res = imresize(e,[P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        f_res = imresize(f,[P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
   end  

end
P.time = toc(tStart);

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
    fignm_param = sprintf([P.outputFolder,'/IDEALfit_%s_steps_%s_param.fig'], file_name, num2str(size(P.Dims_steps, 1)));
    savefig(gcf, fignm_param);
    close(gcf);
end

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


filenm = sprintf([P.outputFolder '/IDEALfit%s_steps_%s.mat'], file_name, num2str(size(P.Dims_steps, 1)));
save(filenm);

end
