function [FitResults,FitQuality,P,ROIstat] = IDEALfitIVIM_tt(DataNii,P,MaskNii,ROINii)

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
    % res : current downsampling step    
    fprintf('Downsampling step no: %s\n', num2str(res));
    
    % Basic ADC Parameters
    S_0 = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % S_0
    D_slow = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % D_slow
    
    % Parameters for Bi and Tri
    if strcmp(P.op.Model,'Biexp') || strcmp(P.op.Model,'Triexp')          
        f_fast = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % f_fast
        D_fast = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % D_fast
    end
    
    % Parameters for Tri
    if strcmp(P.op.Model,'Triexp')
        f_inter = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % f_inter
        D_inter = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % D_inter   
    end        
    
    % Downsample Matrix and Data
    if res < size(P.Dims_steps, 1)
        % Downsampling Image Matrix to desired size
        Mask_res = imresize(Mask, [P.Dims_steps(res,2) P.Dims_steps(res,1)],'bilinear');
        Data_res = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2),16);  
        
        for i = 1:numel(P.b_values) 
             Data_res(:,:,i) = imresize(squeeze(Data(:, :, i)),[P.Dims_steps(res,2) P.Dims_steps(res,1)],'bilinear');               
        end
        
        Mask_res = abs(Mask_res);
        % Thresholding edges
        Mask_res(Mask_res < 0.025) = 0;
    else
        % Using original resolution
        % What if imagesize ~= 176:176???
        Data_res = Data_masked;
        Mask_res = Mask;
    end
    
    FitResults = cell(size(Mask_res, 1), size(Mask_res, 2));
    gof = FitResults;
    output = FitResults;    
    
    % Load starting values
    if res == 1
        op.Lower =      P.op.Lower;
        op.StartPoint = P.op.StartPoint;
        op.Upper =      P.op.Upper;
    end 
    
    % Voxelvise fitting iterators
    for x = 1 : size(Mask_res, 1)
        for y = 1: size(Mask_res, 2)           
              
            if res > 1
                % Setup Parameter for corresponding res
                switch P.op.Model 
                    case 'ADC'
                       [op.Lower,op.StartPoint,op.Upper] = ...
                            set_fitting_boundries(P.Tol,P.op.Model,S_0_res,...
                            D_slow_res); 
                    case 'Biexp'
                        [op.Lower,op.StartPoint,op.Upper] = ...
                            set_fitting_boundries(P.Tol,P.op.Model,S_0_res,...
                            D_slow_res,f_fast_res,D_fast_res);
                    case 'Triexp'
                        [op.Lower,op.StartPoint,op.Upper] = ...
                            set_fitting_boundries(P.Tol,P.op.Model,S_0_res,...
                            D_slow_res,f_fast_res,D_fast_res,...
                            f_inter_res,D_inter_res);
                end
            end
            if Mask_res(x,y)
                if ~isnan(Data_res(x,y,1))
                    switch P.op.Model
                        case 'ADC'
                        case 'Biexp'
                            [FitResults{x,y}, gof{x,y}, output{x,y}] = BiexpFit(...
                                   P.b_values, squeeze(Data_res(x,y,:))', op);
                            f_fast(x,y)  = FitResults{x,y}.a; % f_fast
                            D_slow(x,y)  = FitResults{x,y}.b; % D_slow
                            D_fast(x,y)  = FitResults{x,y}.c; % D_fast
                            S_0(x,y)     = FitResults{x,y}.d; % S_0
                               
                        case 'Triexp'
                            [FitResults{x,y}, gof{x,y}, output{x,y}] = TriexpFit(...
                                   P.b_values, squeeze(Data_res(x,y,:))', op);
                            f_inter(x,y) = FitResults{x,y}.a; % f_inter
                            f_fast(x,y)  = FitResults{x,y}.b; % f_fast
                            D_slow(x,y)  = FitResults{x,y}.c; % D_slow
                            D_inter(x,y) = FitResults{x,y}.d; % D_inter
                            D_fast(x,y)  = FitResults{x,y}.e; % D_fast
                            S_0(x,y)     = FitResults{x,y}.f; % S_0
                    end                   
                end

            end
        end

    end
    
   % Interpolate parameters
   if res < size(P.Dims_steps, 1)       
        S_0_res = imresize(S_0,...
            [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');       
        D_slow_res = imresize(D_slow,...
            [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        
        if strcmp(P.op.Model,'Biexp') || strcmp(P.op.Model,'Triexp')  
            f_fast_res = imresize(f_fast,...
                [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');        
            D_fast_res = imresize(D_fast,...
                [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        end
        
        if strcmp(P.op.Model,'Triexp')
            f_inter_res = imresize(f_inter,...
                [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
            D_inter_res = imresize(D_inter,...
                [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
        end
   end  

end
P.time = toc(tStart);


end


%% NESTED FUNCTIONS

function [Lower,Upper,StartPoint] = set_fitting_boundries(Tol,Model,S_0,D1,varargin)
%%
%
%
% varargin : {f_slow, D_slow, f_inter, D_inter}

    switch Model 
        case 'ADC'
            Lower = zeros(1,2);
            StartPoint = zeros(1,2);
            Upper = zeros(1,2);
            
            % D1 = D_slow
            Lower(1) = D1*(1-Tol(2));
            StartPoint(1) = D1;
            Upper(1) = D1*(1+Tol(2));
        case 'Bi'
            if nargin < 2
                error('Missing f_slow and/or D_slow'); end
            Lower = zeros(1,4);
            StartPoint = zeros(1,4);
            Upper = zeros(1,4);
            
            % slow
            f2 = varargin{1};
            D2 = varargin{2};            
            
            % f_fast
            Lower(1) = f2*(1-Tol(1));
            StartPoint(1) = f2;
            Upper(1) = f2*(1+Tol(1));
            
            % D_slow
            Lower(2) = D1*(1-Tol(2));
            StartPoint(2) = D1;
            Upper(2) = D1*(1+Tol(2));
            
            % D_fast
            Lower(3) = D2*(1-Tol(2));
            StartPoint(3) = D2;
            Upper(3) = D2*(1+Tol(2));
            
        case 'Tri'                        
            Lower = zeros(1,6);
            StartPoint = zeros(1,6);
            Upper = zeros(1,6);
            
            if nargin < 4
                error('Missing f_inter and/or D_inter'); end
            if nargin < 2
                error('Missing f_slow and/or D_slow'); end
            
            % slow
            f2 = varargin{1};
            D2 = varargin{2};
            % inter
            f3 = varargin{3};
            D3 = varargin{4};
            
            % f_fast
            Lower(1) = f2*(1-Tol(1));
            StartPoint(1) = f2;
            Upper(1) = f2*(1+Tol(1));
            
            % f_inter
            Lower(2) = f3*(1-Tol(1));
            StartPoint(2) = f3;
            Upper(2) = f3*(1+Tol(1));
            
            % D_slow
            Lower(3) = D1*(1-Tol(2));
            StartPoint(3) = D1;
            Upper(3) = D1*(1+Tol(2));
            
            % D_inter
            Lower(4) = D3*(1-Tol(2));
            StartPoint(4) = D3;
            Upper(4) = D3_res*(1+Tol(2));
            
            % D_fast
            Lower(5) = D2*(1-Tol(21));
            StartPoint(5) = D2;
            Upper(5) = D2*(1+Tol(2));
    end
    
    % add S_0
    Lower(end)     = S_0*(1-Tol(3));
    StartPoint(end)= f_res;
    Upper(end)     = S_0*(1+Tol(3));
    
end
function IDEALevalTri(Mask,FitResults,gof,output,DataNii,P,ROIs,ROInii)

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