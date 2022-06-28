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

% check if ROIs were supplied
if nargin > 3 
    load_rois = true;
else
    load_rois = false;
end

% loop over selceted slices

% Load Data
[Data_raw, Mask_raw, Data_raw_masked] = load_files(DataNii,MaskNii);
ROIs = load_ROIS(MaskNii,ROINii,load_rois); % maybe move down to eval section


%% Perform IDEAL fitting

tStart = tic;
% clearvars a b c d e f fitresults gof output
for slice = 1:(size(Data_raw,3))
    
    % Select single slice for image resizing
    Data = squeeze(Data_raw(:,:,slice,:));
    Data_masked = squeeze(Data_raw_masked(:,:,slice,:));
    Mask = squeeze(Mask_raw(:,:,slice,:));
    
    
    % prepare Output structs
    FitQuality = cell(size(Data_raw,3),1,1);
    ROIstat = cell(size(Data_raw,3),1,1);    
    
    % Check if sclice contains data after masking. If not skip slice
    if isempty(Data_masked(~isnan(Data_masked)))
        FitQuality{slice} = {};
        ROIstat{slice} = {};
        fprintf('Slice %s does not contain any data after masking and will be skipped.\n',...
            num2str(slice));
        continue
    end
    
    for res = 1 : size(P.Dims_steps, 1)
        % res : current resampling step    
        fprintf('Downsampling step no: %s of slice %s\n', num2str(res), num2str(slice));

        % Basic ADC Parameters
        S_0 = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % S_0
        D_slow = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % D_slow

        % Parameters for Bi and Tri
        if strcmp(P.Model,'Biexp') || strcmp(P.Model,'Triexp')          
            f_fast = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % f_fast
            D_fast = zeros(P.Dims_steps(res,1), P.Dims_steps(res,2)); % D_fast
        end

        % Parameters for Tri
        if strcmp(P.Model,'Triexp')
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
            op = P.op;
        end 

        % Voxelvise fitting iterators
        for x = 1 : size(Mask_res, 1)
            for y = 1: size(Mask_res, 2)
                if res > 1
                    % Updating Start Values and Boundries for next run
                    switch P.Model 
                        case 'ADC'
                           [op.Lower,op.StartPoint,op.Upper] = ...
                                set_fitting_boundries(P.Tol,P.Model,...
                                S_0_res(x,y),D_slow_res(x,y)); 
                        case 'Biexp'
                            [op.Lower,op.StartPoint,op.Upper] = ...
                                set_fitting_boundries(P.Tol,P.Model,...
                                S_0_res(x,y),D_slow_res(x,y),...
                                f_fast_res(x,y),D_fast_res(x,y));
                        case 'Triexp'
                            [op.Lower,op.StartPoint,op.Upper] = ...
                                set_fitting_boundries(P.Tol,P.Model,...
                                S_0_res(x,y),D_slow_res(x,y),...
                                f_fast_res(x,y),D_fast_res(x,y),...
                                f_inter_res(x,y),D_inter_res(x,y));
                    end
                end
                if Mask_res(x,y)
                    if ~isnan(Data_res(x,y,1))
                        switch P.Model
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

            if strcmp(P.Model,'Biexp') || strcmp(P.Model,'Triexp')  
                f_fast_res = imresize(f_fast,...
                    [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');        
                D_fast_res = imresize(D_fast,...
                    [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
            end

            if strcmp(P.Model,'Triexp')
                f_inter_res = imresize(f_inter,...
                    [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
                D_inter_res = imresize(D_inter,...
                    [P.Dims_steps(res+1,2) P.Dims_steps(res+1,1)], 'bilinear');
            end
       end  

    end
    P.time = toc(tStart);
    
    
    % prepare ROI
    ROI{1} = squeeze(ROIs{1}(:,:,slice));
    switch p.Model
        case 'ADC'
        case 'Biexp'
            [FitQuality{slice},ROIstat{slice}] = IDEALevalBiexp(...
                FitResults,gof,output,DataNii,P,ROIs,ROINii,Data,MaskNii,slice);
        case 'Triexp'
            [FitQuality{slice},ROIstat{slice}] = IDEALevalTriexp(Mask,...
                FitResults,gof,output,DataNii,P,ROI,ROINii,Data,MaskNii,slice);
    end
end
end


%% NESTED FUNCTIONS

function [Lower,Upper,StartPoint] = set_fitting_boundries(Tol,Model,S_0,D1,varargin)
%% [Lower,Upper,StartPoint] = set_fitting_boundries(Tol,Model,S_0,D1,varargin)
% recalculate fitting starting points and boundries from previous iteration
% depending on Model choosen. 
% Basic ADC parameters have to be provided
% 
% varargin : {f_slow, D_slow, f_inter, D_inter}
%            {f_2,    D_2,    f3,      D3)

    switch Model 
        case 'ADC'
            Lower = zeros(1,2);
            StartPoint = zeros(1,2);
            Upper = zeros(1,2);
            
            % D1 = D_slow
            Lower(1) = D1*(1-Tol(2));
            StartPoint(1) = D1;
            Upper(1) = D1*(1+Tol(2));
        case 'Biexp'
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
            
        case 'Triexp'                        
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
            Upper(4) = D3*(1+Tol(2));
            
            % D_fast
            Lower(5) = D2*(1-Tol(2));
            StartPoint(5) = D2;
            Upper(5) = D2*(1+Tol(2));
    end
    
    % add S_0
    Lower(end)     = S_0*(1-Tol(3));
    StartPoint(end)= S_0;
    Upper(end)     = S_0*(1+Tol(3));
    
end