function [FitResults,FitQuality,Params,ROIstat] = IDEALfitIVIM_exp(path_Data,Params,path_Mask,path_ROI)
%IDEALfitIVIM_exp(path_Data 'string', Params 'struct', path_Mask 'string', path_ROI 'cell')
% 
%   Fit a bi-exponential or tri-exponential model to DWI data using the
%   Iterative Downsampling Adaptive Least-squares (IDEAL) approach
%  
%
%   Input:    path_DataNii                         - 4D-DWI data (*.nii.gz)
%               (x_pos, y_pos, z_pos, b-value)
%             Params                               - structure with experiment and fit parameters
%             path_MaskNii                         - mask (e.g. both kidneys)  (*.nii.gz)
%             path_ROINii                          - regions of interest (e.g. cortex, medulla, pelvis) (*.nii.gz) as cell 
%             
%
%   Output:   FitResults          - structure with fit results
%             FitQuality          - structure with goodnes-of-fit measures (SSE, R2, adjR2, stdError, CVs)
%             Params              - updated structure with experiment and fit parameters
%
%
%   Authors: Julia Stabinska (jstabin3@jhmi.edu)
%            Helge Jörn Zöllner (hzoelln2@jhmi.edu)
%            Thomas Andreas Thiel (thomas.thiel@hhu.de)
%

%% Parse input arguments and load data,mask and ROIs if supplied
    if nargin < 3
        error('You need to supply at least a 3D-DWI data matrix as *.nii.gz, a struct with settings (see IDEAL script), and a mask volume of the volume of interest as *.nii.gz');
    end
    
    if nargin < 4
        % if no ROIs are deployed
        [Data_raw, Mask_raw, Data_raw_masked,ROIs,Params] = ...
            load_data(path_Data,path_Mask,Params);
    elseif nargin == 4
        % else load ROIs to cell array
        [Data_raw, Mask_raw, Data_raw_masked,ROIs,Params] = ...
            load_data(path_Data,path_Mask,Params,path_ROI);
    end
    
    
    %% Perform IDEAL fitting
    tStart = tic;
%     for slice = 1:size(Data_raw,3)
    for slice = Params.slice
    
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
            
        % Loop Resampling steps
        for res_step = 1:size(Params.Dims_steps,1) 
            fprintf('Downsampling step no: %s of slice %s\n', num2str(res_step), num2str(slice));
            % Resample Data if not in finale step
            if res_step < size(Params.Dims_steps, 1)
                [Mask_res, Data_res] = ...
                    resample_data(Mask,Data, ...
                    [Params.Dims_steps(res_step,1), Params.Dims_steps(res_step,2)],...
                    Params.b_values);
            else
                Data_res = Data_masked;
                Mask_res = Mask;
            end            
            
            % Prepare fitting parameters
            [fit, FitResults, gof, output, op] = setup_fitting(Params,res_step,size(Mask_res));

            % Start Voxelvise-Fitting
            for x = 1:size(Mask_res,1)
                for y = 1:size(Mask_res,2)
                    % if matrix size is not 1x1, update parameters
                    if res_step ~= 1
                        op = update_fitting(op, fit_res, Params, x, y);
                    end
                    if Mask_res(x,y) 
                        if ~isnan(Data_res(x,y,1))
                            % only run if voxel is part of mask and data 
                            % is present else skip
                            switch Params.Model
                                case {"biexp","Biexp"}
                                    [FitResults{x,y}, gof{x,y}, ...
                                        output{x,y}] = ...
                                        BiexpFit(Params.b_values, ...
                                        squeeze(Data_res(x,y,:))', op);
                                case {"triexp","Triexp"}
                                    [FitResults{x,y}, gof{x,y}, ...
                                        output{x,y}] = ...
                                        TriexpFit(Params.b_values, ...
                                        squeeze(Data_res(x,y,:))', op);
                            end                            
                            fit = extract_fit(fit,FitResults, ...
                                            Params.Model,x,y);
                       
                        end
                    end
                end
            end 
            % interpolate parameters
            if res_step < size(Params.Dims_steps,1)
                fit_res = interpolate_fit(fit,Params.Dims_steps, ...
                                            res_step,Params.Model);
            end
        end
        fprintf("Fitting Completed!\nStarting Plotting...\n");
        Params.time = toc(tStart);
        
        % Extract Quality Parameters for final fit 
        [fit, FitQuality] = summarize_results(fit, gof, output, size(Mask_res), Params.Model);
        
        % Plot Figures and Save 
        [~] = plot_params_figs(Params, fit, path_Data);
        fprintf("Plotting Completted!\n");
    
        % Save fit struct
        [~,file_name,~] = fileparts(path_Data);
        fname = Params.outputFolder + filesep + "IDEALfit_" + ...
                file_name + "_" + string(Params.Model) + "_steps_" + ...
                "_sl_" + slice + ...
                num2str(size(Params.Dims_steps,1)) + "_fit.mat";
        save(fname,"fit");
    
        % evaluate ROIs
        if nargin == 4
            ROIstat = eval_ROIS(ROIs,fit,Params.Model);    
            fname = Params.outputFolder + filesep + "IDEALfit_" + ...
                    file_name + "_" + string(Params.Model) + "_steps_" + ...
                    "_sl_" + slice + ...
                    num2str(size(Params.Dims_steps,1)) + "_ROIstat.mat";
            save(fname,"ROIstat");
            fprintf("ROIS evaluated!\n");
        end
    end
end