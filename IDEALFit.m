function [FitResults, FitQuality, Params, ROIstat] = IDEALFit(path_Data, Params, path_Mask, path_ROI)
%IDEALFit(path_Data, Params, path_Mask, path_ROI)
%
%   Orchestrator for the IDEAL fitting pipeline.
%   Handles all file I/O, the slice loop, result saving, plotting, and
%   optional ROI evaluation.  Pure computation is delegated to IDEALFitSlice.
%
%   Input:
%       path_Data  - path to 4D DWI NIfTI file  (.nii or .nii.gz)
%       Params     - struct with experiment and fit parameters
%       path_Mask  - path to binary mask NIfTI
%       path_ROI   - cell array of ROI NIfTI paths (optional)
%
%   Output:
%       FitResults  - cell array of cfit objects at final resolution
%       FitQuality  - cell array of goodness-of-fit structs (one per slice)
%       Params      - updated parameter struct (includes Params.time)
%       ROIstat     - ROI statistics cell array (empty when nargin < 4)
%
%
%   Authors: Julia Stabinska (jstabin3@jhmi.edu)
%            Helge Jörn Zöllner (hzoelln2@jhmi.edu)
%            Thomas Andreas Thiel (thomas.thiel@hhu.de)

%% Input validation
if nargin < 3
    error('IDEALFit requires at least path_Data, Params, and path_Mask.');
end

%% Load data
if nargin < 4
    [Data_raw, Mask_raw, Data_raw_masked, ROIs, ~, Params] = ...
        load_data(path_Data, path_Mask, Params);
else
    [Data_raw, Mask_raw, Data_raw_masked, ROIs, ~, Params] = ...
        load_data(path_Data, path_Mask, Params, path_ROI);
end

%% Initialise output cell arrays
FitQuality = cell(size(Data_raw, 3), 1, 1);
ROIstat    = cell(size(Data_raw, 3), 1, 1);

tStart = tic;

%% Slice loop
for slice = Params.slice

    % Extract single slice (squeeze to 3-D / 2-D)
    Data        = squeeze(Data_raw(:,:,slice,:));
    Data_masked = squeeze(Data_raw_masked(:,:,slice,:));
    Mask        = squeeze(Mask_raw(:,:,slice,:));

    % Skip slices that contain no data after masking
    if isempty(Data_masked(~isnan(Data_masked)))
        FitQuality{slice} = {};
        ROIstat{slice}    = {};
        fprintf('Slice %s does not contain any data after masking and will be skipped.\n', ...
            num2str(slice));
        continue
    end

    % --- Pure computation ---
    [FitMaps, FitResults, gof, output] = IDEALFitSlice(Data, Mask, Data_masked, Params);

    % --- Build Experiment struct ---
    Experiment.data = Data;
    Experiment.mask = Mask;
    Experiment.ROIs = ROIs;

    % --- Summarise goodness-of-fit ---
    [FitMaps, FitQuality{slice}] = summarize_results(FitMaps, gof, output, size(Mask), Params.Model);

    % --- Plot ---
    if Params.plot
        plot_params_figs(Params, FitMaps, path_Data);
    end

    % --- Save fit maps ---
    [~, file_name, ~] = fileparts(path_Data);
    fname = Params.outputFolder + filesep + "IDEALfit_" + ...
            file_name + "_" + string(Params.Model) + "_steps_" + ...
            num2str(size(Params.Dims_steps, 1)) + ...
            "_sl_" + string(slice) + "_fit.mat";
    save(fname, "FitMaps");

    % --- Save experiment struct ---
    fname = Params.outputFolder + filesep + "IDEALfit_" + ...
            file_name + "_" + string(Params.Model) + "_steps_" + ...
            num2str(size(Params.Dims_steps, 1)) + ...
            "_sl_" + string(slice) + "_data.mat";
    save(fname, "Experiment");

    % --- ROI evaluation ---
    if nargin == 4
        ROIstat{slice} = eval_ROIS(ROIs, FitMaps, Params.Model);
        fname = Params.outputFolder + filesep + "IDEALfit_" + ...
                file_name + "_" + string(Params.Model) + "_steps_" + ...
                num2str(size(Params.Dims_steps, 1)) + ...
                "_sl_" + string(slice) + "_ROIstat.mat";
        save(fname, "ROIstat");
        fprintf("ROIs evaluated!\n");
    end

    fprintf("Fitting Completed!\n");
end

Params.time = toc(tStart);
