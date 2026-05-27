function [FitMaps, FitResults, gof, output] = IDEALFitSlice(Data, Mask, Data_masked, Params)
%IDEALFitSlice(Data, Mask, Data_masked, Params)
%
%   Pure-computation core of the IDEAL fitting pipeline.
%   Runs the iterative resampling loop and voxel-wise fitting for a
%   single slice.  Contains no file I/O, no plotting, and no timing.
%   The calling orchestrator (IDEALFit) is responsible for those concerns.
%
%   Input:
%       Data        - 3D signal array for one slice  (X, Y, n_bvalues)
%       Mask        - 2D mask for the slice           (X, Y)
%       Data_masked - 3D masked signal array          (X, Y, n_bvalues)
%       Params      - struct with experiment and fit parameters
%
%   Output:
%       FitMaps     - struct with final-resolution parameter maps
%       FitResults  - cell array of cfit objects at final resolution
%       gof         - cell array of goodness-of-fit structs
%       output      - cell array of fit output structs


%% Resampling loop
for res_step = 1:size(Params.Dims_steps, 1)
    fprintf('Downsampling step no: %s\n', num2str(res_step));

    if res_step == 1 && Params.mean
        % Average all masked voxels in the first step
        Mask_res = ones(1, 1);
        Data_res = zeros(1, 1, numel(Params.b_values));
        for i = 1:numel(Params.b_values)
            Data_res(:,:,i) = nanmean(Data_masked(:,:,i), 'all');
        end

    elseif res_step < size(Params.Dims_steps, 1)
        % Downsample to intermediate resolution
        [Mask_res, Data_res] = resample_data( ...
            Mask, Data, ...
            [Params.Dims_steps(res_step,1), Params.Dims_steps(res_step,2)], ...
            Params.b_values);

    else
        % Final step: use full-resolution masked data
        Data_res = Data_masked;
        Mask_res = Mask;
    end

    % Prepare fitting parameters (builds op.ft once for the whole step)
    [fit, FitResults, gof, output, op] = setup_fitting(Params, res_step, size(Mask_res));

    %% Voxel-wise fitting
    for x = 1:size(Mask_res, 1)
        for y = 1:size(Mask_res, 2)

            % Update start points from the previous resampling step
            if res_step ~= 1
                op = update_fitting(op, fit_res, Params, x, y);
            end

            if Mask_res(x,y)
                if ~isnan(Data_res(x,y,1))
                    % Fit only masked voxels with valid signal
                    switch Params.Model
                        case {"biexp","Biexp"}
                            [FitResults{x,y}, gof{x,y}, output{x,y}] = ...
                                BiexpFit(Params.b_values, ...
                                         squeeze(Data_res(x,y,:))', op);

                        case {"biexp_T1corr","Biexp_T1corr"}
                            [FitResults{x,y}, gof{x,y}, output{x,y}] = ...
                                BiexpFit_T1corr(Params.b_values, ...
                                                squeeze(Data_res(x,y,:))', op, ...
                                                Params.TM, Params.TR, Params.TE);

                        case {"triexp","Triexp"}
                            [FitResults{x,y}, gof{x,y}, output{x,y}] = ...
                                TriexpFit(Params.b_values, ...
                                          squeeze(Data_res(x,y,:))', op);
                    end

                    fit = extract_fit(fit, FitResults, Params.Model, x, y);
                end
            end
        end
    end

    % Interpolate parameter maps for use as start points in the next step
    if res_step < size(Params.Dims_steps, 1)
        fit_res = interpolate_fit(fit, Params.Dims_steps, res_step, Params.Model);
    end
end

%% Return final-resolution parameter maps
FitMaps = fit;
