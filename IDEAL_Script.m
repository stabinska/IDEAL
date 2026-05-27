%%  IDEAL Script

%% === Paths (edit before running) ===
path_data   = '';   % 4D DWI NIfTI  (.nii or .nii.gz)
path_mask   = '';   % binary mask NIfTI
path_output = '';   % output folder  (created automatically if absent)
path_rois   = {};   % cell of ROI NIfTIs — leave {} to skip ROI statistics

%% === Acquisition and model parameters ===

P.op.Display   = 'Off';
P.op.Algorithm = 'Trust-Region';
P.op.MaxIter   = 600;

P.Model = 'Biexp_T1corr'; % 'Biexp'          — two-compartment IVIM
                          % 'Biexp_T1corr'   — two-compartment IVIM with T1 correction
                          % 'Triexp'         — three-compartment

%   Start point, lower/upper bounds, and per-category tolerance levels.
%   Tol(n+1) is applied as ±Tol(n)*param_n around the estimate from step n.

switch P.Model
    case {"biexp","Biexp"}
%                         f_fast D_slow D_fast S_0
        P.op.Lower =      [0.01  0.0011 0.003  10];
        P.op.StartPoint = [0.1   0.0015 0.015  100];
        P.op.Upper =      [0.9   0.003  0.3    2000];

    case {"biexp_T1corr","Biexp_T1corr"}
%                         f_fast D_slow D_fast S_0  T1
        P.op.Lower =      [0.01  0.0010 0.01   10   600];
        P.op.StartPoint = [0.1   0.0015 0.015  100  1500];
        P.op.Upper =      [0.7   0.01   0.5    2000 4000];

    case {"triexp","Triexp"}
%                         finterm ffast Dslow  Dinterm Dfast  S0
        P.op.Lower =      [0.1    0.01  0.0011 0.003   0.01   10];
        P.op.StartPoint = [0.2    0.1   0.0015 0.005   0.1    210];
        P.op.Upper =      [0.7    0.7   0.003  0.01    0.5    1000];
end

%   Tolerances: one value per parameter category (fraction, D, S0 [, T1])
P.Tol = [0.2 0.1 0.5 0.2];

%   Resampling steps — final row must equal the original image matrix size
P.Dims_steps = [1   1;
                2   2;
                4   4;
                8   8;
                16  16;
                32  32;
                64  64;
                96  96;
                128 128;
                152 152;
                176 176];

%   Acquisition parameters
P.b_values = [0, 5, 10, 15, 20, 25, 30, 50, 70, 100, 150, 250, 350, 450, 550, 750];
P.slice    = 1;
P.TR       = 1900;   % ms
P.TE       = 49;     % ms
P.TM       = 9.8;    % ms  (mixing time for DW-STEAM; set 0 if not applicable)

P.b_threshold  = 1;      % index of the first b-value to include
P.mean         = 1;      % 1 = average masked voxels in the first resampling step
P.plot         = 1;      % 1 = generate output figures
P.outputFolder = path_output;

%% === Run ===
if isempty(path_rois)
    [FitResults, FitQuality, P, ROIstat] = IDEALFit(path_data, P, path_mask);
else
    [FitResults, FitQuality, P, ROIstat] = IDEALFit(path_data, P, path_mask, path_rois);
end
