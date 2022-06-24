%%  IDEAL Script
%
%   Define the parameters for IDEAL fitting
%
%
%%   Struct with settings
%   NLLS fit settings
    P.op.Display = 'Off';
    P.op.Algorithm = 'Trust-Region';
    P.op.MaxIter = 600;
    P.op.Model = 'Triexp';

%   Start parameter, upper and lower points, and tolerance level during the
%   IDEAL resampling. The tolerance level for step n + 1 is given in percent of the
%   paramter estimated during the resampling step n, e.g.,
%   finterm_step_n +- 0.2*finterm_step_n is the startpoint and lower/upper
%   limit for step n + 1.
%                   finterm ffast Dslow Dinterm Dfast S0  Triexp
    P.op.Lower =    [0.1  0.01  0.0011 0.003   0.01  10];  
    P.op.StartPoint = [0.2  0.1  0.0015  0.005  0.1  210];
    P.op.Upper =      [0.7 0.7 0.003  0.01   0.5  1000]; 
    P.Tol = [0.2 0.1 0.5];
%   Size of the image matrix for each resampling step. Should end with the
%   original image matrix size
    P.Dims_steps = [1 1;
                    2 2;
                    4 4;
                    8 8;
                    16 16;
                    32 32;
                    64 64;
                    96 96;
                    128 128;
                    152 152;
                    176 176];
%   b-vlaues set during the acquisition
    P.b_values = [0,10,20,30,50,70,100,150,200,250,300,350,450,550,650,750];
    P.slice = 1;
%   Output folder for the plots (will be created automatically)
    P.outputFolder = 'E:\home\Thomas\Sciebo\Projekte\Kidney_IDEAL\test_data_julia\output';
    P.plot = 1; %Plot flag

%% Call IDEAL fit
% Include path to 3D DWI nifti file, parameter struct, nifit mask, and a
% cell array of nifit ROIs for the statistics (optional)

Data = 'E:\home\Thomas\Sciebo\Projekte\Kidney_IDEAL\test_data_julia\20200804_164104DTIPGSETE71Delta20delta62s006a1001.nii'; %Path to data
MaskNii = 'E:\home\Thomas\Sciebo\Projekte\Kidney_IDEAL\test_data_julia\Delta_Niere06_sl1_niererechts.nii'; %Path to mask
ROIsNiiCell = {'E:\home\Thomas\Sciebo\Projekte\Kidney_IDEAL\test_data_julia\Delta_Niere06_sl1_cortexrechts.nii'}; %Cell with paths to ROIs

[FitResults,FitQuality,P,ROIstat] = IDEALfitIVIM(Data,P,MaskNii,ROIsNiiCell);