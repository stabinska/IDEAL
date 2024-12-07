%%  IDEAL Script
%
%   Define the parameters for IDEAL fitting
%   Struct with settings
%   NLLS fit settings
    P.op.Display = 'Off';
    P.op.Algorithm = 'Trust-Region';
    P.op.MaxIter = 600;
    P.Model = 'Biexp_T1corr'; % 'Biexp' -  two-compartment IVIM fitting
                              % 'Biexp_T1corr' - two-compartment IVIM-T1 fitting
                              % 'Triexp' - three-compartment fitting

%   Start parameter, upper and lower points, and tolerance level during the
%   IDEAL resampling. The tolerance level for step n + 1 is given in percent of the
%   paramter estimated during the resampling step n, e.g.,
%   finterm_step_n +- 0.2*finterm_step_n is the startpoint and lower/upper
%   limit for step n + 1.

switch P.Model
    case {"biexp","Biexp"}        
%                         f_fast D_slow D_fast S_0         
        P.op.Lower =      [0.01  0.0011 0.003  10];  
        P.op.StartPoint = [0.1   0.0015 0.015  100];
        P.op.Upper =      [0.9   0.003  0.3    2000]; 
     
    case {"biexp_T1corr","Biexp_T1corr"}   % rev2 final
%                         f_fast D_slow D_fast S_0  T1
        P.op.Lower =      [0.01  0.0010 0.01   10   600];  
        P.op.StartPoint = [0.1   0.0015 0.015  100  1500];
        P.op.Upper =      [0.7   0.01   0.5    2000 4000];      

    case {"biexp_T1corrFixed","Biexp_T1corrFixed"}   
%                         f_fast D_slow D_fast S_0
        P.op.Lower =      [0.01  0.0011 0.003  10];  
        P.op.StartPoint = [0.1   0.0015 0.015  100];
        P.op.Upper =      [0.9   0.003  0.3    2000];
        
    case {"triexp","Triexp"}        
%                       finterm ffast Dslow Dinterm Dfast S0  
        P.op.Lower =      [0.1  0.01  0.0011 0.003  0.01  10];  
        P.op.StartPoint = [0.2  0.1   0.0015 0.005  0.1  210];
        P.op.Upper =      [0.7  0.7   0.003  0.01   0.5  1000]; 

end

  P.Tol = [0.2 0.1 0.5 0.2]; % Tolerances for f, Dslow, Dfast, S0 (and T1)
  
  
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

%  Nominal b-vlaues set during the acquisition (use effective b-values for DW-STEAM data)  
   P.b_values = [0,5,10,15,20,25,30,50,70,100,150,250,350,450,550,750];
   P.slice = 1;
   P.TR = 1900;
   P.TE = 49;
   P.b_threshold = 1; % index of the first considered b-value
   
%   Output folder for the plots (will be created automatically)
   P.outputFolder = '/Users/julia/Documents/IDEALgit/test_tri';
   P.plot = 1; % Plot flag
   P.mean = 1; % First step mean flag, if P.mean = 1 average the masked data in the first step

% Call IDEAL fit
% Include path to 3D DWI nifti file, parameter struct, nifit mask, and a
% cell array of nifit ROIs for the statistics (optional)

P.TM = 9.8; % Mixing time (if DW-STEAM was used)
Data = '/Users/julia/Documents/DeltaNiereIDEAL/NIfTIs/DeltaNiere00_Delta26_sl3.nii'; %Path to data
MaskNii = '/Users/julia/Documents/DeltaNiereIDEAL/Masks/original/DeltaNiere00_sl3_mask.nii.gz'; %Path to mask
ROIsNiiCell = {'/Users/julia/Documents/DeltaNiereIDEAL/Masks/original/DeltaNiere00_sl3_MR.nii.gz', '/Users/julia/Documents/DeltaNiereIDEAL/Masks/original/DeltaNiere00_sl3_CR.nii.gz'}; %Cell with paths to ROIs
[FitResults,FitQuality,P,ROIstat] = IDEALfitIVIM(Data,P,MaskNii,ROIsNiiCell);    


