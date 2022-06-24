function [Data, Mask, Data_masked] = load_files(DataNii,MaskNii,slice)
%% Load data and mask from Nifti


data_nifti = double(rot90(niftiread(DataNii)));
Data = squeeze(data_nifti(:,:,slice,:));

Mask_nifti = double(rot90(niftiread(MaskNii)));
Mask = Mask_nifti(:,:,slice); 
Mask_nan = Mask;
Mask_nan(Mask == 0) = NaN;

Data_masked = Data.*Mask_nan;
end