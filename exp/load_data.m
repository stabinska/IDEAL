function [Data,Mask,Data_masked,ROIs,Params] = load_data(path_Data,path_Mask,Params,varargin)

%% Function load_Data
% 
% load different niftis for further processing
% if the ROI path cell is not supplied the function will fill it with the
% mask
%
% path_Data         - string or char pointing to the DWI image
% path_Mask         - string or char pointing to the Mask image
% varargin{1}       - cell array containing strings or chars to ROIs


% Load Data
if isstring(path_Data) || ischar(path_Data)
    if ischar(path_Data)
        path_Data = string(path_Data);
    end
    data_nii = double(rot90(niftiread(path_Data)));
    Data = squeeze(data_nii(:,:,:,:));
else
    error("Path_Data has to be a string or a char!")
end

% Check if Data is structured as expected
if size(Data,4) ~= length(Params.b_values)
    if size(Data,3) == length(Params.b_values)
        warning("Data dimensions seem to be swaped trying to fix...")
        Data = permute(Data,[1,2,4,3]);
    end
end 


% Load Mask
if isstring(path_Mask) || ischar(path_Mask)
    if ischar(path_Mask)
        path_Mask = string(path_Mask);
    end
    Mask = double(rot90(niftiread(path_Mask)));
    % weird code below - rework
    Mask_nan = Mask;
    Mask_nan(Mask==0) = NaN;
    if size(Mask_nan,3) ~= size(Data,3)
        warning("Number of slices is different for Data and Mask!");
    end
    Data_masked = Data.*Mask_nan;
else
    error("Path_Mask has to be a string or a char!")
end

% Load ROI
if nargin > 2 
    path_ROIs = varargin{1};
    ROIs = cell(length(path_ROIs),1);
    for nROI = 1:length(path_ROIs)
        ROIs{nROI} = double(rot90(niftiread(path_ROIs{nROI})));
        ROIs{nROI}(ROIs{nROI}==0) =NaN;
    end
else
    ROIs{1} = double(rot90(niftiread(path_Mask)));
    ROIs{1} = ROIs{1}(:,:,:); 
    ROIs{1}(ROIs{1}==0) =NaN; 
    sprintf('No ROI supplied. We will do statistics on the VOI mask.')
end

% check weather outputFolder is a string
if ~isstring(Params.outputFolder)
    Params.outputFolder = string(Params.outputFolder);
end

end