function [Mask_res, Data_res] = resample_data(Mask,Data,steps,b_values)
% resample Matrix and Data for downsampling

    % Downsampling Image Matrix to desired size
    Mask_res = imresize(Mask, [steps(2), steps(1)],'bilinear');
    Data_res = zeros(steps(1), steps(2),16);  

    for i = 1:numel(b_values) 
         Data_res(:,:,i) = imresize(squeeze(Data(:, :, i)),[steps(2) steps(1)],'bilinear');               
    end

    Mask_res = abs(Mask_res);
    % Thresholding edges
    Mask_res(Mask_res < 0.025) = 0;
end