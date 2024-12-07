function [Mask_res, Data_res] = resample_data(Mask,Data,steps,b_values)
%resample_data(Mask 'double', Data 'double', steps 'double', b_values 'double')
% resample Matrix and Data for downsampling
%
%   Mask    : array containing masking 
%   Data    : diffusion data
%   steps   : current resampling step matrix
%   b_values: array containing coresponding b_values

    % Downsampling Image Matrix to desired size
    Mask_res = imresize(Mask, [steps(2), steps(1)],'bilinear');
    Data_res = zeros(steps(1), steps(2), numel(b_values)); 

    for i = 1:numel(b_values) 
         Data_res(:,:,i) = imresize(squeeze(Data(:, :, i)),[steps(2) steps(1)],'bilinear');               
    end
    

    for i = 1:size(b_values,2) 
        temp = Data(:, :, i);
        temp(Mask == 0) = NaN;
        dx = size(Data, 1)/steps(1);
        [r,c] = ndgrid(1:size(temp,1), 1:size(temp,2));
        [~, ibin] = histc(r(:), 0.5:dx:size(temp,1)+0.5);
        [~, jbin] = histc(c(:), 0.5:dx:size(temp,2)+0.5);
        nr = max(ibin);
        nc = max(jbin);
        idx = sub2ind([nr nc], ibin, jbin); 
        temp2 = accumarray(idx, temp(:), [nr*nc 1], @nanmean);
        Data_res(:,:,i) = reshape(temp2, nr, nc);             
    end


    Mask_res = abs(Mask_res);
    
    % Thresholding edges
    Mask_res(Mask_res < 0.025) = 0;

    % figure, imagesc(squeeze(Data_res(:,:,1)).*Mask_res);
end