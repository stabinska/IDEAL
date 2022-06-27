function ROIs = load_ROIS(MaskNii,ROINii,load_rois)

if load_rois
    ROIs{1} = double(rot90(niftiread(MaskNii)));
    ROIs{1} = ROIs{1}(:,:,:); 
    ROIs{1}(ROIs{1}==0) =NaN; 
    for rois = 2 : length(ROINii)+1
        ROIs{rois} = double(rot90(niftiread(ROINii{rois-1})));
        ROIs{rois}(ROIs{rois}==0) =NaN;
    end
else
   ROIs{1} = double(rot90(niftiread(MaskNii)));
   ROIs{1} = ROIs{1}(:,:,:); 
   ROIs{1}(ROIs{1}==0) =NaN; 
   sprintf('No ROI supplied. We will do statistics on the VOI mask.')
end
end