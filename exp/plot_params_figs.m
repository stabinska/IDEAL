function fig = plot_params_figs(Params, fit, path_Data)
% plot and save figures for triexp model


    switch Params.Model
        % check model
        case {"biexp","Biexp"}
            fit_params = ["f_slow";"f_fast";"D_slow";"D_fast";"S_0"];
            titles = ["f_{slow}";"f_{fast}";"D_{slow}";"D_{fast}";"S_{0}"];
            nx = 3; ny = 2;
        case {"triexp","Triexp"}
            fit_params = ["f_slow";"f_inter";"f_fast";
                            "D_slow";"D_inter";"D_fast";"S_0"];
            titles = ["f_{slow}";"f_{inter}";"f_{fast}";
                            "D_{slow}";"D_{inter}";"D_{fast}";"S_{0}"]; 
            nx = 3; ny = 3;
    end
    % detect number of slices
    nslices = size(fit.S_0,3);
    for nslice = 1:nslices

        fig = figure('Visible','on');    
    
        % start ploting
        for nparam = 1:length(fit_params)
            subplot(nx,ny,nparam)
            imagesc(fit.(fit_params(nparam))(:,:,nslice));
            title(titles(nparam));
            colormap gray;
            axis off;
            if contains(fit_params(nparam),"f_")
                caxis(gca,[0,1]);
            end
        end
        % check if output-folder exists
        [~,file_name,~] = fileparts(path_Data);
        if ~exist(Params.outputFolder,"dir")
            mkdir(Params.outputFolder);
        end
        % save figures
        fig_name = Params.outputFolder + filesep + "IDEALfit_" + ...
                        file_name + "_" + string(Params.Model) + "_steps_" + ...
                        num2str(size(Params.Dims_steps,1)) + ...
                        "_sl_" + string(nslice) + "_param.fig";
        savefig(gcf, fig_name);
        close(gcf);
    end 
end