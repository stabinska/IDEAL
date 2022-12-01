function fig = plot_params_figs(Params, fit, path_Data)
% plot and save figures for triexp model
    fig = figure('Visible','on');
    switch Params.Model
        case {"biexp","Biexp"}
            fit_params = ["f_slow";"f_fast";"D_slow","D_fast","S_0"];
        case {"triexp","Triexp"}
            fit_params = ["f_slow";"f_inter";"f_fast";
                            "D_slow";"D_inter";"D_fast";"S_0"];
            nx = 3; ny = 3;
    end
    for nparam = 1:length(fit_params)
        subplot(nx,ny,nparam)
        imagesc(fit.(fit_params(nparam)));
        title(fit_params(nparam));
        colormap gray;
        axis off;
        if contains(fit_params(nparam),"f_")
            caxis(gca,[0,1]);
        end
    end
    [~,file_name,~] = fileparts(path_Data);
    if ~exist(Params.outputFolder,"dir")
        mkdir(Params.outputFolder);
    end
    fig_name = Params.outputFolder+filesep+"IDEALfit_"+ ...
                    file_name+ "_steps_"+ ...
                    num2str(size(Params.Dims_steps,1))+"_param.fig";
    savefig(gcf, fig_name);
    close(gcf);
end