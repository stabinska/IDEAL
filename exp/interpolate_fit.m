function fit_res = interpolate_fit(fit,Dims_steps,res_step,Model)
    % interpolate parametermaps to new size

    switch Model
        case {"biexp","Biexp"}
            params = ["f_fast";"D_slow";"D_fast";"S_0"];
        case {"triexp","Triexp"}
            params = ["f_inter";"f_fast";"D_slow";"D_inter";"D_fast";"S_0"];
    end

    for nparam = 1:length(params)
        fit_res.(params(nparam)) = imresize(fit.(params(nparam)),...
            [Dims_steps(res_step+1,2), Dims_steps(res_step+1,1)], 'bilinear');
    end
end