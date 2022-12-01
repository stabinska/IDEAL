function [fit, FitResults, gof, output, op] = setup_fitting(Params,res_step,sz)
%% function setup_fitting(Params)
% setup variables for fitted parameters
% sz is the image/mask size

    fit = struct();
    
    % Basic ADC Parameters
    fit.S_0 = zeros(Params.Dims_steps(res_step,1), ...
                    Params.Dims_steps(res_step,2));
    fit.D_slow = zeros(Params.Dims_steps(res_step,1), ...
                    Params.Dims_steps(res_step,2));
    
    % Parameters for Bi and Tri
    if strcmp(Params.Model,'biexp') || strcmp(Params.Model,'triexp')          
        fit.f_fast = zeros(Params.Dims_steps(res_step,1), ...
                    Params.Dims_steps(res_step,2)); % f_fast
        fit.D_fast = zeros(Params.Dims_steps(res_step,1), ...
                    Params.Dims_steps(res_step,2)); % D_fast
    end
    % Parameters for Tri
    if strcmp(Params.Model,'triexp')
        fit.f_inter = zeros(Params.Dims_steps(res_step,1), ...
                            Params.Dims_steps(res_step,2)); % f_inter
        fit.D_inter = zeros(Params.Dims_steps(res_step,1), ...
                            Params.Dims_steps(res_step,2)); % D_inter   
    end  
    
    FitResults = cell(sz(1), sz(2));
    gof = FitResults;
    output = FitResults;
    
    op = Params.op;
end