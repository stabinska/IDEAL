function op_out = update_fitting(op, fit_res, Params, x, y)
%update_fitting(op, fit_res, Params, x, y)
% Update fitting parameters according to downsampling
%
%   op      : fitting options
%   fit_res : resampled fitting parameters
%   Params  : Parameter struct
%   x,y     : current x and y coordinate

    op_out = op;
    switch Params.Model
        case {"biexp","Biexp"}
            op_out.Lower = zeros(1,4);
            op_out.StartPoint = zeros(1,4);
            op_out.Upper = zeros(1,4);
    
            % f_fast
            op_out.Lower(1) = fit_res.f_fast(x,y)*(1-Params.Tol(1));
            op_out.StartPoint(1) = fit_res.f_fast(x,y);
            op_out.Upper(1) = fit_res.f_fast(x,y)*(1+Params.Tol(1));
    
            % D_slow
            op_out.Lower(2) = fit_res.D_slow(x,y)*(1-Params.Tol(2));
            op_out.StartPoint(2) = fit_res.D_slow(x,y);
            op_out.Upper(2) = fit_res.D_slow(x,y)*(1+Params.Tol(2));
    
            % D_fast
            op_out.Lower(3) = fit_res.D_fast(x,y)*(1-Params.Tol(2));
            op_out.StartPoint(3) = fit_res.D_fast(x,y);
            op_out.Upper(3) = fit_res.D_fast(x,y)*(1+Params.Tol(2));
    
            % S_0
            op_out.Lower(4) = fit_res.S_0(x,y)*(1-Params.Tol(3));
            op_out.StartPoint(4) = fit_res.S_0(x,y); 
            op_out.Upper(4) = fit_res.S_0(x,y)*(1+Params.Tol(3));
    
        case {"triexp","Triexp"}
            % f_inter f_fast D_slow D_inter D_fast S0
            op_out.Lower = zeros(1,6);
            op_out.StartPoint = zeros(1,6);
            op_out.Upper = zeros(1,6);
    
            % f_inter
            op_out.Lower(1) = fit_res.f_inter(x,y)*(1-Params.Tol(1));
            op_out.StartPoint(1) = fit_res.f_inter(x,y);
            op_out.Upper(1) = fit_res.f_inter(x,y)*(1+Params.Tol(1));
            
            % f_fast
            op_out.Lower(2) = fit_res.f_fast(x,y)*(1-Params.Tol(1));
            op_out.StartPoint(2) = fit_res.f_fast(x,y);
            op_out.Upper(2) = fit_res.f_fast(x,y)*(1+Params.Tol(1));
    
            % D_slow
            op_out.Lower(3) = fit_res.D_slow(x,y)*(1-Params.Tol(2));
            op_out.StartPoint(3) = fit_res.D_slow(x,y);
            op_out.Upper(3) = fit_res.D_slow(x,y)*(1+Params.Tol(2));
    
            % D_inter
            op_out.Lower(4) = fit_res.D_inter(x,y)*(1-Params.Tol(2));
            op_out.StartPoint(4) = fit_res.D_inter(x,y);
            op_out.Upper(4) = fit_res.D_inter(x,y)*(1+Params.Tol(2));
    
            % D_fast
            op_out.Lower(5) = fit_res.D_fast(x,y)*(1-Params.Tol(2));
            op_out.StartPoint(5) = fit_res.D_fast(x,y);
            op_out.Upper(5) = fit_res.D_fast(x,y)*(1+Params.Tol(2));
    
            % S_0
            op_out.Lower(6) = fit_res.S_0(x,y)*(1-Params.Tol(3));
            op_out.StartPoint(6) = fit_res.S_0(x,y); 
            op_out.Upper(6) = fit_res.S_0(x,y)*(1+Params.Tol(3));
    end
end