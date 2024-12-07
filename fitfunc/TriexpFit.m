function [fitresult, gof, output] = TriexpFit(b, data, op)
%CREATEFIT(B,MEDULLA_1)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : b
%      Y Output: Medulla_1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 19-Aug-2020 15:23:59


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( b, data );

% Set up fittype and options.
ft = fittype( 'f*((1-a-b)*exp(-c*x)  + a*exp(-d*x)  + b*exp(-e*x))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%opts.Algorithm = 'Levenberg-Marquardt';
%opts.Algorithm = 'Trust-Region';
%opts.DiffMinChange = 1e-07;
opts.Display = 'Off';

opts.Algorithm = op.Algorithm;
opts.Display = op.Display;
opts.StartPoint = op.StartPoint ;

if strcmp(opts.Algorithm, 'Trust-Region')
    opts.Lower =      op.Lower;
    opts.Upper =      op.Upper;
end


% Fit model to data.
[fitresult, gof, output] = fit( xData, yData, ft, opts );

% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'data vs. b', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel b
% ylabel data
% grid on

