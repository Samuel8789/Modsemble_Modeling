function [fitresult, gof] = createFit_C(pb, LL, bool_plot)
%CREATEFIT(PB,LL)
%  Create a fit.
%
%  Data for 'untitled fit 2' fit:
%      X Input : pb
%      Y Output: LL
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 15-Apr-2020 12:42:21


%% Fit: 'untitled fit 2'.
[xData, yData] = prepareCurveData( pb, LL );

% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [-108.227429502131 3.8811786611033e-10 105.507282293749 -3.70476491091439e-05];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

if bool_plot == 1
    % Plot fit with data.
    figure( 'Name', 'untitled fit 2' );
    h = plot( fitresult, xData, yData );
    legend( h, 'LL vs. pb', 'untitled fit 2', 'Location', 'NorthEast', 'Interpreter', 'none' );
    % Label axes
    xlabel( 'pb', 'Interpreter', 'none' );
    ylabel( 'LL', 'Interpreter', 'none' );
    grid on
end

