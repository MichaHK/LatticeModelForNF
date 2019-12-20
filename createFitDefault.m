function [NewDistances,NewPressures,bulk] = createFitDefault(VarName2, VarName1,SmoothingParam)
%createFitDefault(Pressure,Distance)
%  Data for fit:
%      X Input : VarName2
%      Y Output: VarName1
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

R_cyl=5;

OriginalPressure=VarName2;
OriginalDistance=VarName1;
VarName2=log10(VarName2);              
[xData, yData] = prepareCurveData( VarName2, VarName1 );
model='SplineParam';    % {SplineParam,SplineDef,ShapePerserve,Exponents}
switch model
    case 'SplineParam'
        % Set up fittype and options.
        ft = fittype( 'smoothingspline' );
        opts = fitoptions( 'Method', 'SmoothingSpline' );
        opts.SmoothingParam = SmoothingParam; 
        [fitresult, gof] = fit( xData, yData, ft, opts );
    case 'SplineDef'
        ft = fittype( 'smoothingspline' );
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft );
    case 'ShapePerserve'
        ft = 'pchipinterp';
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, 'Normalize', 'on' );
    case 'Exponents'
        ft = fittype( 'exp2' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
%         opts.Display = 'Off';
%         opts.StartPoint = [228.134652004577 -1.43033949823438e-05 233.016002327094 -1.18163144127127e-07];
        % Fit model to data.
        [fitresult, gof] = fit( xData, yData, ft, opts );
end
    
NewEXP_Pressures=linspace(min(xData),max(xData),600);
NewDistances=feval(fitresult,NewEXP_Pressures);
NewPressures=10.^NewEXP_Pressures;
diffX=diff(NewDistances);
diffY=diff(NewPressures');
D=NewDistances(2:end);
bulk_cyl=-0.5*D.*(diffY./diffX);    % 0.5 because the volume is pi*R^2, not constant*R.
bulk_Hex=-(0.25*sqrt(3)*D.^2-0.5*pi*R_cyl^2)*2/sqrt(3)./D.*(diffY./diffX);    % 0.5 because the volume is pi*R^2, not constant*R.
bulk=bulk_Hex;

if (0)
figure( 'Name', 'Brrrr' );

subplot(2,1,1)
semilogy(NewDistances,NewPressures);
hold on; 
scatter(OriginalDistance,OriginalPressure,'r','o')
% legend( h, 'VarName1 vs. VarName2', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel( 'VarName2' );
% ylabel( 'VarName1' );
subplot(2,1,2)
semilogy(NewDistances(1:end-1),bulk)
grid on
end

