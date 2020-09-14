function [TD, maximum, transient_norm_square] =...
    transient_duration_step_changes_additional_output(...
    y,Fs,ED,Td,DD,TL,SW_Plot,nD)


%--------------------------------------------------------------------------
%                   Benchmark on Adaptive Regulation:
%
%   Rejection of unknown/time-varying multiple narrow band disturbances
%
%--------------------------------------------------------------------------
%
%   function TD = transient_duration_step_changes(y,Fs,ED,Td,DD,TL)
%
%   This function computes the transient duration for narrow band disturbance
%   rejection (Application to a Benchmark problem) for frequency step
%   changes protocol.
%
%   The function computes the transient duration for each step change in
%   frequency during a sequence of step changes. The user has to specify
%   what step change in the sequence has to be analysed (first, second, etc)
%   by specifying the varaible nD.
%
%   The transient duration is defined as the time necessary to the output
%   to reach the value +/-2.17sigma , where sigma is the standard deviation
%   for the residual output (once the algo has converged) computed during
%   the 3 last seconds of the experiment.
%
%   sigma (x) = sqrt ( var(x))
%   where "x" is the vector of data corresponding to the 3 last seconds of
%   the experiment.
%
%   In addition the output after has reached for
%   the first time 2.17 sigma should remain within +/- 4 sigma.
%   The test has been developped assuming that the residual output is a
%   gaussian process.(ie 97% of values are between +/-2.17 sigma, and
%   99.99% of values are between +/- 4 sigma).
%
%   However it is tolerated that 0.1% of the measurements be outside +/- 4
%   sigma. The percentage is denoted Dp. (to take in account that the
%   residual error is not a pure gaussian variable).
%
%
%  TD = transient_duration_step_changes(y,Fs,ED,Td,T1,DD,TL)
%
%  TD : transient duration (in seconds)
%
%  y : is vector of data
%
%  Fs : is the sampling frequency ( Fs = 800 Hz by default)
%
%  Ed : is the exepriment duration ( ED = 32 sec by default)
%       The experiment duration should be at least equal to 10 sec.
%
%  Td : Start disturbance time of the first disturbance (Td = 5 sec by default)
%
%  nD : Considered disturbance
%  nD can be equal to 1 for the first step change in frequencies, 2 for the
%  second ... and 5 for the fifth and last one.
%
%  DD : Duration of the application of the disturbance (DD = 3 sec by default)
%
%  TL : Last disturbance duration (TL = 15 sec by default)
%
%  threshold 1 = 2.17 * standard deviation (after convergence)
%  threshold 2 = 4 * standard deviation (after convergence)
%
%  Important Remark :
%  If you use the default values : Fs=800; ED=32; Td=5; DD=3; TL=15; you can
%  run the function with the syntax :
%
%  TD = transient_duration_step_changes(y)
%
%   Written by M. Alma and I.D. Landau (GIPSA-LAB)
%   Version 1, September 16, 2010
%
%   Updated by Xu Chen (xuchen@cal.berkeley.edu)
%   Version 2, 2010-09-26
%              Added residual norm calculation
%              Added a switch to disable the plots

if nargin<2,
    Fs = 800; ED = 32; Td = 5;nD =1; DD = 3; TL = 15; SW_Plot = 'PlotOff';
end
if nargin<3,
    ED = 32; Td = 5;nD =1; DD = 3; TL = 15; SW_Plot = 'PlotOff';
end
if nargin<4,
    Td = 5;nD =1; DD = 3; TL = 15; SW_Plot = 'PlotOff';
end
if nargin<5,
    nD = 1; DD = 3; TL = 15; SW_Plot = 'PlotOff';
end
if nargin<6,
    TL = 15; SW_Plot = 'PlotOff'; nD = 1;
end
if nargin<7,
    SW_Plot = 'PlotOff'; nD = 1;
end
if nargin<8,
    nD = 1;
end
%  y,Fs,ED,Td,DD,TL,SW_Plot,nD
if strcmp(SW_Plot,'PlotOn');
    SW_PLOT = 1;
else
    SW_PLOT = 0;
end

% fprintf (...
%     'Choose the considered disturbance (1- First 2- Second 3- Third 4- Fourth 5- Fifth )')
% fprintf ('\n')
% nD=input('nD=')

if Td>=ED, error('Disturbance must be injected before the end of the experiment'),end

T1 = Td + (nD-1)*DD;

y2=y(T1*Fs+1:(T1+DD)*Fs);
y2=y2-mean(y2);

y=y-mean(y);

perc = 99.9;

yres = y((ED-3)*Fs+1:ED*Fs);
sigma = sqrt(var(yres));

if SW_PLOT
    sigma
end

threshold1 = 2.17*sigma;
threshold2 = 4*sigma;
l=0;
for n = 1:length(y2)
    g=0;
    if abs(y(n)) <= threshold1
        for k = n:length(y2)
            if abs(y2(k)) > threshold2
                g= g+1;
            end
        end
        if g<(length(y2)-n+1)*((100-perc)/100)
            l=n;
            break
        end
    end
end
% g;
Dp = g/(length(y2)-n+1) *100;
maximum = max(abs(y2));
TD=(l+1)/Fs;

% T1 = Td + (nD-1)*DD;

% y2=y(T1*Fs+1:(T1+DD)*Fs);
% y2=y2-mean(y2);

transient_norm_square = sum(y2(1:TD*Fs).^2);

if l==0, error('Algorithm does not converge'),end

if SW_PLOT
    fprintf ('Percentage of points above +/- 4 sigma is :')
    Dp
    fprintf ('The transient maximum value is:')
    maximum
    fprintf ('Transient duration is :')
    TD
end

t = 0:1/800:DD-1/800;

threshold_vec1=threshold1*ones(1,length(y2));

threshold_vec2=threshold2*ones(1,length(y2));

vec_ver = -max(abs(y)):0.001:max(abs(y2));
TD_ver =  (l+1)/Fs*ones(1,length(vec_ver));
if SW_PLOT
    figure
    plot(t,y2);
    hold on
    plot(t,threshold_vec1,'r','LineWidth',2)
    hold on
    plot(t,-threshold_vec1,'r','LineWidth',2)
    hold on
    plot(t,threshold_vec2,'k','LineWidth',2)
    hold on
    plot(t,-threshold_vec2,'k','LineWidth',2)
    grid on
    XLABEL('Time [sec]')
    YLABEL('Residual force [V]')
    title('Adaptive disturbance rejection')
    
    
    figure
    %hold on
    plot(t,y2);
    hold on
    plot(t,threshold_vec1,'r','LineWidth',2)
    hold on
    plot(t,-threshold_vec1,'r','LineWidth',2)
    hold on
    plot(t,threshold_vec2,'k','LineWidth',2)
    hold on
    plot(t,-threshold_vec2,'k','LineWidth',2)
    hold on
    plot(TD_ver,vec_ver,'-.k','LineWidth',2)
    hold on
    text(TD,max(abs(y2))*3/4,...
        ['Transient duration = ',num2str(TD)],...
        'HorizontalAlignment','center',...
        'BackgroundColor',[.7 .9 .7],...
        'FontSize',16);
    hold on
    text(TD,-max(abs(y2))*3/4,...
        ['Maximum transient value = ',num2str(maximum)],...
        'HorizontalAlignment','center',...
        'BackgroundColor',[.7 .9 .7],...
        'FontSize',16);
    grid on
    XLABEL('Time [sec]')
    YLABEL('Residual force [V]')
    title('Adaptive disturbance rejection')
    xlim([0 2*TD]);
    ylim([-max(abs(y2)) max(abs(y2))]);
end