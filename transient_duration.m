function TD = transient_duration(y,Fs,ED,Td)


%--------------------------------------------------------------------------
%                   Benchmark on Adaptive Regulation:
%
%   Rejection of unknown/time-varying multiple narrow band disturbances
%
%--------------------------------------------------------------------------
%
%   This function computes the transient duration for narrow band disturbance
%   rejection (Application to a Benchmark problem)
%   
%   The transient is the time necessary to the output to reach the value
%   +/-2.17sigma , where sigma is the standard deviation for the residual output
%   (once the algo has converged) computed during the 3 last seconds of the experiment. 
%
%   sigma (x) = sqrt ( var(x))
%   where "x" is the vector of data corresponding to the 3 last seconds of
%   the experiment.
%
%   In addition the output after has reached for
%   the first time 2.17 sigma should remain within +/- 4 sigma. 
%   The test has been developped assuming that the residual output is a gaussian
%   process.(ie 97% of values are between +/-2.17 sigma, and 99.99% of values are
%   between +/- 4 sigma).
%
%   However it is tolerated that 0.1% of the measurements be outside +/- 4
%   sigma. The percentage is denoted Dp. (to take in account that it is not
%   a pure gaussian variable.
%
%
%   Written by M. Alma and I.D. Landau (GIPSA-LAB)
%   Version 1, September 6, 2010
%
%
%  TD = transient_duration(y,Fs,ED,Td)
%  TD : transient duration (in seconds)
%
%  y : is vector of data
%  Fs : is the sampling frequency ( Fs = 800 Hz by default)
%  Ed : is the exepriment duration ( ED = 15 sec by default)
%       The experiment duration should be at least equal to 10 sec. 
%  Td : Start disturbance time (Td = 5 sec by default)
%  threshold 1 = 2.17 * standard deviation (after convergence)
%  threshold 2 = 4 * standard deviation (after convergence)

if nargin<2, Fs = 800; ED = 15; Td = 5;end
if nargin<3, ED = 15; Td = 5;end
if nargin<4, Td = 5;end

if Td>=ED, error('Disturbance must be injected before the end of the experiment'),end

y=y-mean(y);

yres = y((ED-3)*Fs+1:ED*Fs);
perc = 99.9;


threshold1 = 2.17*sqrt(var(yres));
threshold2 = 4*sqrt(var(yres));
l=0;
for n = 1:length(y)
    if abs(y(n)) <= threshold1
        g=0;
        for k = n:length(y)
        if abs(y(k)) > threshold2
        g= g+1;
        end
        end
        if g<(length(y)-n+1)*((100-perc)/100)
            l=n;
        break
        end
    end
end
g;
fprintf ('Percentage of points above +/- 4 sigma is :')
Dp = g/(length(y)-n+1) *100


if l==0, error('Algorithm does not converge'),end

fprintf ('The transient maximum value is:')
maximum = max(abs(y))

fprintf ('Transient duration is :')
TD=(l+1)/Fs - Td;

t= 0:1/800:ED;

threshold_vec1=threshold1*ones(1,length(y));

threshold_vec2=threshold2*ones(1,length(y));

vec_ver = -max(abs(y)):0.001:max(abs(y));
TD_ver =  (l+1)/Fs*ones(1,length(vec_ver));

figure
plot(t,y);
hold on
plot(t,threshold_vec1,'r','LineWidth',2)
hold on
plot(t,-threshold_vec1,'r','LineWidth',2)
hold on
plot(t,threshold_vec2,'k','LineWidth',2)
hold on
plot(t,-threshold_vec2,'k','LineWidth',2)
grid on
xlabel('Time [sec]')
ylabel('Residual force [V]')
 title('Adaptive disturbance rejection')


figure

%hold on
plot(t,y);
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
text(Td+TD,max(abs(y))*3/4,...
	['Transient duration = ',num2str(TD)],...
	'HorizontalAlignment','center',... 
	'BackgroundColor',[.7 .9 .7],...
    'FontSize',16);
hold on
text(Td+TD,-max(abs(y))*3/4,...
	['Maximum transient value = ',num2str(maximum)],...
	'HorizontalAlignment','center',... 
	'BackgroundColor',[.7 .9 .7],...
    'FontSize',16);
grid on
xlabel('Time [sec]')
ylabel('Residual force [V]')
title('Adaptive disturbance rejection')
xlim([Td-0.1 Td+2*TD]);
ylim([-max(abs(y)) max(abs(y))]);

 