function [norm_residual,norm_after_transient,ratio,percentage]=...
    transient_evaluation(y,Fs,DR,Td)
%--------------------------------------------------------------------------
%                   Benchmark on Adaptive Regulation:
%
%   Rejection of unknown/time-varying multiple narrow band disturbances
%
%--------------------------------------------------------------------------
%
%     This function evaluates the transient behavior for narrow band 
%   disturbance rejection (Application to a Benchmark problem).
%   
%   The transient behavior is considered as the necessary time for the
%  square of the truncated two norm to reach a value
%    <or = 1.21 times the square of the truncated 
%   two norm in steady state. The square of the truncated two norm in  
%  steady state is obtained using the last 3 seconds of the test 
%   before the disturbance is removed (once the algorithm has converged).
%   Since in the benchmark specifications  a transient duration
%   equal or less than 2 seconds is imposed, the square of the truncated
%   norm is evaluated over an horizon of 3 sec starting after the first 2
%   seconds of the test  and compared with the square of
%   the truncated two norm correponsing to the steday state behavior.
%    The quadratic
%   truncated norm will be used for this analysis. 
%
%   quadratic truncated norm = norm (x) ^ 2;
%   where x is a vector containing the values of the samples.
%
%   The time periods to be used are:
%   Begining of the test is at 5 seconds
%   For the transient behavior from 7 to 10 seconds
%   For the residual behavior (steady state) from 17 to 20 seconds
%   In addition, the ratio between these quantities will be expressed as a 
%    percentage with respect
%   to the optimum behavior (a ratio equal or less than 1.21 will
%   correspond to 100%). 
%
%   Written by A. Castellanos Silva and I.D. Landau (GIPSA-LAB)
%   Version 1, March 12, 2013
%
%
%   y : is the vector of data
%   Fs : is the sampling frequency ( Fs = 800 Hz by default)
%   DR : is the disturbance remotion time ( DR = 20 sec by default)
%   The experiment duration should be at least equal to 20 sec. 
%   Td : Start disturbance time (Td = 5 sec by default)
%   
%   norm_residual : the square of the truncated two norm of the time period
%   between 17 and 20 seconds.
%   norm_after_transient : the square of the truncated two norm of the time
%   period between 7 and 10 seconds.
%   ratio : the ratio between the two norms above defined.
%   percentage : the percentage of BSI defined for this criterion.
%
if nargin<2, Fs = 800; DR = 20; Td = 5;end
if nargin<3, DR = 20; Td = 5;end
if nargin<4, Td = 5;end

if Td>=DR, error('Disturbance must be injected before the end of the experiment'),end

y=y-mean(y);

y_residual = y((DR-3)*Fs:(DR*Fs)-1);
norm_residual = norm (y_residual) ^ 2;

y_after_transient = y((Td+2)*Fs:((Td+5)*Fs)-1);
norm_after_transient = norm (y_after_transient) ^ 2;

ratio = norm_after_transient /  norm_residual; 
limit_benchmark = 2.42;
if ratio >= limit_benchmark
    percentage = 0;
elseif ratio < limit_benchmark && ratio > 1.21
    Delta_ratio = ratio - 1.21;
    percentage = ((1.21 - Delta_ratio)/1.21)*100;
else
    percentage = 100;
end