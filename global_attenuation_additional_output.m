function [G_A, norm_square_ol_in_dB, norm_square_cl_in_dB] =...
    global_attenuation_additional_output (y_ol,y_cl,Fs,ED,Td)


%--------------------------------------------------------------------------
%                   Benchmark on Adaptive Regulation:
%
%   Rejection of unknown/time-varying multiple narrow band disturbances
%
%--------------------------------------------------------------------------
%
%   This function computes the global attenuation for narrow band disturbance
%   rejection (Application to a Benchmark problem) in closed loop with
%   respect to the open loop.
%
%
%
%  G_A = global_attenuation (y_ol,y_cl,Fs,ED,Td)
%
%  G_A : Global attenuation (in dB)
%
%  y_ol : is vector of "open loop" data
%
%  y_cl : is vector of "closed loop" data
%
%  Fs : is the sampling frequency ( Fs = 800 Hz by default)
%
%  Ed : is the experiment duration ( ED = 15 sec by default)
%       The experiment duration should be at least equal to 10 sec. 
%
%  Td : Start disturbance time (Td = 5 sec by default) 
%
%  To compute the global attenuation, the function uses the last three
%  seconds of the records of data in open and in closed loop.
%  
%  The global attenuation is computed with the formula : 
%
%  G_A = 20*log10(norm_ol^2/norm_cl^2);
%  
%  where : 
%  norm_ol^2 is the square of the norm of the output in open loop (last
%  three seconds)
%
%  norm_cl^2 is the square of the norm of the output in closed loop (last
%  three seconds)
%
%  Important Remark : 
%  If you use the default values : Fs=800; ED=15; Td=5; you can
%  run the function with the syntax : 
%
%  G_A = global_attenuation (y_ol,y_cl)
%
%   Written by M. Alma and I.D. Landau (GIPSA-LAB)
%   Version 1, September 22, 2010

if nargin<2, error('Missed open or closed loop data');end
if nargin<3, Fs = 800; ED = 32; Td = 5;end
if nargin<4, Td = 5;ED = 15;end
if nargin<5, Td=5;end


if ED-Td<3, error('Disturbance duration must be at least equal to 3 sec'),end
if Td>=ED, error('Disturbance must be injected before the end of the experiment'),end

y_ol2=y_ol-mean(y_ol);
y_cl2=y_cl-mean(y_cl);

norm_ol = norm(y_ol2(Fs*(ED-3)+1:Fs*ED))^2;
norm_cl = norm(y_cl2(Fs*(ED-3)+1:Fs*ED))^2;

norm_square_ol_in_dB = 20*log10(norm_ol);
norm_square_cl_in_dB = 20*log10(norm_cl);

% norm_ol = var(y_ol2(Fs*(ED-3)+1:Fs*ED));
% norm_cl = var(y_cl2(Fs*(ED-3)+1:Fs*ED));

fprintf ('Global attenuation expressed in "dB" is :')
G_A = 20*log10(norm_ol/norm_cl);

