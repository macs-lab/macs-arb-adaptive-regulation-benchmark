function varargout = narrowBPDefn_imp(alpha,omega_Hz,Ts,plotflag)
% ==============================================================
% Multiple narrow band-pass filter design via internal model principle
%
% [Q, IMP] = PeakFilterDefn_imp(alpha,omega_Hz,Ts,plotflag)
% 
% Inputs:
%   alpha: shaping coefficient
%   omega_Hz: center frequencies in Hz
%   Ts: sampling time
%   plotflag: 
% 
% Outputs:
%   Q: Desired filter
%   IMP: corresponding internal model
% ============================================================
%   Copyright (c) 2008-, Xu Chen, University of Washington
%   Author(s): Xu Chen
%   Initial version: 2010-09-27
% ============================================================
NBn = length(omega_Hz); % number of passbands

omega_rad = 2*pi*Ts*omega_Hz;

% the internal model
A_element = zeros(NBn,3);
IMP = 1;
for ii = 1:NBn
    A_element(ii,:) = [1 -2*cos(omega_rad(ii)) 1];
    IMP = conv(IMP,A_element(ii,:));
end
% Q = Bq/Aq
Aq = zeros(size(IMP));
for ii = 1:length(IMP)
    Aq(ii) = alpha^(ii-1)*IMP(ii);
end

Bq = Aq - alpha^NBn * IMP;
Q = tf(Bq,Aq,Ts);
    
if nargin == 4
    if strcmp(plotflag, 'plot')
        figure;
        h = bodeplot(Q);
        xlim([1,0.5/Ts]);
        title 'Multiple narrow bandpass filter';
        setoptions(h, 'FreqUnits', 'Hz');grid
    end
end
if nargout == 0
    try
        figure, xbodeplot({Q})
        figure, xbodeplot({1-Q})
    catch
        h = bodeplot(Q);
        xlim([1,0.5/Ts]);
        title 'Multiple narrow bandpass filter';
        setoptions(h, 'FreqUnits', 'Hz');grid 
    end
elseif nargout == 2
    varargout{1} = Bq;
    varargout{2} = Aq;
elseif nargout == 3
    varargout{1} = Bq;
    varargout{2} = Aq;
    varargout{3} = IMP;
elseif nargout == 4
    varargout{1} = Bq;
    varargout{2} = Aq;
    varargout{3} = IMP;
    varargout{4} = Q;
else
    error 'error in the number of outputs' 
end
