function list_properties = component_properties(EEG,blink_chans,lpf_band)

% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


list_properties = [];
%
if isempty(EEG.icaweights)
    fprintf('No ICA data.\n');
    return;
end

if ~exist('lpf_band','var') || length(lpf_band)~=2 || ~any(lpf_band)
    ignore_lpf=1;
else
    ignore_lpf=0;
end

delete_activations_after=0;
if ~isfield(EEG,'icaact') || isempty(EEG.icaact)
    delete_activations_after=1;
    EEG.icaact = eeg_getica(EEG);
end
try
    checkfunctionmatlab('pwelch', 'signal_toolbox')
end
for u = 1:size(EEG.icaact,1)
    [spectra(u,:) freqs] = pwelch(EEG.icaact(u,:),[],[],(EEG.srate),EEG.srate);
end

list_properties = zeros(size(EEG.icaact,1),5); %This 5 corresponds to number of measurements made.

for u=1:size(EEG.icaact,1)
    measure = 1;
    % TEMPORAL PROPERTIES
    
    % 1 Median gradient value, for high frequency stuff
    list_properties(u,measure) = median(diff(EEG.icaact(u,:)));
    measure = measure + 1;
    
    % 2 Mean slope around the LPF band (spectral)
    if ignore_lpf
        list_properties(u,measure) = 0;
    else
        list_properties(u,measure) = mean(diff(10*log10(spectra(u,find(freqs>=lpf_band(1),1):find(freqs<=lpf_band(2),1,'last')))));
    end
    measure = measure + 1;
    
    % SPATIAL PROPERTIES
    
    % 3 Kurtosis of spatial map (if v peaky, i.e. one or two points high
    % and everywhere else low, then it's probably noise on a single
    % channel)
    list_properties(u,measure) = kurt(EEG.icawinv(:,u));
    measure = measure + 1;
    
    % OTHER PROPERTIES
    
    % 4 Hurst exponent
    list_properties(u,measure) = hurst_exponent(EEG.icaact(u,:));
    measure = measure + 1;
    
    % 10 Eyeblink correlations
    if (exist('blink_chans','var') && ~isempty(blink_chans))
        for v = 1:length(blink_chans)
            if ~(max(EEG.data(blink_chans(v),:))==0 && min(EEG.data(blink_chans(v),:))==0);
                f = corrcoef(EEG.icaact(u,:),EEG.data(blink_chans(v),:));
                x(v) = abs(f(1,2));
            else
                x(v) = v;
            end
        end
        list_properties(u,measure) = max(x);
        measure = measure + 1;
    end
end

for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
    list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end

if delete_activations_after
    EEG.icaact=[];
end


% The Hurst exponent
%--------------------------------------------------------------------------
% This function does dispersional analysis on a data series, then does a
% Matlab polyfit to a log-log plot to estimate the Hurst exponent of the
% series.
%
% This algorithm is far faster than a full-blown implementation of Hurst's
% algorithm.  I got the idea from a 2000 PhD dissertation by Hendrik J
% Blok, and I make no guarantees whatsoever about the rigor of this approach
% or the accuracy of results.  Use it at your own risk.
%
% Bill Davidson
% 21 Oct 2003

function [hurst] = hurst_exponent(data0)   % data set

data=data0;         % make a local copy

[M,npoints]=size(data0);

yvals=zeros(1,npoints);
xvals=zeros(1,npoints);
data2=zeros(1,npoints);

index=0;
binsize=1;

while npoints>4
    
    y=std(data);
    index=index+1;
    xvals(index)=binsize;
    yvals(index)=binsize*y;
    
    npoints=fix(npoints/2);
    binsize=binsize*2;
    for ipoints=1:npoints % average adjacent points in pairs
        data2(ipoints)=(data(2*ipoints)+data((2*ipoints)-1))*0.5;
    end
    data=data2(1:npoints);
    
end % while

xvals=xvals(1:index);
yvals=yvals(1:index);

logx=log(xvals);
logy=log(yvals);

p2=polyfit(logx,logy,1);
hurst=p2(1); % Hurst exponent is the slope of the linear fit of log-log plot

return;

% NANSUM provides a replacement for MATLAB's nanmean.
%
% For usage see SUM.

function y = nansum(x, dim)

x(isnan(x)) = 0;
if nargin==1
    y = sum(x);
else
    y = sum(x,dim);
end


% NANMEAN provides a replacement for MATLAB's nanmean.
%
% For usage see MEAN.

function y = nanmean(x, dim)

if nargin<2
    N = sum(~isnan(x));
    y = nansum(x) ./ N;
else
    N = sum(~isnan(x), dim);
    y = nansum(x, dim) ./ N;
end

