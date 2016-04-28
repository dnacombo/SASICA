% computeSED_NOnorm() - Computes Spatial Eye Difference feature
% without normalization
%
% Usage:
%   >> [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   out        - SED values
%   medie_left - Left Eye area average values
%   medie_right- Right Eye area average values
%
%
% Copyright (C) 2009-2014 Andrea Mognon (1) and Marco Buiatti (2),
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
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

function [out,medie_left,medie_right]=computeSED_NOnorm(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Left Eye area (LE)
dimleft=0; %number of LE electrodes
index1=zeros(1,nchannels); %indexes of LE electrodes

for k=1:nchannels
    if (-61<chanlocs(1,k).theta) && (chanlocs(1,k).theta<-35) && (chanlocs(1,k).radius>0.30) %electrodes are in LE
        dimleft=dimleft+1; %count electrodes
        index1(1,dimleft)=k;
    end
end

% Find electrodes in Right Eye area (RE)
dimright=0; %number of RE electrodes
index2=zeros(1,nchannels); %indexes of RE electrodes
for g=1:nchannels
    if (34<chanlocs(1,g).theta) && (chanlocs(1,g).theta<61) && (chanlocs(1,g).radius>0.30) %electrodes are in RE
        dimright=dimright+1; %count electrodes
        index2(1,dimright)=g;
    end
end

% Find electrodes in Posterior Area (PA)
dimback=0;
index3=zeros(1,nchannels);
for h=1:nchannels
    if (abs(chanlocs(1,h).theta)>110)
        dimback=dimback+1;
        index3(1,dimback)=h;
    end
end

if dimleft*dimright*dimback==0
    disp('ERROR: no channels included in some scalp areas.')
    disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
    disp('ADJUST session aborted.')
    return
end

%% Outputs

out=zeros(1,n); %memorizes SED
medie_left=zeros(1,n); %memorizes LE mean value
medie_right=zeros(1,n); %memorizes RE mean value

%% Output computation

for i=1:n  % for each topography
    %create LE electrodes vector
    left=zeros(1,dimleft);
    for h=1:dimleft
        left(1,h)=topog(i,index1(1,h));
    end

    %create RE electrodes vector
    right=zeros(1,dimright);
    for h=1:dimright
        right(1,h)=topog(i,index2(1,h));
    end

    %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end



    %compute features
    out1=abs(mean(left)-mean(right));
    out2=var(back);
    out(1,i)=out1; % SED not notmalized
    medie_left(1,i)=mean(left);
    medie_right(1,i)=mean(right);


end


