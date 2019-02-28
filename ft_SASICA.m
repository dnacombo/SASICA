% [cfg] = ft_SASICA(cfg,comp,data)
%
% Use SASICA selection methods to detect components with specific signal
% properties, as described in
%   Chaumon M, Bishop DV, Busch NA. A Practical Guide to the Selection of
%   Independent Components of the Electroencephalogram for Artifact
%   Correction. Journal of neuroscience methods. 2015
%
% Inputs: - cfg: a structure with field
%               layout: that will be passed to ft_prepare_layout
%              as well as any of the following fields:
%               autocorr:   detect components with low autocorrelation
%               focalcomp   detect focal components in sensor space
%               trialfoc    detect focal components in trial space
%               resvar      Not implemented.
%               SNR         detect components with low signal to noise
%                           ratio across trials between two time windows.
%               EOGcorr     detect components with high correlation with
%                           vertical and horizontal EOG channels
%               chancorr    detect components with high correlation with
%                           any channel
%               FASTER      use FASTER (Nolan et al. 2010) detection
%                           methods.
%               ADJUST      use ADJUST (Mongon et al. 2011) detection
%                           methods
%               MARA        use MARA (Winkler et al. 2011) detection
%                           methods
%               opts        set various options: noplot, nocompute, FontSize
%
%           For more detailed information, see doc eeg_SASICA
%
%       For an example cfg structure, run cfg = ft_SASICA('getdefs')
%
%        - comp: the output of ft_componentanalysis
%        - data: the output of ft_preprocessing
%               it is important that data.elec is populated
%
% output: - cfg structure with additional field reject containing the
%           results of the selection methods
%
%
%   v0 Maximilien Chaumon Feb 2016
%
%   SASICA is a software that helps select independent components of
%   the electroencephalogram based on various signal measures.
%     Copyright (C) 2014  Maximilien Chaumon
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [cfg, EEG] = ft_SASICA(cfg,comp,data)

% removing any eeglab from the path to avoid interference.
rm_frompath eeglab
% adding the local minimal copy of eeglab
addpath(genpath(fullfile(fileparts(mfilename('fullpath')),'eeglab')));

if ischar(cfg)
    cfg = eeg_SASICA([],['EEG = ' cfg ';']);
    return
end
    
% create a minimal EEG structure to pass to eeg_SASICA.
EEG = eeg_emptyset;
EEG.setname = 'internal';
EEG.nbchan = numel(comp.topolabel);
EEG.trials = numel(comp.trial);
EEG.pnts = size(comp.trial{1},2);
EEG.srate = comp.fsample;
EEG.xmin = comp.time{1}(1);
EEG.xmax = comp.time{end}(end);
EEG.times = comp.time{1}*1000;
if numel(unique(cellfun(@(x) size(x,2),comp.trial))) == 1
    EEG.icaact = cat(3,comp.trial{:});
else
    warning('Trials have unequal length. Catenating for display')
    EEG.icaact = cat(2,comp.trial{:});
end
EEG.icawinv = comp.topo;
EEG.icaweights = pinv(EEG.icawinv);
EEG.icasphere  = eye(size(EEG.icaweights,2));


EEG.chanlocs = struct();
if isfield(cfg,'layout');
    cfg.layout = ft_prepare_layout(cfg);
end
for i = 1:EEG.nbchan
    EEG.chanlocs(i).labels = comp.topolabel{i};
    % attempt to create a chanlocs
    if isfield(cfg,'layout')
        ichan = chnb(comp.topolabel{i},cfg.layout.label);
        if ~isempty(ichan)
            [EEG.chanlocs(i).X] = cfg.layout.pos(ichan,1);
            [EEG.chanlocs(i).Y] = cfg.layout.pos(ichan,2);
            [EEG.chanlocs(i).Z] = 1;
        end
    elseif exist('data','var') && isfield(data,'elec')
        ichan = chnb(comp.topolabel{i},data.elec.label);
        if isfield(data.elec,'pnt')
            EEG.chanlocs(i).X = data.elec.pnt(ichan,1);
            EEG.chanlocs(i).Y = data.elec.pnt(ichan,2);
            EEG.chanlocs(i).Z = data.elec.pnt(ichan,3);
        elseif isfield(data.elec,'chanpos')
            EEG.chanlocs(i).X = data.elec.chanpos(ichan,1);
            EEG.chanlocs(i).Y = data.elec.chanpos(ichan,2);
            EEG.chanlocs(i).Z = data.elec.chanpos(ichan,3);
        end            
    
    end
end
EEG.chanlocs = convertlocs(EEG.chanlocs,'cart2all');
mr = max([EEG.chanlocs.radius]);
for i = 1:numel(EEG.chanlocs)
    EEG.chanlocs(i).radius = EEG.chanlocs(i).radius * .5/ mr;
end
EEG.chaninfo.nosedir = '+Y';

if not(exist('data','var'))
    EEG.data = reshape(EEG.icawinv * EEG.icaact(:,:),EEG.nbchan,EEG.pnts,EEG.trials);
else
    if numel(unique(cellfun(@(x) size(x,2),comp.trial))) == 1
        EEG.data = cat(3,data.trial{:});
    else
        EEG.data = cat(2,data.trial{:});
        EEG.trials = 1;
    end
end
EEG.pnts = size(EEG.data,2);
EEG.icachansind = 1:size(EEG.data,1);

if isfield(cfg,'reject')
    EEG.reject = cfg.reject;
end
assignin('base','EEG',EEG);

EEG = eeg_SASICA(EEG,cfg);

assignin('base','EEG',EEG);

cfg.reject = EEG.reject;


