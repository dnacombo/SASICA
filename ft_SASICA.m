% [cfg] = ft_SASICA(cfg,comp)
%
% This is a wrapper to eeg_SASICA. See help eeg_SASICA for more info.


function [cfg] = ft_SASICA(cfg,comp)


% create a minimal EEG structure to pass to eeg_SASICA.
EEG = [];
EEG.setname = 'work in progress';
EEG.nbchan = numel(comp.topolabel);
EEG.trials = numel(comp.trial);
EEG.pnts = size(comp.trial{1},2);
EEG.srate = comp.fsample;
EEG.xmin = comp.time{1}(1);
EEG.xmax = comp.time{end}(end);
EEG.times = comp.time{1}*1000;
EEG.icaact = cat(3,comp.trial{:});

cfg.layout = ft_prepare_layout(cfg);
EEG.chanlocs = struct('labels',cfg.layout.label);
[EEG.chanlocs.X] = rep2struct(cfg.layout.pos(:,1));
[EEG.chanlocs.Y] = rep2struct(cfg.layout.pos(:,2));
[EEG.chanlocs.Z] = rep2struct(zeros(size(EEG.chanlocs)));

EEG = eeg_SASICA(EEG,cfg);

cfg.reject = EEG.reject;



function [varargout] = rep2struct(varargin)

% [s.target] = rep2struct(dat)
% replicate the value dat into each element of structure s in field target.
% if dat has same nb or elements as s, each element of dat goes into one
% element of s. if dat is more dimensional and doesn't have the same number
% of elements as s, and has same size along dimension 1, then pass each
% slice into s.target.

if numel(varargin) == 1
    dat = varargin{1};
    if numel(dat) == nargout
        for i = 1:nargout
            varargout{i} = dat(i);
        end
    elseif size(dat,1) == nargout
        for i = 1:nargout
            varargout{i} = dat(i,:);
        end
    else
        for i = 1:nargout
            varargout{i} = dat;
        end
    end
elseif numel(varargin) == nargout
    for i = 1:nargout
        varargout{i} = varargin{i};
    end
else
    error('Wrong number of arguments');
end

