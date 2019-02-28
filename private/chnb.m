function [nb,channame,strnames] = chnb(channame, varargin)

% chnb() - return channel number corresponding to channel names in an EEG
%           structure
%
%   Optional: initialization (beginning of a script/function):
%   >>                        chnb('labels',labels)
%   Usage:
%   >> [nb]                 = chnb(channameornb);
%   >> [nb,names]           = chnb(channameornb,...);
%   >> [nb,names,strnames]  = chnb(channameornb,...);
%   >> [nb]                 = chnb(channameornb, EEG);
%   >> [nb]                 = chnb(channameornb, labels);
%
% Input:
%   'labels',labels - initialize labels to the list in variable labels for
%                   subsequent executions until clear all is called.
%                   This uses a persistent variable
%   channameornb  - If a string or cell array of strings, it is assumed to
%                   be (part of) the name of channels to search. Either a
%                   string with space separated channel names, or a cell
%                   array of strings. 
%                   Note that regular expressions can be used to match
%                   several channels. See regexp.
%                   If only one channame pattern is given and the string
%                   'inv' is attached to it, the channels NOT matching the
%                   pattern are returned.
%   labels        - Channel names as found in {EEG.chanlocs.labels}.
%
% Output:
%   nb            - Channel numbers in labels, or in the EEG structure
%                   found in the caller workspace (i.e. where the function
%                   is called from) or in the base workspace, if no EEG
%                   structure exists in the caller workspace.
%   names         - Channel names, cell array of strings.
%   strnames      - Channel names, one line character array.
narginchk(1,2);
persistent labels
if nargin == 2
    if ischar(channame)
        labels = varargin{1};
        if strcmp(channame,'labels')
            return
        end
    elseif isstruct(varargin{1}) && isfield(varargin{1},'setname')
        % assume it's an EEG dataset
        labels = {varargin{1}(1).chanlocs.labels};
    else
        labels = varargin{1};
    end
else
    if isempty(labels)
        try
            EEG = evalin('caller','EEG');
        catch
            try
                EEG = evalin('base','EEG');
            catch
                error('Could not find EEG structure');
            end
        end
        if not(isfield(EEG,'chanlocs'))
            error('No channel list found');
        end
        EEG = EEG(1);
        labels = {EEG.chanlocs.labels};
    end
end
if iscell(channame) || ischar(channame)
    
    if ischar(channame) || iscellstr(channame)
        if iscellstr(channame) && numel(channame) == 1 && isempty(channame{1})
            channame = '';
        end
        tmp = regexp(channame,'(\S*) ?','tokens');
        channame = {};
        for i = 1:numel(tmp)
            if iscellstr(tmp{i}{1})
                channame{i} = tmp{i}{1}{1};
            else
                channame{i} = tmp{i}{1};
            end
        end
        if isempty(channame)
            nb = [];
            channame = {};
            strnames = '';
            return
        end
    end
    if not(isempty(strmatch('inv',channame{1})))
        cmd = 'exactinv';
        channame{1} = strrep(channame{1},'inv','');
    else
        channame{1} = channame{1};
        cmd = 'exact';
    end
    % use fieldtrip's shortcuts
    shortcuts = {'meg','megmag','meggrad'};
    isshortcut = ismember(channame,shortcuts);
    if any(isshortcut)
        tmp = {};
        for i_short = 1:numel(isshortcut)
            if ~isshortcut(i_short)
                tmp{end+1} = channame{i_short};
            else
                tmp2 = ft_channelselection(channame{i_short},labels);
                % remove duplicates if any (sorting alphabetically)
                tmp2 = unique(tmp2);
                % undo the sorting, make the order identical to that of the original labels
                [~, indx] = match_str(labels, tmp2);
                tmp2 = tmp2(indx);
                tmp = [tmp; tmp2];
            end
        end
        channame = tmp;
    end
    
    nb = regexpcell(labels,channame,[cmd 'ignorecase']);
    
elseif isnumeric(channame)
    nb = channame;
    if any(nb > numel(labels))
        nb = [];
    end
end
channame = labels(nb);
strnames = sprintf('%s ',channame{:});
if not(isempty(strnames))
    strnames(end) = [];
end
if nargout == 0
    disp(channame)
    disp(nb)
    clear
end


function idx = regexpcell(c,pat, cmds)

% idx = regexpcell(c,pat, cmds)
%
% Return indices idx of cells in c that match pattern(s) pat (regular expression).
% Pattern pat can be char or cellstr. In the later case regexpcell returns
% indexes of cells that match any pattern in pat.
%
% cmds is a string that can contain one or several of these commands:
% 'inv' return indexes that do not match the pattern.
% 'ignorecase' will use regexpi instead of regexp
% 'exact' performs an exact match (regular expression should match the whole strings in c).
% 'all' (default) returns all indices, including repeats (if several pat match a single cell in c).
% 'unique' will return unique sorted indices.
% 'intersect' will return only indices in c that match ALL the patterns in pat.
% 
% v1 Maximilien Chaumon 01/05/09
% v1.1 Maximilien Chaumon 24/05/09 - added ignorecase
% v2 Maximilien Chaumon 02/03/2010 changed input method.
%       inv,ignorecase,exact,combine are replaced by cmds

narginchk(2,3)
if not(iscellstr(c))
    error('input c must be a cell array of strings');
end
if nargin == 2
    cmds = '';
end
if not(isempty(regexpi(cmds,'inv', 'once' )))
    inv = true;
else
    inv = false;
end
if not(isempty(regexpi(cmds,'ignorecase', 'once' )))
    ignorecase = true;
else
    ignorecase = false;
end
if not(isempty(regexpi(cmds,'exact', 'once' )))
    exact = true;
else
    exact = false;
end
if not(isempty(regexpi(cmds,'unique', 'once' )))
    combine = 2;
elseif not(isempty(regexpi(cmds,'intersect', 'once' )))
    combine = 3;
else
    combine = 1;
end

if ischar(pat)
    pat = cellstr(pat);
end

if exact
    for i_pat = 1:numel(pat)
        pat{i_pat} = ['^' pat{i_pat} '$'];
    end
end
    
for i_pat = 1:length(pat)
    if ignorecase
        trouv = regexpi(c,pat{i_pat}); % apply regexp on each pattern
    else
        trouv = regexp(c,pat{i_pat}); % apply regexp on each pattern
    end
    idx{i_pat} = [];
    for i = 1:numel(trouv)
        if not(isempty(trouv{i}))% if there is a match, store index
            idx{i_pat}(end+1) = i;
        end
    end
end
switch combine
    case 1
        idx = [idx{:}];
    case 2
        idx = unique([idx{:}]);
    case 3
        for i_pat = 2:length(pat)
            idx{1} = intersect(idx{1},idx{i_pat});
        end
        idx = idx{1};
end
if inv % if we want to invert result, then do so.
    others = 1:numel(trouv);
    others(idx) = [];
    idx = others;
end
