% internal function of SASICA (adapted from eeglab)

% pop_subcomp() - remove specified components from an EEG dataset.
%                 and subtract their activities from the data. Else,
%                 remove components already marked for rejection.
% Usage:
%   >> OUTEEG = pop_subcomp( INEEG ); % pop-up window mode
%   >> OUTEEG = pop_subcomp( INEEG, components, confirm);
%
% Pop-up window interface:
%   "Component(s) to remove ..." - [edit box] Array of components to
%                remove from the data. Sets the 'components' parameter
%                in the command line call (see below).
%   "Component(s) to retain ..." - [edit box] Array of components to
%                to retain in the data. Sets the 'components' parameter in
%                the command line call. Then, comp_to_remove = ...
%                    setdiff([1:size(EEG.icaweights,1)], comp_to_keep)
%                Overwrites "Component(s) to remove" (above).
% Command line inputs:
%   INEEG      - Input EEG dataset.
%   components - Array of components to remove from the data. If empty,
%                 remove components previously marked for rejection (e.g.,
%                 EEG.reject.gcompreject).
%   confirm    - [0|1] Display the difference between original and processed
%                dataset. 1 = Ask for confirmation. 0 = Do not ask. {Default: 0}
% Outputs:
%   OUTEEG     - output dataset.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: compvar()

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% 01-25-02 reformated help & license -ad
% 02-15-02 propagate ica weight matrix -ad sm jorn

function [EEG, com] = pop_subcomp( EEG, components, plotag )

com='';


component_keep = setdiff_bc(1:size(EEG.icaweights,1), components);
compproj = EEG.icawinv(:, component_keep)*eeg_getdatact(EEG, 'component', component_keep, 'reshape', '2d');
compproj = reshape(compproj, size(compproj,1), EEG.pnts, EEG.trials);

%fprintf( 'The ICA projection accounts for %2.2f percent of the data\n', 100*varegg);

ButtonName = 'continue';
while ~strcmpi(ButtonName, 'Close')
    ButtonName=questdlg( [ 'What would you like to plot?' ], ...
        'Choose', 'Close', 'ERPs', 'single trials','single trials');
    if strcmpi(ButtonName, 'ERPs')
        if EEG.trials > 1
            tracing  = [ squeeze(mean(EEG.data(EEG.icachansind,:,:),3)) squeeze(mean(compproj,3))];
            figure;
            plotdata(tracing, EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000 0 0], ...
                'Trial ERPs with (red) and without (blue) these components');
        else
            uiwait(warndlg('Cannot plot ERPs for continuous data'));
        end;
    elseif strcmpi(ButtonName, 'single trials')
        eegplot( EEG.data(EEG.icachansind,:,:), 'srate', EEG.srate, 'title', 'Black = channel before rejection; red = after rejection -- eegplot()', ...
            'limits', [EEG.xmin EEG.xmax]*1000, 'data2', compproj);
        uiwait(gcf);
    end;
end;

return;
