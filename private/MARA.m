%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   END FASTER CODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   BEGIN MARA CODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MARA() - Automatic classification of multiple artifact components
%          Classies artifactual ICs based on 6 features from the time domain,
%           the frequency domain, and the pattern
%
% Usage:
%   >> [artcomps, info] = MARA(EEG);
%
% Inputs:
%   EEG         - input EEG structure
%
% Outputs:
%   artcomps    - array containing the numbers of the artifactual
%                 components
%   info        - struct containing more information about MARA classification
%                   .posterior_artefactprob : posterior probability for each
%                            IC of being an artefact
%                   .normfeats : <6 x nIC > features computed by MARA for each IC,
%                            normalized by the training data
%                      The features are: (1) Current Density Norm, (2) Range
%                      in Pattern, (3) Local Skewness of the Time Series,
%                      (4) Lambda, (5) 8-13 Hz, (6) FitError.
%
%  For more information see:
%  I. Winkler, S. Haufe, and M. Tangermann, Automatic classification of artifactual ICA-components
%  for artifact removal in EEG signals, Behavioral and Brain Functions, 7, 2011.
%
% See also: processMARA()

% Copyright (C) 2013 Irene Winkler and Eric Waldburger
% Berlin Institute of Technology, Germany
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
function [artcomps, info] = MARA(EEG)
%%%%%%%%%%%%%%%%%%%%
%%  Calculate features from the pattern (component map)
%%%%%%%%%%%%%%%%%%%%
% extract channel labels
clab = {};
for i=1:length(EEG.icachansind) %i=1:length(EEG.chanlocs), NF edit
    clab{i} = EEG.chanlocs(EEG.icachansind(i)).labels; %EEG.chanlocs(i).labels; NF edit
end

% cut to channel labels common with training data
load('fv_training_MARA'); %load struct fv_tr
[clab_common i_te i_tr ] = intersect(upper(clab), upper(fv_tr.clab));
clab_common = fv_tr.clab(i_tr);
if length(clab_common) == 0
    error(['There were no matching channeldescriptions found.' , ...
        'MARA needs channel labels of the form Cz, Oz, F3, F4, Fz, etc. Aborting.'])
end
patterns = (EEG.icawinv(i_te,:));
[M100 idx] = get_M100_ADE(clab_common); %needed for Current Density Norm

disp('MARA is computing features. Please wait');
%standardize patterns
patterns = patterns./repmat(std(patterns,0,1),length(patterns(:,1)),1);

%compute current density norm
feats(1,:) = log(sqrt(sum((M100*patterns(idx,:)).^2)));
%compute spatial range
feats(2,:) = log(max(patterns) - min(patterns));

%%%%%%%%%%%%%%%%%%%%
%%  Calculate time and frequency features
%%%%%%%%%%%%%%%%%%%%
%compute time and frequency features (Current Density Norm, Range Within Pattern,
%Average Local Skewness, Band Power 8 - 13 Hz)
feats(3:6,:) = extract_time_freq_features(EEG);
disp('Features ready');


%%%%%%%%%%%%%%%%%%%%%%
%%  Adapt train features to clab
%%%%%%%%%%%%%%%%%%%%
fv_tr.pattern = fv_tr.pattern(i_tr, :);
fv_tr.pattern = fv_tr.pattern./repmat(std(fv_tr.pattern,0,1),length(fv_tr.pattern(:,1)),1);
fv_tr.x(2,:) = log(max(fv_tr.pattern) - min(fv_tr.pattern));
fv_tr.x(1,:) = log(sqrt(sum((M100 * fv_tr.pattern).^2)));

%%%%%%%%%%%%%%%%%%%%
%%  Classification
%%%%%%%%%%%%%%%%%%%%
[C, foo, posterior] = classify(feats',fv_tr.x',fv_tr.labels(1,:));
artcomps = find(C == 0)';
info.posterior_artefactprob = posterior(:, 1)';
info.normfeats = (feats - repmat(mean(fv_tr.x, 2), 1, size(feats, 2)))./ ...
    repmat(std(fv_tr.x,0, 2), 1, size(feats, 2));

