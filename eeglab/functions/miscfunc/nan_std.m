% nan_std() - std, not considering NaN values
%
% Usage: std across the first dimension

% Author: Arnaud Delorme, CNL / Salk Institute, Sept 2003

% Copyright (C) 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

function out = nan_std(in)
    
    if nargin < 1
      help nan_std;
      return;
    end;
    
    nans = find(isnan(in));
    in(nans) = 0;
   
    nonnans = ones(size(in));
    nonnans(nans) = 0;
    nonnans = sum(nonnans);
    nononnans = find(nonnans==0);
    nonnans(nononnans) = NaN;
   
    out = sqrt((sum(in.^2)-sum(in).^2./nonnans)./(nonnans-1));
    out(nononnans) = NaN;
