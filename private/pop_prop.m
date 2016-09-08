% pop_prop() - plot the properties of a channel or of an independent
%              component.
% Usage:
%   >> pop_prop( EEG);           % pops up a query window
%   >> pop_prop( EEG, typecomp); % pops up a query window
%   >> pop_prop( EEG, typecomp, chanorcomp, winhandle,spectopo_options);
%
% Inputs:
%   EEG        - EEGLAB dataset structure (see EEGGLOBAL)
%
% Optional inputs:
%   typecomp   - [0|1] 1 -> display channel properties
%                0 -> component properties {default: 1 = channel}
%   chanorcomp - channel or component number[s] to display {default: 1}
%
%   winhandle  - if this parameter is present or non-NaN, buttons
%                allowing the rejection of the component are drawn.
%                If non-zero, this parameter is used to back-propagate
%                the color of the rejection button.
%   spectopo_options - [cell array] optional cell arry of options for
%                the spectopo() function.
%                For example { 'freqrange' [2 50] }
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: pop_runica(), eeglab()

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

% hidden parameter winhandle

% 01-25-02 reformated help & license -ad
% 02-17-02 removed event index option -ad
% 03-17-02 debugging -ad & sm
% 03-18-02 text settings -ad & sm
% 03-18-02 added title -ad & sm

function pop_prop(EEG, typecomp, chanorcomp, winhandle, spec_opt)


% assumed input is chanorcomp
% -------------------------
try, icadefs;
catch,
    BACKCOLOR = [0.8 0.8 0.8];
    GUIBUTTONCOLOR   = [0.8 0.8 0.8];
end;
basename = ['Component ' int2str(chanorcomp) ];

fh = figure('name', ['pop_prop() - ' basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'on');
pos = get(gcf,'position');
set(gcf,'Position', [pos(1) pos(2)-700+pos(4) 500 700], 'visible', 'on');
pos = get(gca,'position'); % plot relative to current axes
hh = gca;
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
delete(gca);
p = panel();
p.margin = [10 10 10 10];
p.pack('v',{.35 []});
p(1).margin = [0 0 0 0];
p(1).pack('h',{.4 [] .01});


% plotting topoplot
p(1,1).select()
topoplot( EEG.icawinv(:,chanorcomp), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
    'shading', 'interp', 'numcontour', 3); axis square;
title(basename, 'fontsize', 14);

% plotting erpimage
p(1,2).margin = [15 15 5 15];
p(1,2).select();
eeglab_options;
if EEG.trials > 1
    % put title at top of erpimage
    axis off
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
        ei_smooth = 1;
    else
        ei_smooth = 3;
    end
    icaacttmp = eeg_getdatact(EEG, 'component', chanorcomp);
    offset = nan_mean(icaacttmp(:));
    era    = nan_mean(squeeze(icaacttmp)')-offset;
    era_limits=get_era_limits(era);
    erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times*1000, ...
        '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','erp_vltg_ticks',era_limits);
    title(sprintf('%s activity \\fontsize{10}(global offset %3.3f)', basename, offset));
else
    % put title at top of erpimage
    EI_TITLE = 'Continous data';
    ERPIMAGELINES = 200; % show 200-line erpimage
    while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
        ERPIMAGELINES = 0.9 * ERPIMAGELINES;
    end
    ERPIMAGELINES = round(ERPIMAGELINES);
    if ERPIMAGELINES > 2   % give up if data too small
        if ERPIMAGELINES < 10
            ei_smooth = 1;
        else
            ei_smooth = 3;
        end
        erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
        erpimageframestot = erpimageframes*ERPIMAGELINES;
        eegtimes = linspace(0, erpimageframes-1, EEG.srate/1000);
        if typecomp == 1 % plot channel
            offset = nan_mean(EEG.data(chanorcomp,:));
            % Note: we don't need to worry about ERP limits, since ERPs
            % aren't visualized for continuous data
            erpimage( reshape(EEG.data(chanorcomp,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar');
        else % plot component
            icaacttmp = eeg_getdatact(EEG, 'component', chanorcomp);
            offset = nan_mean(icaacttmp(:));
            erpimage(reshape(icaacttmp(:,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset,ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel', '');
        end
    else
        axis off;
        text(0.1, 0.3, [ 'No erpimage plotted' 10 'for small continuous data']);
    end;
end;

% plotting spectrum
% -----------------
if ~exist('winhandle','var') || isempty(winhandle) || ~ishandle(winhandle)
    winhandle = NaN;
    p(2).pack('v',{.3 [] })
else
    p(2).pack('v',{.3 [] .1})
end;
p(2,1).pack('h',{.01,[],.01});
p(2,1).margin = [15 15 0 55];
p(2,1,1).margin = 0;
p(2,1,3).margin = 0;
p(2,1,2).pack('v',{.01 []});
p(2,1,2,1).margin = 0;
p(2,1,2,2).margintop = 5;
p(2,1,2,2).select();
try
    spectopo( EEG.icaact(chanorcomp,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
    % set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 12);
    % set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 12);
    axis on
    xlabel('Frequency (Hz)')
    h = title('Activity power spectrum', 'fontsize', 10);
    %     set(h,'position',get(h,'position')+[-15 -7 0]);
    set(gca,'fontSize',10)
catch
    axis off;
    lasterror
    text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;

%%%% Add SASICA measures.
%               eye      muscle/noise    channel     ~ok
colors = { [0 .75 .75]      [0 0 1]      [0 .5 0] [.2 .2 .2]};
% C={[1 0 0],[.6 0 .2],[1 1 0],[0 1 0], [0 1 1]};% colors used in ADJ
computed = fieldnames(EEG.reject.SASICA);
computed = computed(regexpcell(computed,'rej|thresh|^var$','inv'));
computedthresh = regexprep(computed,'ica','icathresh');
computedrej = regexprep(computed,'ica','icarej');
toPlot = {};
toPlot_axprops = {};
toPlot_title = {}; SXticks = {};co = [];
for i = 1:numel(computed)
    if strcmp(computed{i},'icaADJUST')
        struct2ws(EEG.reject.SASICA.icaADJUST)
        toPlot{end+1}{1} = (SAD(chanorcomp)-med2_SAD)/(soglia_SAD-med2_SAD);
        toPlot{end}{2} = (SED(chanorcomp)-med2_SED)/(soglia_SED-med2_SED);
        toPlot{end}{3} = (GDSF(chanorcomp)-med2_GDSF)/(soglia_GDSF-med2_GDSF);
        toPlot{end}{4} = (nuovaV(chanorcomp)-med2_V)/(soglia_V-med2_V);
        toPlot{end}{5} = (meanK(chanorcomp)-med2_K)/(soglia_K-med2_K);
        ADJis = '';
        aco = repmat(colors{4},numel(toPlot{end}),1);
        if ismember(chanorcomp,horiz)
            ADJis = [ADJis 'HEM/'];
            aco([2 4],:) = repmat(colors{1},2,1);
        end
        if ismember(chanorcomp,vert)
            ADJis = [ADJis 'VEM/'];
            aco([1 4],:) = repmat(colors{1},2,1);
        end
        if ismember(chanorcomp,blink)
            ADJis = [ADJis 'Blink/'];
            aco([1 4 5],:) = repmat(colors{1},3,1);
        end
        if ismember(chanorcomp,disc)
            ADJis = [ADJis 'Disc/'];
            aco([3 4],:) = repmat(colors{3},2,1);
        end
        if isempty(ADJis)
            ADJis = 'OK';
        else
            ADJis(end) = [];
        end
        toPlot_title{end+1} = ['ADJUST: ' ADJis];
        toPlot_axprops{end+1} = {'ColorOrder' aco,...
            'ylim' [0 2],...
            'ytick' [1 2],...
            'yticklabel' {'Th' '2*Th'},...
            'xtick' 1:numel(toPlot{end}),...
            'xticklabel' {'SAD' 'SED' 'GDSF' 'MEV' 'TK'}};
    elseif strcmp(computed{i},'icaFASTER')
        listprops = EEG.reject.SASICA.icaFASTER.listprops;
        str='FASTER: ';
        FASTER_reasons = {'HighFreq ' 'FlatSpectrum ' 'SpatialKurtosis ' 'HurstExponent ' 'EOGCorrel '};
        %                     1 Median gradient value, for high frequency stuff
        %                     2 Mean slope around the LPF band (spectral)
        %                     3 Kurtosis of spatial map
        %                     4 Hurst exponent
        %                     5 Eyeblink correlations
        zlist = zscore(listprops);
        for i = 1:size(listprops,2)
            fst(:,i) = min_z(listprops(:,i));
        end
        reasons = FASTER_reasons(fst(chanorcomp,:));
        if isempty(reasons)
            str = [str 'OK'];
        else
            str = [str reasons{:}];
        end
        FSTis = str;
        toPlot{end+1} = {};
        for ip = 1:numel(zlist(chanorcomp,:))
            toPlot{end}{ip} = abs(zlist(chanorcomp,ip))/3;% normalized by threshold
        end
        toPlot_title{end+1} = FSTis;
        toPlot_axprops{end+1} = {'ColorOrder' [colors{2};colors{2};colors{3};colors{2};colors{1}],...
            'ylim' [0 2],...
            'ytick' [1 2],...
            'yticklabel' {'Th' '2*Th'},...
            'xtick',1:numel(toPlot{end}),...
            'xticklabel',{'MedGrad' 'SpecSl' 'SK' 'HE' 'EOGCorr'}};
    elseif strcmp(computed{i},'icaMARA')
        info = EEG.reject.SASICA.icaMARA.info;
        str='MARA: ';
        MARA_meas = {'CurrDensNorm ' 'SpatRange ' 'AvgLocSkew ' '\lambda ' '8-13 Pow' '1/F Fit '};
        %                     1 Current Density Norm
        %                     2 Spatial Range
        %                     3 Average Local Skewness
        %                     4 lambda
        %                     5 Band Power (8-13 Hz)
        %                     6 Fit Error
        if ~ EEG.reject.SASICA.icarejMARA(chanorcomp)
            str = [str 'OK       '];
        else
            str = [str 'Reject    '];
        end
        MARAis = [str '(' num2str(round(100*info.posterior_artefactprob(chanorcomp)),'%g') '%)'];
        toPlot{end+1} = {};
        for ip = 1:numel(info.normfeats(:,chanorcomp))
            toPlot{end}{ip} = info.normfeats(ip,chanorcomp) ;
        end
        toPlot_title{end+1} = MARAis;
        toPlot_axprops{end+1} = {'ColorOrder' repmat(colors{4},numel(MARA_meas),1),...
            'ylimmode' 'auto',...
            'xtick',1:numel(toPlot{end}),...
            'xticklabel',{'CDN' 'SpRg' 'AvLocSkw' 'lambda' '8-13 Hz' '1/F Fit'}
            };
    else
        rejfields = {
            'icaautocorr'       'LoAC'   colors{2}
            'icafocalcomp'      'FocCh'       colors{3}
            'icatrialfoc'       'FocTr'        colors{3}
            'icaSNR'            'LoSNR'        colors{2}
            'icaresvar'         'ResV'          colors{2}
            'icachancorrVEOG'   'CorrV'         colors{1}
            'icachancorrHEOG'   'CorrH'         colors{1}
            'icachancorrchans'  'CorrC'         colors{3}
            };
        if isempty(toPlot)
            toPlot{1} = {};
            toPlot_axprops{1} = {};
            toPlot_title{1} = 'SASICA';
        end
        switch computed{i}
            case 'icaautocorr'
                toPlot{1}{end+1} = 2 - (EEG.reject.SASICA.(computed{i})(chanorcomp) +1)/(EEG.reject.SASICA.(computedthresh{i}) +1);
            case 'icaSNR'
                toPlot{1}{end+1} = EEG.reject.SASICA.(computedthresh{i})/EEG.reject.SASICA.(computed{i})(chanorcomp);
            otherwise
                toPlot{1}{end+1} = EEG.reject.SASICA.(computed{i})(:,chanorcomp)/EEG.reject.SASICA.(computedthresh{i});
        end
        SXticks{end+1} = rejfields{strcmp(computed{i},rejfields(:,1)),2};
        co(end+1,:) = rejfields{strcmp(computed{i},rejfields(:,1)),3};
    end
end
if not(isempty(SXticks))
    toPlot_axprops{1} = {toPlot_axprops{1}{:} 'ylim' [0 2]...
        'ytick' [1 2] ...
        'yticklabel' {'Th' '2*Th'} ...
        'xtick' 1:numel(SXticks) ...
        'Xticklabel' SXticks...
        'xlim',[.5 numel(SXticks)+.5],...
        'colororder',co};
end

p(2,2).pack('v',numel(toPlot));
p(2,2).de.margintop = 0;
for i = 1:numel(toPlot)
    p(2,2,i).pack('h',{.2 []});
    p(2,2,i,1).select();
    text(1.1,0.5,strjust(strwrap(toPlot_title{i},15),'right'),'horizontalalignment','right');
    axis off
    p(2,2,i,2).select()
    hold on
    set(gca,toPlot_axprops{i}{:});
    cs = get(gca,'colorOrder');
    for j = 1:numel(toPlot{i})
        xj = linspace(j-(numel(toPlot{i}{j})>1)*.3,j+(numel(toPlot{i}{j})>1)*.3,numel(toPlot{i}{j}));
        bar(xj,toPlot{i}{j},'facecolor',cs(rem(j-1,size(cs,1))+1,:));
    end
    hline(1,':k')
end


% display buttons
% ---------------
if ishandle(winhandle)
    COLREJ = '[1 0.6 0.6]';
    COLACC = '[0.75 1 0.75]';
    % CANCEL button
    % -------------
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

    %     % VALUE button
    %     % -------------
    %     hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

    % REJECT button
    % -------------
    if ~isempty(EEG.reject.gcompreject)
        status = EEG.reject.gcompreject(chanorcomp);
    else
        status = 0;
    end;
    hr = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
        'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
    command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
        'if get(gcbo, ''userdata''),' ...
        '     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
        'else ' ...
        '     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
        'end;' ];
    set( hr, 'callback', command);

    %     % HELP button
    %     % -------------
    %     h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop'');');

    % OK button
    % ---------
    command = [ 'global EEG;' ...
        'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
        'EEG.reject.gcompreject(' num2str(chanorcomp) ') = tmpstatus;' ];
    if winhandle ~= 0
        command = [ command ...
            sprintf('if tmpstatus set(gcbo, ''backgroundcolor'', %s); else set(gcbo, ''backgroundcolor'', %s); end;', ...
            COLREJ, COLACC) ...
            ['obj = findobj(''-regexp'',''name'',''pop_selectcomps.* -- SASICA''); obj = fastif(isempty(obj),[],findobj(obj,''tag'',''comp' num2str(chanorcomp) '''));'] ...
            sprintf('if ~isempty(obj) && tmpstatus set(obj, ''backgroundcolor'', %s); else set(obj, ''backgroundcolor'', %s); end;', ...
            COLREJ, COLACC)];
    end;
    command = [ command 'close(gcf); clear tmpstatus' ];
    h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

    %     % draw the figure for statistical values
    %     % --------------------------------------
    %     index = num2str( chanorcomp );
    %     command = [ ...
    %         'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
    %         '' ...
    %         'pos = get(gcf,''Position'');' ...
    %         'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
    %         'pos = get(gca,''position'');' ...
    %         'q = [pos(1) pos(2) 0 0];' ...
    %         's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
    %         'axis off;' ...
    %         ''  ...
    %         'txt1 = sprintf(''(\n' ...
    %         'Entropy of component activity\t\t%2.2f\n' ...
    %         '> Rejection threshold \t\t%2.2f\n\n' ...
    %         ' AND                 \t\t\t----\n\n' ...
    %         'Kurtosis of component activity\t\t%2.2f\n' ...
    %         '> Rejection threshold \t\t%2.2f\n\n' ...
    %         ') OR                  \t\t\t----\n\n' ...
    %         'Kurtosis distibution \t\t\t%2.2f\n' ...
    %         '> Rejection threhold\t\t\t%2.2f\n\n' ...
    %         '\n' ...
    %         'Current thesholds sujest to %s the component\n\n' ...
    %         '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
    %         'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
    %         'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
    %         '' ...
    %         'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
    %         'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
    %         'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
    %         'clear txt1 q s h pos;' ];
    %     set( hval, 'callback', command);
    %     if isempty( EEG.stats.compenta )
    %         set(hval, 'enable', 'off');
    %     end;

    com = sprintf('pop_prop( %s, %d, %d, 0, %s);', inputname(1), typecomp, chanorcomp, vararg2str( { spec_opt } ) );
else
    com = sprintf('pop_prop( %s, %d, %d, NaN, %s);', inputname(1), typecomp, chanorcomp, vararg2str( { spec_opt } ) );
end;

return;

function era_limits=get_era_limits(era)
%function era_limits=get_era_limits(era)
%
% Returns the minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)
%
% Inputs:
% era - [vector] Event related activation or potential
%
% Output:
% era_limits - [min max] minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)

mn=min(era);
mx=max(era);
mn=orderofmag(mn)*round(mn/orderofmag(mn));
mx=orderofmag(mx)*round(mx/orderofmag(mx));
era_limits=[mn mx];


function [lengths]  =  min_z(list_properties, rejection_options)
if (~exist('rejection_options', 'var'))
    rejection_options.measure = ones(1, size(list_properties, 2));
    rejection_options.z = 3*ones(1, size(list_properties, 2));
end

rejection_options.measure = logical(rejection_options.measure);
zs = list_properties - repmat(mean(list_properties, 1), size(list_properties, 1), 1);
zs = zs./repmat(std(zs, [], 1), size(list_properties, 1), 1);
zs(isnan(zs)) = 0;
all_l  =  abs(zs) > repmat(rejection_options.z, size(list_properties, 1), 1);
lengths  =  any(all_l(:, rejection_options.measure), 2);

function out = nan_mean(in)

nans = find(isnan(in));
in(nans) = 0;
sums = sum(in);
nonnans = ones(size(in));
nonnans(nans) = 0;
nonnans = sum(nonnans);
nononnans = find(nonnans==0);
nonnans(nononnans) = 1;
out = sum(in)./nonnans;
out(nononnans) = NaN;


function ord=orderofmag(val)
%function ord=orderofmag(val)
%
% Returns the order of magnitude of the value of 'val' in multiples of 10
% (e.g., 10^-1, 10^0, 10^1, 10^2, etc ...)
% used for computing erpimage trial axis tick labels as an alternative for
% plotting sorting variable

val=abs(val);
if val>=1
    ord=1;
    val=floor(val/10);
    while val>=1,
        ord=ord*10;
        val=floor(val/10);
    end
    return;
else
    ord=1/10;
    val=val*10;
    while val<1,
        ord=ord/10;
        val=val*10;
    end
    return;
end


function tw = strwrap(t,n)

% tw = strwrap(t,n)
%
% wrap text array t at n characters taking non alphanumeric characters as
% breaking characters (i.e. not cutting words strangely).

t = deblank(t(:)');
seps = '[\s-]';
tw = '';
while not(isempty(t))
    breaks = regexp(t,seps);
    breaks(end+1) = numel(t);
    idx = 1:min(n,breaks(find(breaks < n, 1,'last')));
    if isempty(idx)
        idx = 1:min(n,numel(t));
    end
    tw(end+1,:) =  [ t( idx ) repmat( char( 32 ) , [1 n - numel( idx ) ] ) ];
    t(idx)= [];
    t = strtrim(t);
end



