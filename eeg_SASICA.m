% [EEG, cfg] = eeg_SASICA(EEG,cfg)
%
% Suggest components to reject from an EEG dataset with ICA decomposition.
%
% Inputs: EEG: EEGlab structure with ICA fields.
%         cfg: structure describing which methods are to use for suggesting
%              bad components (see structure called def, in the code below)
%              Available methods are:
%              Autocorrelation: detects noisy components with weak
%                               autocorrelation (muscle artifacts usually)
%              Focal components: detects components that are too focal and
%                               thus unlikely to correspond to neural
%                               activity (bad channel or muscle usually).
%              Focal trial activity: detects components with focal trial
%                               activity, with same algorhithm as focal
%                               components above. Results similar to trial
%                               variability.
%              Signal to noise ratio: detects components with weak signal
%                               to noise ratio between arbitrary baseline
%                               and interest time windows.
%              Dipole fit residual variance: detects components with high
%                               residual variance after subtraction of the
%                               forward dipole model. Note that the inverse
%                               dipole modeling using DIPFIT2 in EEGLAB
%                               must have been computed to use this
%                               measure.
%              EOG correlation: detects components whose time course
%                               correlates with EOG channels.
%              Bad channel correlation: detects components whose time course
%                               correlates with any channel(s).
%              ADJUST selection: use ADJUST routines to select components
%                               (see Mognon, A., Jovicich, J., Bruzzone,
%                               L., & Buiatti, M. (2011). ADJUST: An
%                               automatic EEG artifact detector based on
%                               the joint use of spatial and temporal
%                               features. Psychophysiology, 48(2), 229-240.
%                               doi:10.1111/j.1469-8986.2010.01061.x)
%              FASTER selection: use FASTER routines to select components
%                               (see Nolan, H., Whelan, R., & Reilly, R. B.
%                               (2010). FASTER: Fully Automated Statistical
%                               Thresholding for EEG artifact Rejection.
%                               Journal of Neuroscience Methods, 192(1),
%                               152-162. doi:16/j.jneumeth.2010.07.015)
%              MARA selection:  use MARA classification engine to select components
%                               (see Winkler I, Haufe S, Tangermann M.
%                               2011. Automatic Classification of
%                               Artifactual ICA-Components for Artifact
%                               Removal in EEG Signals. Behavioral and
%                               Brain Functions. 7:30.)
%
%              Options: noplot: just compute and store result in EEG. Do
%                           not make any plots.
%
% If you use this program in your research, please cite the following
% article:
%   Chaumon M, Bishop DV, Busch NA. A Practical Guide to the Selection of
%   Independent Components of the Electroencephalogram for Artifact
%   Correction. Journal of neuroscience methods. 2015
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


function [EEG, cfg] = eeg_SASICA(EEG,varargin)

if nargin < 1
    error('Need at least one input argument')
end
if numel(varargin) == 1
    cfg = varargin{1};
elseif numel(varargin) == 0
    cfg = struct;
else
    cfg = vararg2struct(varargin);
end
% deal with calling pop_prop here
if ischar(cfg) && strncmp(cfg,'pop_',4)
    try
        eval(cfg);
    catch ME
        disp('================================');
        disp('================================');
        disp('ERROR. Please send the entire error message below to max.chaumon@gmail.com. Thanks for your help!');
        disp('================================');
        disp(['This is ' eegplugin_SASICA])
        disp(['This is MATLAB ' version])
        disp(['Running on ' computer])
        rethrow(ME)
    end
    return
end
%
PLOTPERFIG = 35;
def = SASICA('getdefs');

cfg = setdef(cfg,def);
v = regexp(version,'^\d+\.\d+','match');
if num2str(v{1}) >= 8.4
    mkersize = 25;
else
    mkersize = 20;
end

if isempty(EEG.icawinv)
    errordlg('No ica weights in the current EEG dataset! Compute ICA on your data first.')
    error('No ica weights! Compute ICA on your data first.')
end
struct2ws(cfg.opts);

rejfields = {'icarejautocorr' 'Autocorrelation' [         0         0    1.0000]
    'icarejfocalcomp' 'Focal components' [         0    0.5000         0]
    'icarejtrialfoc' 'Focal trial activity' [    0.7500         0    0.7500]
    'icarejSNR' 'Signal to noise ' [    0.8000         0         0]
    'icarejresvar' 'Residual variance' [    0     0.7500    0.7500]
    'icarejchancorr' 'Correlation with channels' [    0.7500    0.7500         0]
    'icarejADJUST' 'ADJUST selections' [    .3 .3 .3]
    'icarejFASTER' 'FASTER selections' [    0 .7 0]
    'icarejMARA' 'MARA selections' [    .5 .5 0]
    };

ncomp= size(EEG.icawinv,2); % ncomp is number of components
rejects = zeros(size(rejfields,1),1);

if numel(noplot) == 1
    noplot = noplot * ones(1,size(rejfields,1));
end

if any(~noplot)
    figure(321541);clf;% just a random number so we always work in the same figure
    BACKCOLOR           = [.93 .96 1];
    set(gcf,'numbertitle', 'off','name','Automatic component rejection measures','color',BACKCOLOR)
    isubplot = 1;
end


if ~nocompute
    icaacts = eeg_getdatact(EEG,'component',1:ncomp);
    EEG.icaact = icaacts;
    EEG.reject.SASICA = [];
    for ifield = 1:size(rejfields,1)
        %     EEG.reject.SASICA.(rejfields{ifield}) = false(1,ncomp);
        EEG.reject.SASICA.([rejfields{ifield} 'col']) = rejfields{ifield,3};
    end
    fprintf('Computing selection methods...\n')
end
if cfg.autocorr.enable
    rejects(1) = 1;
    disp('Autocorrelation.')
    %% Autocorrelation
    % Identifying noisy components
    %----------------------------------------------------------------
    struct2ws(cfg.autocorr);

    if ~nocompute
        Ncorrint=round(autocorrint/(1000/EEG.srate)); % number of samples for lag
        rej = false(1,ncomp);
        for k=1:ncomp
            y=icaacts(k,:,:);
            yy=xcorr(mean(y,3),Ncorrint,'coeff');
            autocorr(k) = yy(1);
        end
        dropautocorr = readauto(dropautocorr,autocorr,'-');
        for k = 1:ncomp
            if autocorr(k) < dropautocorr
                rej(k)=true;
            end
        end
        EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','')) = autocorr;
        EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','thresh')) = dropautocorr;
        EEG.reject.SASICA.(rejfields{1,1}) = logical(rej);
    else
        autocorr = EEG.reject.SASICA.(strrep(rejfields{1,1},'rej',''));
        dropautocorr = EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{1,1});
    end
    %----------------------------------------------------------------
    if ~noplot(1)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        plot(autocorr,'k','linestyle','none');

        hold on
        xlim([0 ncomp+1]);
        s = std(autocorr);
        m = mean(autocorr);
        yl = ylim;xl = xlim;
        [x,y] = meshgrid(xl(1):.1:xl(2),yl(1):.1:yl(2));
        galpha = 1./(s*(2*pi)^.5).*exp(-(y-m).^2./(2.*s^2));
        %     h = surf(x,y,-ones(size(y)));shading flat
        %     color = [ 0 0 0]';
        %     C = repmat(color,[1,size(y)]);
        %     C = permute(C,[2 3 1]);
        %     set(h,'alphadata',1-galpha,'alphadatamapping','scaled','facealpha','interp',...
        %         'CData',C,'CDataMapping','direct')
        hline(dropautocorr,'r');

        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{1,3},'markersize',40)
        xlabel('Components')
        ylabel('Autocorrelation')
        title(['Autocorrelation at ' num2str(autocorrint) ' ms.'])
        toplot = autocorr;
        toplot(toplot > dropautocorr) = NaN;
        plot(toplot,'o','color',rejfields{1,3})
        for i = 1:numel(autocorr)
            h = scatter(i,autocorr(i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end

end
if cfg.focalcomp.enable
    rejects(2) = 1;
    disp('Focal components.')
    %% Focal activity
    %----------------------------------------------------------------
    struct2ws(cfg.focalcomp);
    if ~nocompute
        rej = false(1,ncomp);
        clear mywt
        for k=1:ncomp
            mywt(:,k) = sort(abs(zscore(EEG.icawinv(:,k))),'descend'); %sorts standardized weights in descending order
        end
        focalICAout = readauto(focalICAout,mywt(1,:),'+');
        for k = 1:ncomp
            if mywt(1,k) > focalICAout
                rej(k)=true;
            end
        end
        EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','')) = mywt(1,:);
        EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','thresh')) = focalICAout;
        EEG.reject.SASICA.(rejfields{2,1}) = logical(rej);
    else
        mywt(1,:) = EEG.reject.SASICA.(strrep(rejfields{2,1},'rej',''));
        focalICAout = EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{2,1});
    end
    %----------------------------------------------------------------
    if ~noplot(2)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        toplot = mywt(1,:);
        plot(toplot,'k','linestyle','none');
        toplot(toplot < focalICAout) = NaN;
        hold on
        hline(focalICAout,'r');
        plot(toplot,'o','color',rejfields{2,3});
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{2,3},'markersize',40)
        xlabel('Components')
        ylabel('Standardized weights')
        title('Components with focal activity')
        for i = 1:numel(mywt(1,:))
            h = scatter(i,mywt(1,i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end

end

if cfg.trialfoc.enable
    rejects(3) = 1;
    disp('Focal trial activity.');
    %% Focal trial activity
    struct2ws(cfg.trialfoc);
    if ~nocompute
        % Find components with focal trial activity (those that have activity
        % on just a few trials and are almost zero on others)
        %----------------------------------------------------------------
        if ndims(icaacts) < 3
            error('This method cannot be used on continuous data (no ''trials''!)');
        end
        myact =sort(abs(zscore(range(icaacts,2),[],3)),3,'descend'); % sorts standardized range of trial activity
        focaltrialout = readauto(focaltrialout,myact(:,:,1)','+');
        % in descending order
        rej = myact(:,:,1) > focaltrialout;
        EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','')) = myact(:,:,1)';
        EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','thresh')) = focaltrialout;
        EEG.reject.SASICA.(rejfields{3,1}) = rej';
    else
        myact = EEG.reject.SASICA.(strrep(rejfields{3,1},'rej',''))';
        focaltrialout = EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{3,1})';
    end
    %----------------------------------------------------------------
    if ~noplot(3)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        if EEG.trials > 1
            set(gca,'fontsize',FontSize)
            toplot = myact(:,:,1);
            plot(toplot,'k','linestyle','none')
            hold on
            toplot(toplot < focaltrialout) = NaN;
            plot(1:ncomp,toplot,'o','color',rejfields{3,3});
            xlim([0 ncomp+1])
            hline(focaltrialout,'r');
            xl = xlim;yl =ylim;
            xlabel('Components')
            ylabel('Standardized peak trial activity')
            plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{3,3},'markersize',40)
            for i = 1:numel(myact(:,:,1))
                h = scatter(i,myact(i),mkersize,'k','filled');
                cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
                set(h,'buttondownfcn',cb);
            end

            title(['Focal trial activity'])
        else
            xl = xlim;yl = ylim;
            text(xl(1)+diff(xl)/2,yl(1)+diff(yl)/2,{'Only one trial.' 'Focal trial' 'activity method'  'is inadequate.'},'horizontalalignment','center');
            axis off
        end
    end
    %----------------------------------------------------------------
end

if cfg.SNR.enable
    rejects(4) = 1;
    disp('Signal to noise ratio.')
    %% Low Signal to noise components
    struct2ws(cfg.SNR);
    if ~nocompute
        rejfields{4,2} = ['Signal to noise Time of interest ' num2str(snrPOI,'%g ') ' and Baseline ' num2str(snrBL,'%g ') ' ms.'];

        POIpts = timepts(snrPOI);
        BLpts = timepts(snrBL);

        zz = zscore(icaacts,[],2);% zscore along time
        av1 = mean(zz(:,POIpts,:),3); % average activity in POI across trials
        av2 = mean(zz(:,BLpts,:),3); % activity in baseline acros trials
        SNR = std(av1,[],2)./std(av2,[],2); % ratio of the standard deviations of activity and baseline
        snrcut = readauto(snrcut,SNR,'-');
        rej = SNR < snrcut;
        EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','')) = SNR';
        EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','thresh')) = snrcut;
        EEG.reject.SASICA.(rejfields{4,1}) = rej';
    else
        SNR = EEG.reject.SASICA.(strrep(rejfields{4,1},'rej',''))';
        snrcut = EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{4,1})';
    end
    %----------------------------------------------------------------
    if ~noplot(4)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        plot(SNR,'k','linestyle','none');
        hold on
        xlim([0 ncomp+1]);
        xl = xlim; yl = ylim;
        hline(snrcut,'r');
        toplot = SNR;
        toplot(toplot > snrcut) = NaN;
        plot(toplot,'o','color',rejfields{4,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{4,3},'markersize',40)
        for i = 1:numel(SNR)
            h = scatter(i,SNR(i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
        title({'Signal to noise ratio between' ['Time of interest ' num2str(snrPOI,'%g ') ' and Baseline ' num2str(snrBL,'%g ') ' ms.']})
        xlabel('Components')
        ylabel('SNR')
    end

    %----------------------------------------------------------------
end

if cfg.resvar.enable
    rejects(5) = 1;
    disp('Residual variance thresholding.')
    %% High residual variance
    struct2ws(cfg.resvar);
    if ~nocompute
        resvar = 100*[EEG.dipfit.model.rv];
        rej = resvar > thresh;

        EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','')) = resvar;
        EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','thresh')) = thresh;
        EEG.reject.SASICA.(rejfields{5,1}) = rej;
    else
        resvar = EEG.reject.SASICA.(strrep(rejfields{5,1},'rej',''));
        thresh = EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','thresh'));
        rej = EEG.reject.SASICA.(rejfields{5,1});
    end
    %----------------------------------------------------------------
    if ~noplot(5)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        plot(resvar,'k','linestyle','none');
        hold on
        xlim([0 ncomp+1]);
        ylim([0 100]);
        xl = xlim; yl = ylim;
        hline(thresh,'r');
        toplot = resvar;
        toplot(toplot < thresh) = NaN;
        plot(toplot,'o','color',rejfields{5,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{5,3},'markersize',40)
        for i = 1:numel(resvar)
            h = scatter(i,resvar(i),mkersize,'k','filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
        title({'Residual variance of dipole fit'})
        xlabel('Components')
        ylabel('RV (%)')
    end

    %----------------------------------------------------------------
end

if cfg.EOGcorr.enable
    rejects(6) = 1;
    disp('Correlation with EOGs.');
    %% Correlation with EOG
    struct2ws(cfg.EOGcorr);
    if ~nocompute
        noV = 0;noH = 0;
        try
            Veogchan = chnb(Veogchannames);
        catch
            Veogchan = [];
        end
        try
            Heogchan = chnb(Heogchannames);
        catch
            Heogchan = [];
        end
        if numel(Veogchan) == 1
            VEOG = EEG.data(Veogchan,:,:);
        elseif numel(Veogchan) == 2
            VEOG = EEG.data(Veogchan(1),:,:) - EEG.data(Veogchan(2),:,:);
        else
            disp('no Vertical EOG channels...');
            noV = 1;
        end
        if numel(Heogchan) == 1
            HEOG = EEG.data(Heogchan,:,:);
        elseif numel(Heogchan) == 2
            HEOG = EEG.data(Heogchan(1),:,:) - EEG.data(Heogchan(2),:,:);
        else
            disp('no Horizontal EOG channels...');
            noH = 1;
        end
        if noV && noH
            waitfor(errordlg({'No EOG channel names entered.' 'Please enter vertical and/or horizontal channel names' 'or disable "Correlation with EOG" measure and start over.'}))
            return
        end
        ICs = icaacts(:,:)';
        if ~noV
            VEOG = VEOG(:);
            cV  = abs(corr(ICs,VEOG))';
            corthreshV = readauto(corthreshV,cV,'+');
            rejV = cV > corthreshV ;
        else
            cV = NaN(1,size(ICs,2));
            corthreshV = NaN;
            rejV = false(size(cV));
        end
        if ~noH
            HEOG = HEOG(:);
            cH  = abs(corr(ICs,HEOG))';
            corthreshH = readauto(corthreshH,cH,'+');
            rejH = cH > corthreshH;
        else
            cH = NaN(1,size(ICs,2));
            corthreshH = NaN;
            rejH = false(size(cH));
        end

        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'VEOG']) = cV;
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'VEOG']) = corthreshV;
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'HEOG']) = cH;
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'HEOG']) = corthreshH;
        EEG.reject.SASICA.(rejfields{6,1}) = [rejV|rejH];
    else
        if existnotempty(EEG.reject.SASICA,[strrep(rejfields{6,1},'rej','') 'VEOG'])
            cV = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'VEOG']);
            corthreshV = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'VEOG']);
        end
        if existnotempty(EEG.reject.SASICA,[strrep(rejfields{6,1},'rej','') 'HEOG'])
            cH = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'HEOG']);
            corthreshH = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'HEOG']);
        end
    end
    %----------------------------------------------------------------
    if ~noplot(6)
        subplot(2,3,isubplot);cla;isubplot = isubplot+1;
        set(gca,'fontsize',FontSize)
        cols = get(gca,'colororder');
        [hplotcorr] = plot([cV;cH]','.','linestyle','none');
        icol = 2;
        hold all
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;set(gca,'ylimMode','manual');
        hline(corthreshV,'color',cols(1,:));
        hline(corthreshH,'color',cols(2,:));

        title(['Correlation with EOG'])
        legstr = {'VEOG' 'HEOG'};
        legidx = cellfun(@(x)~all(isnan(x)),get(hplotcorr,'ydata'));
        legstr = legstr(legidx);
        ylabel('Correlation coef (r)');
        xlabel('Components');
        toplot = cV;
        toplot(toplot < corthreshV) = NaN;
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        toplot = cH;
        toplot(toplot < corthreshH) = NaN;
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
        if not(sum(legidx) == 0)
            legend(hplotcorr(legidx),legstr,'fontsize',10, 'location', 'best');
        end
        for i = 1:numel(cH)
            h(1) = scatter(i,cV(i),mkersize,cols(1,:),'filled');
            h(2) = scatter(i,cH(i),mkersize,cols(2,:),'filled');
            cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end
    %----------------------------------------------------------------
end

if cfg.chancorr.enable
    rejects(6) = 1;
    disp('Correlation with other channels.')
    %% Correlation with other channels
    struct2ws(cfg.chancorr);
    try
        [chan cellchannames channames] = chnb(channames);
    end
    if ~nocompute
        if ~cfg.EOGcorr.enable
            rejH = false(1,ncomp);
            rejV = false(1,ncomp);
        end
        if ~isempty(channames)
            chanEEG = EEG.data(chan,:)';
            ICs = icaacts(:,:)';
            c  = abs(corr(ICs,chanEEG))';
            corthresh = mean(readauto(corthresh,c,'+'));
            rej = c > corthresh ;
            if size(rej,1) > 1
                rej = sum(rej)>=1;
            end
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']) = c;
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'chans']) = corthresh;
            EEG.reject.SASICA.(rejfields{6,1}) = [rej|rejH|rejV];
        else
            noplot(6) = 1;
            disp('Could not find the channels to compute correlation.');
            c = NaN(1,ncomp);
            corthresh = mean(readauto(corthresh,c,'+'));
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']) = c;
            EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'chans']) = corthresh;
            rej = false(1,ncomp);
            EEG.reject.SASICA.(rejfields{6,1}) = [rej|rejV|rejH];
        end
    else
        c = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']);
        corthresh = EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','thresh') 'chans']);
    end
    %----------------------------------------------------------------
    if ~noplot(6);
        if exist('hplotcorr','var')
            isubplot = isubplot-1;
        end
        subplot(2,3,isubplot);
        if ~cfg.EOGcorr.enable
            cla;
            set(gca,'fontsize',FontSize);
            cols = get(gca,'colororder');
        end
        hold all
        if not(exist('hplotcorr','var'))
            hplotcorr = [];
        end
        icol = numel(hplotcorr);
        for ichan = 1:numel(chan)
            [hplotcorr(end+1)] = plot([c(ichan,:)]','.','linestyle','none','color',cols(rem(icol+ichan-1,size(cols,1))+1,:));
        end
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;set(gca,'ylimMode','manual');
        hline(corthresh,'r');
        title(['Correlation with channels'])
        if cfg.EOGcorr.enable
            legstr = {'VEOG' 'HEOG' cellchannames{:}};
        else
            legstr = {cellchannames{:}};
        end
        if not(isempty(hplotcorr))
            legidx = cellfun(@(x)~all(isnan(x)),get(hplotcorr,'ydata'));
        else legidx = [];
        end
        legstr = legstr(legidx);
        
        ylabel('Correlation coef (r)');
        xlabel('Components');
        toplot = c;
        for i = 1:size(toplot,1)
            toplot(i,toplot(i,:) < corthresh) = NaN;
        end
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
        if not(sum(legidx) == 0)
            legend(hplotcorr(legidx),legstr,'fontsize',10, 'location', 'best');
        end
        for ichan = 1:size(c,1)
            for i = 1:size(c,2)
                h = scatter(i,c(ichan,i),mkersize,cols(rem(icol+ichan-1,size(cols,1))+1,:),'filled');
                cb = sprintf('eeg_SASICA(EEG, ''pop_prop( %s, 0, %d, findobj(''''tag'''',''''comp%d''''), { ''''freqrange'''', [1 50] })'');', inputname(1), i, i);
                set(h,'buttondownfcn',cb);
            end
        end

    end
    %----------------------------------------------------------------
end
if cfg.ADJUST.enable
    rejects(7) = 1;
    disp('ADJUST methods selection')
    %% ADJUST
    struct2ws(cfg.ADJUST);
    if ~nocompute
        [art, horiz, vert, blink, disc,...
            soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
            soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin] = ADJUST (EEG);

        ADJ.art = art;ADJ.horiz = horiz;ADJ.vert = vert;ADJ.blink = blink;ADJ.disc = disc;

        ADJ.soglia_DV = soglia_DV; ADJ.diff_var = diff_var;
        ADJ.soglia_K = soglia_K;ADJ.med2_K = med2_K; ADJ.meanK = meanK;
        ADJ.soglia_SED = soglia_SED; ADJ.med2_SED = med2_SED;ADJ.SED = SED;
        ADJ.med2_SAD = med2_SAD;ADJ.soglia_SAD = soglia_SAD;ADJ.SAD = SAD;
        ADJ.soglia_GDSF = soglia_GDSF; ADJ.med2_GDSF = med2_GDSF;ADJ.GDSF = GDSF;
        ADJ.soglia_V = soglia_V;ADJ.med2_V = med2_V;ADJ.nuovaV = nuovaV;
        ADJ.soglia_D = soglia_D; ADJ.maxdin = maxdin;

        rej = false(1,size(EEG.icaact,1));
        rej([ADJ.art ADJ.horiz ADJ.vert ADJ.blink ADJ.disc]) = true;

        EEG.reject.SASICA.(strrep(rejfields{7,1},'rej','')) = ADJ;
        EEG.reject.SASICA.(rejfields{7,1}) = rej;
    else
        ADJ = EEG.reject.SASICA.(strrep(rejfields{7,1},'rej',''));
    end
    %----------------------------------------------------------------
end
if cfg.FASTER.enable
    rejects(8) = 1;
    disp('FASTER methods selection')
    %% FASTER
    struct2ws(cfg.FASTER);
    if ~nocompute
        blinkchanname = chnb(blinkchanname);
        if isempty(blinkchanname)
            waitfor(errordlg({'If you use the FASTER method,' 'it is highly recommended to provide a "blink channel".' 'The method is unreliable without it.' 'Please provide a blink channel or disable the "FASTER" method.'}))
            return
        end
        listprops = component_properties(EEG,blinkchanname);
        FST.rej = min_z(listprops)' ~= 0;
        FST.listprops = listprops;

        EEG.reject.SASICA.(strrep(rejfields{8,1},'rej','')) = FST;
        EEG.reject.SASICA.(rejfields{8,1}) = FST.rej;
    else
        FST = EEG.reject.SASICA.(strrep(rejfields{8,1},'rej',''));
    end
    %----------------------------------------------------------------
end
if cfg.MARA.enable
    rejects(9) = 1;
    disp('MARA methods selection')
    %% MARA
    struct2ws(cfg.MARA);
    if ~nocompute
        [rej info] = MARA(EEG);
        MR.rej = false(1,size(EEG.icaact,1));
        MR.rej(rej) = true;
        MR.info = info;

        EEG.reject.SASICA.(strrep(rejfields{9,1},'rej','')) = MR;
        EEG.reject.SASICA.(rejfields{9,1}) = MR.rej;
    else
        MR = EEG.reject.SASICA.(strrep(rejfields{9,1},'rej',''));
    end
    %----------------------------------------------------------------
end

EEG.reject.SASICA.var = var(EEG.icaact(:,:),[],2);% variance of each component

if (cfg.ADJUST.enable||cfg.FASTER.enable) && any(~noplot)
    h = uicontrol('style','text','string','for ADJUST or FASTER results, click on component buttons in the other window(s)','units','normalized','position',[0 0 1 .05],'backgroundcolor',get(gcf,'color'));
    uistack(h,'bottom')
end
fprintf('... Done.\n')

drawnow

%% Final computations
% combine in gcompreject field and pass to pop_selectcomps
EEG.reject.gcompreject = false(1,ncomp);
for ifield = 1:size(rejfields,1)
    if rejects(ifield)
        EEG.reject.gcompreject = [EEG.reject.gcompreject ; EEG.reject.SASICA.(rejfields{ifield})];
    end
end
EEG.reject.gcompreject = sum(EEG.reject.gcompreject) >= 1;

%% plotting
try
    delete(findobj('-regexp','name','pop_selectcomps'))
    drawnow
end
if any(~noplot)
    if ~isempty([EEG.chanlocs.radius])% assume we have sensor locations...
        clear hfig
        delete(findobj('tag','waitcomp'))
        textprogressbar;
        textprogressbar('Drawing topos...');
        for ifig = 1:ceil((ncomp)/PLOTPERFIG)
            cmps = [1+(ifig-1)*PLOTPERFIG:min([ncomp,ifig*PLOTPERFIG])];
            eeg_SASICA(EEG,['pop_selectcomps(EEG, [' num2str(cmps) '],' num2str(ncomp) ');']);
            hfig(ifig) = gcf;
            set(hfig(ifig),'name',[get(hfig(ifig),'name') ' -- SASICA ' num2str(ifig)]);
            % find the ok button and change its callback fcn
            okbutt = findobj(hfig(ifig),'string','OK');
            set(okbutt,'callback',['delete(findobj(''-regexp'',''name'',''pop_selectcomps.* -- SASICA''));delete(findobj(''-regexp'',''name'',''Automatic component rejection measures''));' ...
                'if exist(''ALLEEG'',''var'') && exist(''EEG'',''var'') && exist(''CURRENTSET'',''var''); [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,CURRENTSET); if not(isempty(findobj(''-regexp'',''name'',''^EEGLAB''))); eeglab(''redraw'');end;end;' ...
                'warndlg({''Remember you need to now subtract the marked components.'' ''Use Tools > Remove components''});']);
            % find the cancel button and change its callback fcn
            cancelbutt = findobj(hfig(ifig),'string','Cancel');
            closecallback = ['try; delete(findobj(''-regexp'',''name'',''pop_selectcomps''));delete(findobj(''-regexp'',''name'',''Automatic component rejection measures''));end;'];
            set(cancelbutt,'callback',[closecallback 'EEG.reject.gcompreject = false(size(EEG.reject.gcompreject));disp(''Operation cancelled. No component is selected for rejection.'');']);
            set(hfig(ifig),'closerequestfcn',closecallback)
            % crazy thing to find and order the axes for the topos.
            ax{ifig} = findobj(hfig(ifig),'type','Axes');
            ax{ifig} = ax{ifig}(end-1:-1:1);% erase pointer to the big axis behind all others and reorder the axes handles.
        end;
        ax = vertcat(ax{:});

        if not(numel(ax) == ncomp) || isempty(okbutt) || ~ishandle(okbutt)
            errordlg('Please do not click while I''m drawing these topos, it''s disturbing. Start over again...')
            error('Please do not click while I''m drawing these topos, it''s disturbing. Start over again...')
        end

        % create markers next to each topoplot showing which threshold has been
        % passed.
        for i_comp = 1:ncomp
            if EEG.reject.gcompreject(i_comp)
                %                 axes(ax(i_comp))
                f = get(ax(i_comp),'parent');
                set(0,'currentFigure',f);
                set(f,'CurrentAxes',ax(i_comp));
                drawnow;
                hold on
                for irej = 1:size(rejfields,1)
                    if isfield(EEG.reject.SASICA,rejfields{irej,1}) && ...
                            EEG.reject.SASICA.(rejfields{irej,1})(i_comp)
                        x = -.5 + (irej > 6);
                        y = .5 - .1*irej-.3*(rem(irej-1,6)+1>3);
                        h = plot(x,y,'markerfacecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'marker','o');
                    end
                end
            end
        end
        set(hfig,'visible','on');
        try
            pop_selectcomps(EEG, [ncomp+1]);
        end
        textprogressbar;
        hlastfig = gcf;
        set(hlastfig,'name',[get(hlastfig,'name') ' -- SASICA']);
        lastax = findobj(hlastfig,'type','Axes');
        set(lastax,'visible','off');
        axes(lastax(end));
        hold on
        for irej = 1:numel(rejects)
            set(gca,'xlimmode','manual');
            if rejects(irej)
                x = 0;
                y = .5 - .2*irej;

                scatter(x,y,'markerfacecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.SASICA.([rejfields{irej} 'col']));
                text(x+.1,y,[rejfields{irej,2} ' (' num2str(sum(EEG.reject.SASICA.(rejfields{irej,1}))) ')']);
            end
        end
        for i = numel(hfig):-1:1
            figure(hfig(i));
            setctxt(hfig(i),EEG,cfg);
        end
        figure(hlastfig);
    else
        disp('No channel locations. I''m not plotting.');
    end
end
if nargout == 0
    assignin('caller','EEG',EEG);
end


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

function thresh = readauto(thresh,dat,comp)
% if thresh starts with 'auto'
% compute auto threshold as mean(dat) +/- N std(dat)
% with N read in the string thresh = 'auto N'
% if not, use thresh as a value
if isstr(thresh) && strncmp(thresh,'auto',4)
    if numel(thresh) > 4
        threshsigma = str2num(thresh(5:end));
    else
        threshsigma = 2;
    end
    thresh = eval(['mean(dat,2)' comp 'threshsigma * std(dat,[],2)']);
end



function setctxt(hfig,EEG,cfg)
COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
buttons = findobj(hfig,'-regexp','tag','^comp\d{1,3}$');
buttonnums = regexp(get(buttons,'tag'),'comp(\d{1,3})','tokens');
if numel(buttonnums)>1
    buttonnums = cellfun(@(x)(str2num(x{1}{1})),buttonnums);
else
    buttonnums = str2num(buttonnums{1}{1});
end
for i = 1:numel(buttonnums)
    hcmenu = uicontextmenu;

    if ~isempty(EEG.reject.gcompreject)
        status = EEG.reject.gcompreject(buttonnums(i));
    else
        status = 0;
    end;

    hcb1 = ['EEG.reject.gcompreject(' num2str(buttonnums(i)) ') = ~EEG.reject.gcompreject(' num2str(buttonnums(i)) ');'...
        'set(gco,''backgroundcolor'',fastif(EEG.reject.gcompreject(' num2str(buttonnums(i)) '), ' COLREJ ',' COLACC '));'...
        'set(findobj(''tag'',''ctxt' num2str(buttonnums(i)) '''), ''Label'',fastif(EEG.reject.gcompreject(' num2str(buttonnums(i)) '),''ACCEPT'',''REJECT''));' ];
    uimenu(hcmenu, 'Label', fastif(status,'ACCEPT','REJECT'), 'Callback', hcb1,'tag',['ctxt' num2str(buttonnums(i))]);

    mycb = strrep(get(buttons(i),'Callback'),'''','''''');
    mycb = regexprep(mycb,'pop_prop','eeg_SASICA(EEG,''pop_prop');
    mycb = [mycb ''');'];
    set(buttons(i),'CallBack',mycb)
    set(buttons(i),'uicontextmenu',hcmenu)
end

function s = setdef(s,d)
% s = setdef(s,d)
% Merges the two structures s and d recursively.
% Adding the default field values from d into s when not present or empty.

if isstruct(s) && not(isempty(s))
    fields = fieldnames(d);
    for i_f = 1:numel(fields)
        if isfield(s,fields{i_f})
            s.(fields{i_f}) = setdef(s.(fields{i_f}),d.(fields{i_f}));
        else
            s.(fields{i_f}) = d.(fields{i_f});
        end
    end
elseif not(isempty(s))
    s = s;
elseif isempty(s);
    s = d;
end

function res = existnotempty(s,f)
res = isfield(s,f) && not(isempty(s.(f)));
