% EEG = eeg_SASICA(EEG,cfg)
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
%
%              Options: NbMin: minimum number of the above methods that
%                           must be passed in order to select component.
%                       noplot: just compute and store result in EEG. Do
%                           not make any plots.
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

% Some of the measures used here are based on http://bishoptechbits.blogspot.com/2011/05/automated-removal-of-independent.html

function [EEG, cfg] = eeg_SASICA(EEG,cfg)

if nargin < 1
    error('Need at least one input argument')
end
% deal with calling pop_prop_ADJ or pop_prop_FST here
if ischar(cfg) && strncmp(cfg,'pop_prop',8)
    eval(cfg);
    return
end
if ~exist('cfg','var')
    cfg = struct;
end
%
PLOTPERFIG = 35;
def = SASICA('getdefs');

cfg = setdef(cfg,def);

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
    };
rejects = zeros(size(rejfields,1),1);

if numel(noplot) == 1
    noplot = noplot * ones(1,size(rejfields,1));
end

if any(~noplot)
    figure(321541);clf;% just a random number so we always work in the same figure
    BACKCOLOR           = [.93 .96 1];
    set(gcf,'numbertitle', 'off','name','Automatic component rejection measures','color',BACKCOLOR)
end


ncomp= size(EEG.icawinv,2); % ncomp is number of components
icaacts = eeg_getdatact(EEG,'component',1:ncomp);
EEG.icaact = icaacts;
EEG.reject.SASICA = [];
for ifield = 1:size(rejfields,1)
%     EEG.reject.SASICA.(rejfields{ifield}) = false(1,ncomp);
    EEG.reject.SASICA.([rejfields{ifield} 'col']) = rejfields{ifield,3};
end
fprintf('Computing selection methods...\n')
if cfg.autocorr.enable
    rejects(1) = 1;
    disp('Autocorrelation.')
    %% Autocorrelation
    % Identifying noisy components
    %----------------------------------------------------------------
    struct2ws(cfg.autocorr);

    Ncorrint=round(autocorrint/(1000/EEG.srate)); % number of samples for lag
    rej = false(1,ncomp);
    for k=1:ncomp
        y=icaacts(k,:,:);
        yy=xcorr(mean(y,3),Ncorrint,'coeff');
        autocorr(k) = yy(1);
        if yy(1) < dropautocorr
            rej(k)=true;
        end
    end
    EEG.reject.SASICA.(strrep(rejfields{1,1},'rej','')) = autocorr;
    EEG.reject.SASICA.(rejfields{1,1}) = logical(rej);
    %----------------------------------------------------------------
    if ~noplot(1)
        subplot(2,3,1);cla
        set(gca,'fontsize',FontSize)
        plot(autocorr,'k');

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
        plot(xl,[dropautocorr dropautocorr],'r');

        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{1,3},'markersize',40)
        xlabel('Components')
        ylabel('Autocorrelation')
        title(['Autocorrelation at ' num2str(autocorrint) ' ms.'])
        toplot = autocorr;
        toplot(toplot > dropautocorr) = NaN;
        plot(toplot,'o','color',rejfields{1,3})
        for i = 1:numel(autocorr)
            h = scatter(i,autocorr(i),20,'k','filled');
            cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
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
    rej = false(1,ncomp);
    clear mywt
    for k=1:ncomp
        mywt(:,k) = sort(abs(zscore(EEG.icawinv(:,k))),'descend'); %sorts standardized weights in descending order
        if mywt(1,k) > focalICAout
            rej(k)=true;
        end
    end
    EEG.reject.SASICA.(strrep(rejfields{2,1},'rej','')) = mywt(1,:);
    EEG.reject.SASICA.(rejfields{2,1}) = logical(rej);
    %----------------------------------------------------------------
    if ~noplot(2)
        subplot(2,3,2);cla
        set(gca,'fontsize',FontSize)
        toplot = mywt(1,:);
        plot(toplot,'k');
        toplot(toplot < focalICAout) = NaN;
        hold on
        plot(xlim,[focalICAout focalICAout],'r');
        plot(toplot,'o','color',rejfields{2,3});
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{2,3},'markersize',40)
        xlabel('Components')
        ylabel('Standardized weights')
        title('Components with focal activity')
        for i = 1:numel(mywt(1,:))
            h = scatter(i,mywt(1,i),20,'k','filled');
            cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
            set(h,'buttondownfcn',cb);
        end
    end

end

if cfg.trialfoc.enable
    rejects(3) = 1;
    disp('Focal trial activity.');
    %% Focal trial activity
    % Find components with focal trial activity (those that have activity
    % on just a few trials and are almost zero on others)
    %----------------------------------------------------------------
    if ndims(icaacts) < 3
        error('This method cannot be used on continuous data (no ''trials''!)');
    end
    struct2ws(cfg.trialfoc);
    myact =sort(abs(zscore(range(icaacts,2),[],3)),3,'descend'); % sorts standardized range of trial activity
    % in descending order
    rej = myact(:,:,1) > focaltrialout;
    EEG.reject.SASICA.(strrep(rejfields{3,1},'rej','')) = myact(:,:,1)';
    EEG.reject.SASICA.(rejfields{3,1}) = rej';

    %----------------------------------------------------------------
    if ~noplot(3)
        subplot(2,3,3);cla
        if EEG.trials > 1
            set(gca,'fontsize',FontSize)
            hold on
            toplot = myact(:,:,1);
            plot(toplot,'k')
            toplot(toplot < focaltrialout) = NaN;
            plot(1:ncomp,toplot,'o','color',rejfields{3,3});
            xlim([0 ncomp+1])
            plot(xlim,[focaltrialout focaltrialout],'r');
            xl = xlim;yl =ylim;
            xlabel('Components')
            ylabel('Standardized peak trial activity')
            plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{3,3},'markersize',40)
            for i = 1:numel(myact(:,:,1))
                h = scatter(i,myact(i),20,'k','filled');
                cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
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
    rejfields{4,2} = ['Signal to noise Time of interest ' num2str(snrPOI,'%g ') ' and Baseline ' num2str(snrBL,'%g ') ' ms.'];

    POIpts = timepts(snrPOI);
    BLpts = timepts(snrBL);

    zz = zscore(icaacts,[],2);% zscore along time
    av1 = mean(zz(:,POIpts,:),3); % average activity in POI across trials
    av2 = mean(zz(:,BLpts,:),3); % activity in baseline acros trials
    SNR = std(av1,[],2)./std(av2,[],2); % ratio of the standard deviations of activity and baseline
    rej = SNR < snrcut;
    EEG.reject.SASICA.(strrep(rejfields{4,1},'rej','')) = SNR';
    EEG.reject.SASICA.(rejfields{4,1}) = rej';

    %----------------------------------------------------------------
    if ~noplot(4)
        subplot(2,3,4);cla
        set(gca,'fontsize',FontSize)
        plot(SNR,'k');
        hold on
        xlim([0 ncomp+1]);
        xl = xlim; yl = ylim;
        plot(xl,[snrcut snrcut],'r');
        toplot = SNR;
        toplot(toplot > snrcut) = NaN;
        plot(toplot,'o','color',rejfields{4,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{4,3},'markersize',40)
        for i = 1:numel(SNR)
            h = scatter(i,SNR(i),20,'k','filled');
            cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
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

    resvar = 100*[EEG.dipfit.model.rv];
    rej = resvar > thresh;

    EEG.reject.SASICA.(strrep(rejfields{5,1},'rej','')) = resvar;
    EEG.reject.SASICA.(rejfields{5,1}) = rej;

    %----------------------------------------------------------------
    if ~noplot(5)
        subplot(2,3,5);cla
        set(gca,'fontsize',FontSize)
        plot(resvar,'k');
        hold on
        xlim([0 ncomp+1]);
        ylim([0 100]);
        xl = xlim; yl = ylim;
        plot(xl,[thresh thresh],'r');
        toplot = resvar;
        toplot(toplot < thresh) = NaN;
        plot(toplot,'o','color',rejfields{5,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{5,3},'markersize',40)
        for i = 1:numel(resvar)
            h = scatter(i,resvar(i),20,'k','filled');
            cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
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
    ICs = icaacts(:,:)';
    if ~noV
        VEOG = VEOG(:);
        cV  = abs(corr(ICs,VEOG))';
        rejV = cV > corthreshV ;
    else
        cV = NaN(1,size(ICs,2));
        rejV = false(size(cV));
    end
    if ~noH
        HEOG = HEOG(:);
        cH  = abs(corr(ICs,HEOG))';
        rejH = cH > corthreshH;
    else
        cH = NaN(1,size(ICs,2));
        rejH = false(size(cH));
    end

    EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'VEOG']) = cV;
    EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'HEOG']) = cH;
    EEG.reject.SASICA.(rejfields{6,1}) = [rejV|rejH];

    %----------------------------------------------------------------
    if ~noplot(6)
        subplot(2,3,6);cla
        set(gca,'fontsize',FontSize)
        [hplotcorr] = plot([cV;cH]');
        hold all
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        plot(xlim,[corthreshV corthreshV ],'r');
        plot(xlim,[corthreshH corthreshH ],'m');

        title(['Correlation with EOG'])
        legstr = {'VEOG' 'HEOG'};
        ylabel('Correlation coef (r)');
        xlabel('Components');
        toplot = cV;
        toplot(toplot < corthreshV) = NaN;
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        toplot = cH;
        toplot(toplot < corthreshH) = NaN;
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
        legend(legstr,'fontsize',10, 'location', 'best');
        for i = 1:numel(cH)
            h(1) = scatter(i,cH(i),20,'k','filled');
            h(2) = scatter(i,cV(i),20,'k','filled');
            cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
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
    if ~cfg.EOGcorr.enable
        rejH = false(1,ncomp);
        rejV = false(1,ncomp);
    end
    if ~isempty(channames)
        try
            [chan cellchannames channames] = chnb(channames);
        end
        chanEEG = EEG.data(chan,:)';
        ICs = icaacts(:,:)';
        c  = abs(corr(ICs,chanEEG))';
        rej = c > corthresh ;
        if size(rej,1) > 1
            rej = sum(rej)>=1;
        end
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']) = c;
        EEG.reject.SASICA.(rejfields{6,1}) = [rej|rejH|rejV];
    else
        noplot(6) = 1;
        disp('Could not find the channels to compute correlation.');
        c = NaN(1,ncomp);
        EEG.reject.SASICA.([strrep(rejfields{6,1},'rej','') 'chans']) = c;
        rej = false(1,ncomp);
        EEG.reject.SASICA.(rejfields{6,1}) = [rej|rejV|rejH];
    end

    %----------------------------------------------------------------
    if ~noplot(6);
        subplot(2,3,6);
        if ~cfg.EOGcorr.enable
            cla;
            set(gca,'fontsize',FontSize);
        end
        hold all
        if not(exist('hplotcorr'))
            hplotcorr = [];
        end
        [hplotcorr(end+1:end+numel(chan))] = plot([c]');
        xlim([0 ncomp+1]);
        xl = xlim;yl = ylim;
        plot(xlim,[corthresh corthresh ],'b');
        title(['Correlation with channels'])
        if cfg.EOGcorr.enable
            legstr = {'VEOG' 'HEOG' cellchannames{:}};
        else
            legstr = {cellchannames{:}};
        end
        ylabel('Correlation coef (r)');
        xlabel('Components');
        toplot = c;
        for i = 1:size(toplot,1)
            toplot(i,toplot(i,:) < corthresh) = NaN;
        end
        plot(1:ncomp,toplot,'o','color',rejfields{6,3})
        plot(xl(2)-diff(xl)/20,yl(2)-diff(yl)/20,'marker','.','color',rejfields{6,3},'markersize',40)
        legend(hplotcorr,legstr,'fontsize',10, 'location', 'best');
        for ichan = 1:size(c,1)
            for i = 1:size(c,2)
                h = scatter(i,c(ichan,i),20,'k','filled');
                cb = sprintf('pop_prop( %s, 0, %d, findobj(''tag'',''comp%d''), { ''freqrange'', [1 50] });', inputname(1), i, i);
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

%----------------------------------------------------------------
end
if cfg.FASTER.enable
    rejects(8) = 1;
    disp('FASTER methods selection')
    %% FASTER
    struct2ws(cfg.FASTER);
    blinkchans = chnb(blinkchans);
    listprops = component_properties(EEG,blinkchans);
    FST.rej = min_z(listprops)' ~= 0;
    FST.listprops = listprops;

    EEG.reject.SASICA.(strrep(rejfields{8,1},'rej','')) = FST;
    EEG.reject.SASICA.(rejfields{8,1}) = FST.rej;


    %----------------------------------------------------------------
end
if cfg.ADJUST.enable||cfg.FASTER.enable
    uicontrol('style','text','string','for ADJUST or FASTER results, right click on component buttons in the other window(s)','units','normalized','position',[0 0 1 .05],'backgroundcolor',get(gcf,'color'));
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
EEG.reject.gcompreject = sum(EEG.reject.gcompreject) >= NbMin;

%% plotting
try
    delete(findobj('-regexp','name','pop_selectcomps'))
end
if any(~noplot)
    if ~isempty([EEG.chanlocs.radius])% assume we have sensor locations...
        clear hfig
        for ifig = 1:ceil((ncomp)/PLOTPERFIG)
            pop_selectcomps(EEG, [1+(ifig-1)*PLOTPERFIG:min([ncomp,ifig*PLOTPERFIG])]);
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
                axes(ax(i_comp))
                hold on
                for irej = 1:numel(rejects)
                    if isfield(EEG.reject.SASICA,rejfields{irej,1}) && ...
                            EEG.reject.SASICA.(rejfields{irej,1})(i_comp)
                        x = -.5 + (irej > 6);
                        y = .5 - .1*irej-.3*(rem(irej-1,6)+1>3);
                        scatter(x,y,'markerfacecolor',EEG.reject.SASICA.([rejfields{irej} 'col']),'markeredgecolor',EEG.reject.SASICA.([rejfields{irej} 'col']));
                    end
                end
            end
        end

        try
            pop_selectcomps(EEG, [ncomp+1]);
        catch
            hlastfig = gcf;
            set(hlastfig,'name',[get(hlastfig,'name') ' -- SASICA']);
            lastax = findobj(hlastfig,'type','Axes');
            lastax(end) = [];
            set(lastax,'visible','off');
        end
        axes(lastax);
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function setctxt(hfig,EEG,cfg)
COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
buttons = findobj(hfig,'-regexp','tag','^comp\d{1,3}$');
buttonnums = regexp(get(buttons,'tag'),'comp(\d{1,3})','tokens');
buttonnums = cellfun(@(x)(str2num(x{1}{1})),buttonnums);
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
    
    if cfg.ADJUST.enable
        % set ADJUST context menu

        cb = regexprep(get(buttons(i),'Callback'),'pop_prop','eeg_SASICA(EEG,''pop_prop_ADJ');
        ADJ = EEG.reject.SASICA.icaADJUST;
        str = num2str(any(ADJ.horiz == buttonnums(i)));
        str = [str ',' num2str(any(ADJ.vert == buttonnums(i)))];
        str = [str ',' num2str(any(ADJ.blink == buttonnums(i)))];
        str = [str ',' num2str(any(ADJ.disc == buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_DV)];
        str = [str ',' num2str(ADJ.diff_var(buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_K)];
        str = [str ',' num2str(ADJ.med2_K)];
        str = [str ',' num2str(ADJ.meanK(buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_SED)];
        str = [str ',' num2str(ADJ.med2_SED)];
        str = [str ',' num2str(ADJ.SED(buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_SAD)];
        str = [str ',' num2str(ADJ.med2_SAD)];
        str = [str ',' num2str(ADJ.SAD(buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_GDSF)];
        str = [str ',' num2str(ADJ.med2_GDSF)];
        str = [str ',' num2str(ADJ.GDSF(buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_V)];
        str = [str ',' num2str(ADJ.med2_V)];
        str = [str ',' num2str(ADJ.nuovaV(buttonnums(i)))];
        str = [str ',' num2str(ADJ.soglia_D)];
        str = [str ',' num2str(ADJ.maxdin(buttonnums(i)))];
        str = [str ');'''];
        cb = [regexprep(cb,'\{.*\}',str)];
        uimenu(hcmenu, 'Label', 'ADJUST results', 'Callback', cb,'tag',['ctxt_ADJ' num2str(buttonnums(i))]);
    end
    if cfg.FASTER.enable
        % set FASTER ctxt menu
        cb = regexprep(get(buttons(i),'Callback'),'pop_prop','eeg_SASICA(EEG,''pop_prop_FST');
        FST = EEG.reject.SASICA.icaFASTER;
        str = num2str(FST.listprops);
        str = [str repmat(';',size(str,1),1)]';
        str = ['['; str(:); ']'; ')' ;';' ;'''']';
        cb = regexprep(cb,'\{.*\}',str) ;
        uimenu(hcmenu, 'Label', 'FASTER results', 'Callback', cb,'tag',['ctxt_FST' num2str(buttonnums(i))]);
    end
    mycb = strrep(get(buttons(i),'Callback'),'''','''''');
    mycb = regexprep(mycb,'pop_prop','eeg_SASICA(EEG,''pop_prop');
    mycb = [mycb ''');'];
    set(buttons(i),'CallBack',mycb)
    set(buttons(i),'uicontextmenu',hcmenu)
end

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
p = panel();
p.margin = [10 10 10 10];
p.pack('v',{.35 []});
p(1).margintop = 10;
p(1).pack('h',{.4 []});

% plotting topoplot
p(1,1).select()
topoplot( EEG.icawinv(:,chanorcomp), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
    'shading', 'interp', 'numcontour', 3); axis square;
title(basename, 'fontsize', 14); 

% plotting erpimage
p(1,2).margin = [15 15 15 15];
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
if ~exist('winhandle')
    winhandle = NaN;
end;
p(2).pack('v',{.3 []})
p(2,1).margin = [25 25 10 10];
p(2,1).select();
try
    spectopo( EEG.icaact(chanorcomp,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
    % set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 12); 
	% set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 12); 
    axis on
    xlabel('')
	title('Activity power spectrum', 'fontsize', 10); 
    set(gca,'fontSize',10)
catch
	axis off;
    lasterror
	text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;

%%%% Add SASICA measures.

computed = fieldnames(EEG.reject.SASICA);
computed = computed(regexpcell(computed,'rej','inv'));
SASICAplot = [];FSTplot = []; ADJplot = [];
SASICAplotX = {};FSTis = []; ADJis = [];
for i = 1:numel(computed)
    if strcmp(computed{i},'icaADJUST')
        struct2ws(EEG.reject.SASICA.icaADJUST)
        ADJplot=[(SAD(chanorcomp)-med2_SAD)/(soglia_SAD-med2_SAD)
            (SED(chanorcomp)-med2_SED)/(soglia_SED-med2_SED)
            (GDSF(chanorcomp)-med2_GDSF)/(soglia_GDSF-med2_GDSF)
            (nuovaV(chanorcomp)-med2_V)/(soglia_V-med2_V)
            (meanK(chanorcomp)-med2_K)/(soglia_K-med2_K)]';
        ADJis = '';
        if ismember(chanorcomp,horiz)
            ADJis = [ADJis 'HEM/'];
        end
        if ismember(chanorcomp,vert)
            ADJis = [ADJis 'VEM/'];
        end
        if ismember(chanorcomp,vert)
            ADJis = [ADJis 'Blink/'];
        end
        if ismember(chanorcomp,vert)
            ADJis = [ADJis 'Disc/'];
        end
        if isempty(ADJis)
            ADJis = 'OK';
        else
            ADJis(end) = [];
        end
    elseif strcmp(computed{i},'icaFASTER')
        listprops = EEG.reject.SASICA.icaFASTER.listprops;
        str='FASTER - Detected as ';
        FASTER_reasons = {'High freq ' 'flat spectrum ' 'spatial Kurtosis ' 'Hurst exponent ' 'EOG correl '};
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
            str = 'FASTER says ok.';
        else
            str = [str reasons{:}];
        end
        FSTis = str;
        FSTplot = zlist(chanorcomp,:);
    else
        rejfields = {'icaautocorr' 'AutoCorr'
            'icafocalcomp' 'FocC'
            'icatrialfoc' 'FocT'
            'icaSNR' 'SNR'
            'icaresvar' 'ResV'
            'icachancorr' 'CorrC'
            'icaADJUST' 'ADJUST selections'
            'icaFASTER' 'FASTER selections'
            };
        SASICAplot(end+1) = EEG.reject.SASICA.(computed{i})(chanorcomp);
        SASICAplotX{end+1} = rejfields{strcmp(computed{i},rejfields(:,1)),2};
    end
end


h = axes('units','normalized', 'position',[40 35 65 20].*s+q);

h = axes('units','normalized', 'position',[40 5 65 20].*s+q);
h = axes('units','normalized', 'position',[-10 5 40 20].*s+q);


	
% display buttons
% ---------------
if ~isnan(winhandle)
	COLREJ = '[1 0.6 0.6]';
	COLACC = '[0.75 1 0.75]';
	% CANCEL button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

	% VALUE button
	% -------------
	hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

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

	% HELP button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop'');');

	% OK button
	% ---------
 	command = [ 'global EEG;' ...
 				'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
 				'EEG.reject.gcompreject(' num2str(chanorcomp) ') = tmpstatus;' ]; 
	if winhandle ~= 0
	 	command = [ command ...
	 				sprintf('if tmpstatus set(%3.15f, ''backgroundcolor'', %s); else set(%3.15f, ''backgroundcolor'', %s); end;', ...
					winhandle, COLREJ, winhandle, COLACC)];
	end;				
	command = [ command 'close(gcf); clear tmpstatus' ];
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

	% draw the figure for statistical values
	% --------------------------------------
	index = num2str( chanorcomp );
	command = [ ...
		'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
		'' ...
		'pos = get(gcf,''Position'');' ...
		'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
		'pos = get(gca,''position'');' ...
		'q = [pos(1) pos(2) 0 0];' ...
		's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
		'axis off;' ...
		''  ...
		'txt1 = sprintf(''(\n' ...
						'Entropy of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ' AND                 \t\t\t----\n\n' ...
					    'Kurtosis of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ') OR                  \t\t\t----\n\n' ...
					    'Kurtosis distibution \t\t\t%2.2f\n' ...
					    '> Rejection threhold\t\t\t%2.2f\n\n' ...
					    '\n' ...
					    'Current thesholds sujest to %s the component\n\n' ...
					    '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
						'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
						'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
		'' ...				
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
		'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
		'clear txt1 q s h pos;' ];
	set( hval, 'callback', command); 
	if isempty( EEG.stats.compenta )
		set(hval, 'enable', 'off');
	end;
	
	com = sprintf('pop_prop( %s, %d, %d, 0, %s);', inputname(1), typecomp, chanorcomp, vararg2str( { spec_opt } ) );
else
	com = sprintf('pop_prop( %s, %d, %d, NaN, %s);', inputname(1), typecomp, chanorcomp, vararg2str( { spec_opt } ) );
end;

return;

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


function [nb,channame,strnames] = chnb(channame, varargin)

% chnb() - return channel number corresponding to channel names in an EEG
%           structure
%
% Usage:
%   >> [nb]                 = chnb(channameornb);
%   >> [nb,names]           = chnb(channameornb,...);
%   >> [nb,names,strnames]  = chnb(channameornb,...);
%   >> [nb]                 = chnb(channameornb, labels);
%
% Input:
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
error(nargchk(1,2,nargin));
if nargin == 2
    labels = varargin{1};
else

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
            return
        end
    end
    if numel(channame) == 1 && not(isempty(strmatch('inv',channame{1})))
        cmd = 'exactinv';
        channame{1} = strrep(channame{1},'inv','');
    else
        channame{1} = channame{1};
        cmd = 'exact';
    end
    nb = regexpcell(labels,channame,[cmd 'ignorecase']);

elseif isnumeric(channame)
    nb = channame;
    if nb > numel(labels)
        nb = [];
    end
end
channame = labels(nb);
strnames = sprintf('%s ',channame{:});
if not(isempty(strnames))
    strnames(end) = [];
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

error(nargchk(2,3,nargin))
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

function struct2ws(s,varargin)

% struct2ws(s,varargin)
%
% Description : This function returns fields of scalar structure s in the
% current workspace
% __________________________________
% Inputs :
%   s (scalar structure array) :    a structure that you want to throw in
%                                   your current workspace.
%   re (string optional) :          a regular expression. Only fields
%                                   matching re will be returned
% Outputs :
%   No output : variables are thrown directly in the caller workspace.
%
%
% _____________________________________
% See also : ws2struct ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if nargin == 0
    cd('d:\Bureau\work')
    s = dir('pathdef.m');
end
if length(s) > 1
    error('Structure should be scalar.');
end
if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = fieldnames(s);
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i})
        continue
    end
    varsmatch(end+1) = i;
end
for i = varsmatch
    assignin('caller',vars{i},s.(vars{i}));
end

function [sortie] = ws2struct(varargin)

% [s] = ws2struct(varargin)
%
% Description : This function returns a structure containing variables
% of the current workspace.
% __________________________________
% Inputs :
%   re (string optional) :  a regular expression matching the variables to
%                           be returned.
% Outputs :
%   s (structure array) :   a structure containing all variables of the
%                           calling workspace. If re input is specified,
%                           only variables matching re are returned.
% _____________________________________
% See also : struct2ws ; regexp
%
% Maximilien Chaumon v1.0 02/2007


if not(isempty(varargin))
    re = varargin{1};
else
    re = '.*';
end

vars = evalin('caller','who');
vmatch = regexp(vars,re);
varsmatch = [];
for i = 1:length(vmatch)
    if isempty(vmatch{i}) || not(vmatch{i} == 1)
        continue
    end
    varsmatch{end+1} = vars{i};
end

for i = 1:length(varsmatch)
    dat{i} = evalin('caller',varsmatch{i});
end

sortie = cell2struct(dat,varsmatch,2);

function [tpts tvals] = timepts(timein, varargin)

% timepts() - return time points corresponding to a certain latency range
%             in an EEG structure.
%
% Usage:
%   >> [tpts] = timepts(timein);
%   >> [tpts tvals] = timepts(timein, times);
%               Note: this last method also works with any type of numeric
%               data entered under times (ex. frequencies, trials...)
%
% Input:
%   timein        - latency range [start stop] (boundaries included). If
%                   second argument 'times' is not provided, EEG.times will
%                   be evaluated from the EEG structure found in the caller
%                   workspace (or base if not in caller).
%   times         - time vector as found in EEG.times
%
% Output:
%   tpts          - index numbers corresponding to the time range.
%   tvals         - values of EEG.times at points tpts
%

error(nargchk(1,2,nargin));
if nargin == 2
    times = varargin{1};
else

    try
        EEG = evalin('caller','EEG');
    catch
        try
            EEG = evalin('base','EEG');
        catch
            error('Could not find EEG structure');
        end
    end
    if not(isfield(EEG,'times'))
        error('No time list found');
    end
    times = EEG.times;
    if isempty(times)
        times = EEG.xmin:1/EEG.srate:EEG.xmax;
    end
end
if isempty(times)
    error('could not find times');
end
if numel(timein) == 1
    [dum tpts] = min(abs(times - timein));% find the closest one
    if tpts == numel(times)
        warning('Strange time is last index of times')
    end
elseif numel(timein) == 2
    tpts = find(times >= timein(1) & times <= timein(2));% find times within bounds
else
    error('timein should be a scalar or a 2 elements vector');
end
tvals = times(tpts);


function tw = strwrap(t,n)

% tw = strwrap(t,n)
%
% wrap text array t at n characters taking non alphanumeric characters as
% breaking characters (i.e. not cutting words strangely).

t = deblank(t(:)');
seps = '\W';
breaks = regexp(t,seps);
breaks(end+1) = numel(t);
tw = '';
while not(isempty(t))
    idx = 1:min(n,breaks(find(breaks < n, 1,'last')));
    if isempty(idx)
        idx = 1:min(n,numel(t));
    end
    tw(end+1,:) = char(padarray(double(t(idx)),[0 n-numel(idx)],32,'post'));
    t(idx)= [];
    t = strtrim(t);
    breaks = breaks-numel(idx);
end


function [z,mu,sigma] = zscore(x,flag,dim)
%ZSCORE Standardized z score.
%   Z = ZSCORE(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-MEAN(X)) ./ STD(X). For
%   matrix X, z-scores are computed using the mean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the mean and standard deviation along the first
%   non-singleton dimension.
%
%   The columns of Z have sample mean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = ZSCORE(X) also returns MEAN(X) in MU and STD(X) in SIGMA.
%
%   [...] = ZSCORE(X,1) normalizes X using STD(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which ZSCORE works.  ZSCORE(X,0) is the same as
%   ZSCORE(X).
%
%   [...] = ZSCORE(X,FLAG,DIM) standardizes X by working along the dimension
%   DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1
%   to use N.
%
%   See also MEAN, STD.

%   Copyright 1993-2006 The MathWorks, Inc.
%   $Revision: 1.1.6.1 $  $Date: 2010/03/16 00:18:33 $

% [] is a special case for std and mean, just handle it out here.
if isequal(x,[]), z = []; return; end

if nargin < 2
    flag = 0;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim);
sigma = std(x,flag,dim);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);

function narginchk(min,max)

n = evalin('caller','nargin');
if  n < min || n > max
    error('number of arguments')
end
function h = vline(x,varargin)

% h = vline(x,varargin)
% add vertical line(s) on the current axes at x
% all varargin arguments are passed to plot...

x = x(:);
ho = ishold;
hold on
h = plot([x x]',repmat(ylim,numel(x),1)',varargin{:});
if not(ho)
    hold off
end
if nargout == 0
    clear h
end
function h = hline(y,varargin)

% h = hline(y,varargin)
% add horizontal line(s) on the current axes at y
% all varargin arguments are passed to plot...

y = y(:);
ho = ishold;
hold on
h = plot(repmat(xlim,numel(y),1)',[y y]',varargin{:});
if not(ho)
    hold off
end
if nargout == 0
    clear h
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% BELOW IS ADJUST CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ADJUST() - Automatic EEG artifact Detector 
% with Joint Use of Spatial and Temporal features
%
% Usage:
%   >> [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
%         soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
% 
% Inputs:
%   EEG        - current dataset structure or structure array (has to be epoched)
%   out        - (string) report file name 
%
% Outputs:
%   art        - List of artifacted ICs
%   horiz      - List of HEM ICs 
%   vert       - List of VEM ICs   
%   blink      - List of EB ICs     
%   disc       - List of GD ICs     
%   soglia_DV  - SVD threshold      
%   diff_var   - SVD feature values
%   soglia_K   - TK threshold      
%   meanK      - TK feature values
%   soglia_SED - SED threshold      
%   SED        - SED feature values
%   soglia_SAD - SAD threshold      
%   SAD        - SAD feature values
%   soglia_GDSF- GDSF threshold      
%   GDSF       - GDSF feature values
%   soglia_V   - MEV threshold      
%   nuovaV     - MEV feature values
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ADJUST
% Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features
% 
% Developed 2007-2014
% Andrea Mognon (1) and Marco Buiatti (2), 
% (1) Center for Mind/Brain Sciences, University of Trento, Italy
% (2) INSERM U992 - Cognitive Neuroimaging Unit, Gif sur Yvette, France
% 
% Last update: 02/05/2014 by Marco Buiatti
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Reference paper:
% Mognon A, Bruzzone L, Jovicich J, Buiatti M, 
% ADJUST: An Automatic EEG artifact Detector based on the Joint Use of Spatial and Temporal features. 
% Psychophysiology 48 (2), 229-240 (2011).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERSIONS LOG
%
% 02/05/14: Modified text in Report.txt (MB).
%
% 30/03/14: Removed 'message to the user' (redundant). (MB)
% 
% 22/03/14: kurtosis is replaced by kurt for compatibility if signal processing
%           toolbox is missing (MB).
%
% V2 (07 OCTOBER 2010) - by Andrea Mognon
% Added input 'nchannels' to compute_SAD and compute_SED_NOnorm;
% this is useful to differentiate the number of ICs (n) and the number of
% sensors (nchannels);
% bug reported by Guido Hesselman on October, 1 2010.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% function [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
%         soglia_GDSF, GDSF, soglia_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out)
function [art, horiz, vert, blink, disc,...
        soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG)

    
%% Settings

% ----------------------------------------------------
% |  Change experimental settings in this section    |
% ----------------------------------------------------

% ----------------------------------------------------
% |  Initial message to user:                        |
% ----------------------------------------------------
% 
% disp(' ')
% disp('Detects Horizontal and Vertical eye movements,')
% disp('Blinks and Discontinuities in dataset:')
% disp([EEG.filename])
% disp(' ')

% ----------------------------------------------------
% |  Collect useful data from EEG structure          |
% ----------------------------------------------------

%number of ICs=size(EEG.icawinv,1);

%number of time points=size(EEG.data,2);

if length(size(EEG.data))==3
    
    num_epoch=size(EEG.data,3);
       
else
    
    num_epoch=0;

end

% Check the presence of ICA activations
EEG.icaact = eeg_getica(EEG);

topografie=EEG.icawinv'; %computes IC topographies

% Topographies and time courses normalization
% 
% disp(' ');
% disp('Normalizing topographies...')
% disp('Scaling time courses...')

for i=1:size(EEG.icawinv,2) % number of ICs
    
    ScalingFactor=norm(topografie(i,:));
    
    topografie(i,:)=topografie(i,:)/ScalingFactor;
 
    if length(size(EEG.data))==3
        EEG.icaact(i,:,:)=ScalingFactor*EEG.icaact(i,:,:);
    else
        EEG.icaact(i,:)=ScalingFactor*EEG.icaact(i,:);
    end
    
end
% 
% disp('Done.')
% disp(' ')

% Variables memorizing artifacted ICs indexes

blink=[];

horiz=[];

vert=[];

disc=[];

%% Check EEG channel position information
nopos_channels=[];
for el=1:length(EEG.chanlocs)
    if(any(isempty(EEG.chanlocs(1,el).X)&isempty(EEG.chanlocs(1,el).Y)&isempty(EEG.chanlocs(1,el).Z)&isempty(EEG.chanlocs(1,el).theta)&isempty(EEG.chanlocs(1,el).radius)))
        nopos_channels=[nopos_channels el];
    end;
end

if ~isempty(nopos_channels)
    disp(['Warning : Channels ' num2str(nopos_channels) ' have incomplete location information. They will NOT be used to compute ADJUST spatial features']);
    disp(' ');
end;

pos_channels=setdiff(1:length(EEG.chanlocs),nopos_channels);

%% Feature extraction

disp(' ')
disp('Features Extraction:')

%GDSF - General Discontinuity Spatial Feature

disp('GDSF - General Discontinuity Spatial Feature...')

GDSF = compute_GD_feat(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SED - Spatial Eye Difference

disp('SED - Spatial Eye Difference...')

[SED,medie_left,medie_right]=computeSED_NOnorm(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2)); 


%SAD - Spatial Average Difference

disp('SAD - Spatial Average Difference...')

[SAD,var_front,var_back,mean_front,mean_back]=computeSAD(topografie,EEG.chanlocs(1,pos_channels),size(EEG.icawinv,2));


%SVD - Spatial Variance Difference between front zone and back zone

diff_var=var_front-var_back;

%epoch dynamic range, variance and kurtosis

K=zeros(num_epoch,size(EEG.icawinv,2)); %kurtosis
Kloc=K;

Vmax=zeros(num_epoch,size(EEG.icawinv,2)); %variance

% disp('Computing variance and kurtosis of all epochs...')

for i=1:size(EEG.icawinv,2) % number of ICs
    
    for j=1:num_epoch              
        Vmax(j,i)=var(EEG.icaact(i,:,j));        
%         Kloc(j,i)=kurtosis(EEG.icaact(i,:,j));
        K(j,i)=kurt(EEG.icaact(i,:,j));
    end  
end

% check that kurt and kurtosis give the same values:
% [a,b]=max(abs(Kloc(:)-K(:)))

%TK - Temporal Kurtosis

disp('Temporal Kurtosis...')

meanK=zeros(1,size(EEG.icawinv,2));

for i=1:size(EEG.icawinv,2)
    if num_epoch>100
    meanK(1,i)=trim_and_mean(K(:,i)); 
    else meanK(1,i)=mean(K(:,i));
    end

end


%MEV - Maximum Epoch Variance

disp('Maximum epoch variance...')

maxvar=zeros(1,size(EEG.icawinv,2));
meanvar=zeros(1,size(EEG.icawinv,2));


for i=1:size(EEG.icawinv,2)
    if num_epoch>100
     maxvar(1,i)=trim_and_max(Vmax(:,i)');
     meanvar(1,i)=trim_and_mean(Vmax(:,i)');
    else 
     maxvar(1,i)=max(Vmax(:,i));
     meanvar(1,i)=mean(Vmax(:,i));
    end
end

% MEV in reviewed formulation:

nuovaV=maxvar./meanvar;



%% Thresholds computation

disp('Computing EM thresholds...')

% soglia_K=EM(meanK);
% 
% soglia_SED=EM(SED);
% 
% soglia_SAD=EM(SAD);
% 
% soglia_GDSF=EM(GDSF);
% 
% soglia_V=EM(nuovaV); 
[soglia_K,med1_K,med2_K]=EM(meanK);

[soglia_SED,med1_SED,med2_SED]=EM(SED);

[soglia_SAD,med1_SAD,med2_SAD]=EM(SAD);

[soglia_GDSF,med1_GDSF,med2_GDSF]=EM(GDSF);

[soglia_V,med1_V,med2_V]=EM(nuovaV); 



%% Horizontal eye movements (HEM)


horiz=intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
    (find(nuovaV>=soglia_V)));




%% Vertical eye movements (VEM)



vert=intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
    intersect(find(diff_var>0),find(nuovaV>=soglia_V)));



%% Eye Blink (EB)


blink=intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
    intersect ( find(meanK>=soglia_K),find(diff_var>0) ));



%% Generic Discontinuities (GD)



disc=intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));


%compute output variable
art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs

% these three are old outputs which are no more necessary in latest ADJUST version.
soglia_D=0;
soglia_DV=0;
maxdin=zeros(1,size(EEG.icawinv,2));

return

% compute_GD_feat() - Computes Generic Discontinuity spatial feature 
%
% Usage:
%   >> res = compute_GD_feat(topografie,canali,num_componenti);
%
% Inputs:
%   topografie - topographies vector
%   canali     - EEG.chanlocs struct
%   num_componenti  - number of components
%
% Outputs:
%   res       - GDSF values

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


function res = compute_GD_feat(topografie,canali,num_componenti)

% Computes GDSF, discontinuity spatial feature
% topografie is the topography weights matrix
% canali is the structure EEG.chanlocs
% num_componenti is the number of ICs
% res is GDSF values

xpos=[canali.X];ypos=[canali.Y];zpos=[canali.Z];
pos=[xpos',ypos',zpos'];

res=zeros(1,num_componenti);   

for ic=1:num_componenti
    
    % consider the vector topografie(ic,:)
    
    aux=[];
    
    for el=1:length(canali)-1
        
        P=pos(el,:); %position of current electrode
        d=pos-repmat(P,length(canali),1);
        %d=pos-repmat(P,62,1);
        dist=sqrt(sum((d.*d),2));
        
        [y,I]=sort(dist);
        repchas=I(2:11); % list of 10 nearest channels to el
        weightchas=exp(-y(2:11)); % respective weights, computed wrt distance
        
        aux=[aux abs(topografie(ic,el)-mean(weightchas.*topografie(ic,repchas)'))];
                % difference between el and the average of 10 neighbors
                % weighted according to weightchas
    end
    
    res(ic)=max(aux);
    
end


% computeSAD() - Computes Spatial Average Difference feature 
%
% Usage:
%   >> [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n);
%
% Inputs:
%   topog      - topographies vector
%   chanlocs   - EEG.chanlocs struct
%   n          - number of ICs
%   nchannels  - number of channels
%
% Outputs:
%   rapp       - SAD values
%   var_front  - Frontal Area variance values
%   var_back   - Posterior Area variance values
%   mean_front - Frontal Area average values
%   mean_back  - Posterior Area average values
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


function [rapp,var_front,var_back,mean_front,mean_back]=computeSAD(topog,chanlocs,n)

nchannels=length(chanlocs);

%% Define scalp zones

% Find electrodes in Frontal Area (FA)
dimfront=0; %number of FA electrodes
index1=zeros(1,nchannels); %indexes of FA electrodes

for k=1:nchannels
    if (abs(chanlocs(1,k).theta)<60) && (chanlocs(1,k).radius>0.40) %electrodes are in FA
        dimfront=dimfront+1; %count electrodes
        index1(1,dimfront)=k; 
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
 
    if dimfront*dimback==0
        disp('ERROR: no channels included in some scalp areas.')
        disp('Check channels distribution and/or change scalp areas definitions in computeSAD.m and computeSED_NOnorm.m')
        disp('ADJUST session aborted.')
        return
    end
    
%% Outputs

     rapp=zeros(1,n); % SAD
      mean_front=zeros(1,n); % FA electrodes mean value
      mean_back=zeros(1,n); % PA electrodes mean value
      var_front=zeros(1,n); % FA electrodes variance value
      var_back=zeros(1,n); % PA electrodes variance value

%% Output computation

for i=1:n % for each topography
    
 %create FA electrodes vector
    front=zeros(1,dimfront);
    for h=1:dimfront
        front(1,h)=topog(i,index1(1,h));
    end
    
  %create PA electrodes vector
    back=zeros(1,dimback);
    for h=1:dimback
        back(1,h)=topog(i,index3(1,h));
    end
    
     
    
   %compute features
    
    rapp(1,i)=abs(mean(front))-abs(mean(back)); % SAD
    mean_front(1,i)=mean(front);
    mean_back(1,i)=mean(back);
    var_back(1,i)=var(back);
    var_front(1,i)=var(front);
    
end


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

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EM - ADJUST package
% 
% Performs automatic threshold on the digital numbers 
% of the input vector 'vec'; based on Expectation - Maximization algorithm

% Reference paper:
% Bruzzone, L., Prieto, D.F., 2000. Automatic analysis of the difference image 
% for unsupervised change detection. 
% IEEE Trans. Geosci. Remote Sensing 38, 1171:1182

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Usage:
%   >> [last,med1,med2,var1,var2,prior1,prior2]=EM(vec);
%
% Input: vec (row vector, to be thresholded)
%
% Outputs: last (threshold value)
%          med1,med2 (mean values of the Gaussian-distributed classes 1,2)
%          var1,var2 (variance of the Gaussian-distributed classes 1,2)
%          prior1,prior2 (prior probabilities of the Gaussian-distributed classes 1,2)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


function [last,med1,med2,var1,var2,prior1,prior2]=EM(vec)

if size(vec,2)>1
	len=size(vec,2); %number of elements
else
	vec=vec';
	len=size(vec,2); 
end

c_FA=1; % False Alarm cost
c_MA=1; % Missed Alarm cost

med=mean(vec);
standard=std(vec);
mediana=(max(vec)+min(vec))/2;

alpha1=0.01*(max(vec)-mediana); % initialization parameter/ righthand side
alpha2=0.01*(mediana-min(vec)); % initialization parameter/ lefthand side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPECTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

train1=[]; % Expectation of class 1
train2=[];
train=[]; % Expectation of 'unlabeled' samples

for i=1:(len)
    if (vec(i)<(mediana-alpha2)) 
        train2=[train2 vec(i)];
    elseif (vec(i)>(mediana+alpha1))
        train1=[train1 vec(i)];
    else
	  train=[train vec(i)];
    end
end

n1=length(train1);
n2=length(train2);

med1=mean(train1);
med2=mean(train2);
prior1=n1/(n1+n2);
prior2=n2/(n1+n2);
var1=var(train1);
var2=var(train2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAXIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count=0;
dif_med_1=1; % difference between current and previous mean
dif_med_2=1;
dif_var_1=1; % difference between current and previous variance
dif_var_2=1;
dif_prior_1=1; % difference between current and previous prior
dif_prior_2=1;
stop=0.0001;

while((dif_med_1>stop)&&(dif_med_2>stop)&&(dif_var_1>stop)&&(dif_var_2>stop)&&(dif_prior_1>stop)&&(dif_prior_2>stop))

    count=count+1;

    med1_old=med1;
    med2_old=med2;
    var1_old=var1;
    var2_old=var2;
    prior1_old=prior1;
    prior2_old=prior2;
	prior1_i=[];
	prior2_i=[];

    % FOLLOWING FORMULATION IS ACCORDING TO REFERENCE PAPER:
    
    for i=1:len
		prior1_i=[prior1_i prior1_old*Bayes(med1_old,var1_old,vec(i))/...
		   (prior1_old*Bayes(med1_old,var1_old,vec(i))+prior2_old*Bayes(med2_old,var2_old,vec(i)))];
		prior2_i=[prior2_i prior2_old*Bayes(med2_old,var2_old,vec(i))/...
		   (prior1_old*Bayes(med1_old,var1_old,vec(i))+prior2_old*Bayes(med2_old,var2_old,vec(i)))];
    end
	
	
	prior1=sum(prior1_i)/len;
	prior2=sum(prior2_i)/len;
	med1=sum(prior1_i.*vec)/(prior1*len);
	med2=sum(prior2_i.*vec)/(prior2*len);
	var1=sum(prior1_i.*((vec-med1_old).^2))/(prior1*len);
	var2=sum(prior2_i.*((vec-med2_old).^2))/(prior2*len);

    dif_med_1=abs(med1-med1_old);
    dif_med_2=abs(med2-med2_old);
    dif_var_1=abs(var1-var1_old);
    dif_var_2=abs(var2-var2_old);
    dif_prior_1=abs(prior1-prior1_old);
    dif_prior_2=abs(prior2-prior2_old);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLDING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k=c_MA/c_FA;
a=(var1-var2)/2;
b= ((var2*med1)-(var1*med2));
c=(log((k*prior1*sqrt(var2))/(prior2*sqrt(var1)))*(var2*var1))+(((((med2)^2)*var1)-(((med1)^2)*var2))/2);
rad=(b^2)-(4*a*c);
if rad<0
    disp('Negative Discriminant!');
    return;
end

soglia1=(-b+sqrt(rad))/(2*a);
soglia2=(-b-sqrt(rad))/(2*a);

if ((soglia1<med2)||(soglia1>med1))
    last=soglia2;
else
    last=soglia1;
end

if isnan(last) % TO PREVENT CRASHES
    last=mediana;
end

return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prob=Bayes(med,var,point)
if var==0
    prob=1;
else
    prob=((1/(sqrt(2*pi*var)))*exp((-1)*((point-med)^2)/(2*var)));
end

% pop_prop_ADJ() - overloaded pop_prop() for ADJUST plugin.
%              plot the properties of a channel or of an independent component. 
%              ADJUST feature values are also shown (normalized wrt threshold).
%              
%
% Usage:
%   >> com = pop_prop_ADJ(EEG, typecomp, numcompo, winhandle, is_H, is_V, is_B, is_D,...
%     soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
%         soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin)
% 
% Inputs:
%   EEG        - current dataset structure or structure array 
%   typecomp   - [0|1] compute electrode property (1) or component 
%                property (0). Default is 1.
%   numcompo   - channel or component number
%   winhandle  - if this parameter is present or non-NaN, buttons for the
%                rejection of the component are drawn. If 
%                non-zero, this parameter is used to backpropagate
%                the color of the rejection button.
%   is_H       - (bool) true if plotted IC is HEM
%   is_V       - (bool) true if plotted IC is VEM
%   is_B       - (bool) true if plotted IC is EB
%   is_D       - (bool) true if plotted IC is GD
%   soglia_DV  - feature1 (SVD) threshold
%   diff_var   - feature1 (SVD) vector
%   soglia_K   - feature2 (TK) threshold
%   meanK      - feature2 (TK) vector 
%   soglia_SED - feature3 (SED) threshold
%   SED        - feature3 (SED) vector
%   soglia_SAD - feature4 (SAD) threshold
%   SAD        - feature4 (SAD) vector 
%   soglia_TDR - feature5 (SDR) threshold
%   topog_DR   - feature5 (SDR) vector
%   soglia_V   - feature6 (MEV) threshold
%   maxvar     - feature6 (MEV) vector 
%   soglia_D   - feature7 (MEDR) threshold
%   maxdin     - feature7 (MEDR) vector 
%

% ORIGINAL HELP:
% pop_prop() - plot the properties of a channel or of an independent
%              component. 
% Usage:
%   >> pop_prop( EEG, typecomp); % pops up a query window 
%   >> pop_prop( EEG, typecomp, chan, winhandle);
%
% Inputs:
%   EEG        - dataset structure (see EEGGLOBAL)
%   typecomp   - [0|1] compute electrode property (1) or component 
%                property (0). Default is 1.
%   chan       - channel or component number
%   winhandle  - if this parameter is present or non-NaN, buttons for the
%                rejection of the component are drawn. If 
%                non-zero, this parameter is used to backpropagate
%                the color of the rejection button.
%   spectral_options - [cell array] cell arry of options for the spectopo()
%                      function.
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
% 
%
%
% VERSIONS LOG
% 
% 26/03/14:
% 

function com = pop_prop_ADJ(EEG, typecomp, numcompo, winhandle, is_H, is_V, is_B, is_D,...
    soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin) %,spec_opt)

com = '';
if nargin < 1
	help pop_prop_ADJ;
	return;   
end;

if nargin == 1
	typecomp = 1;
end;
if typecomp == 0 & isempty(EEG.icaweights)
   error('No ICA weights recorded for this set, first run ICA');
end;   
if nargin == 2
	promptstr    = { fastif(typecomp,'Channel number to plot:','Component number to plot:') ...
                     'Spectral options (see spectopo help):' };
	inistr       = { '1' '''freqrange'', [2 50]' };
	result       = inputdlg2( promptstr, 'Component properties - pop_prop_ADJ()', 1,  inistr, 'pop_prop_ADJ');
	if size( result, 1 ) == 0 return; end;
   
	numcompo   = eval( [ '[' result{1} ']' ] );
    spec_opt   = eval( [ '{' result{2} '}' ] );
end;

% plotting several component properties - STILL TO CHANGE
% -------------------------------------
if length(numcompo) > 1
    for index = numcompo
        pop_prop(EEG, typecomp, index);
    end;
	com = sprintf('pop_prop( %s, %d, [%s]);', inputname(1), typecomp, int2str(numcompo));
    return;
end;

if numcompo < 1 | numcompo > EEG.nbchan
   error('Component index out of range');
end;   

% assumed input is numcompo
% -------------------------
try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;
basename = [fastif(typecomp,'Channel ', 'Component ') int2str(numcompo) ];

fh = figure('name', ['pop_prop_ADJ() - ' basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'off');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2)-500+pos(4) 500 500], 'visible', 'on');
pos = get(gca,'position'); % plot relative to current axes
hh = gca;
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
axis off;

% plotting topoplot
% -----------------
h = axes('Units','Normalized', 'Position',[-10 65 40 35].*s+q); 

%topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs); axis square; 

if typecomp == 1 % plot single channel locations
	topoplot( numcompo, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
             'electrodes','off', 'style', 'blank', 'emarkersize1chan', 12); axis square;
else             % plot component map
	topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
             'shading', 'interp', 'numcontour', 3); axis square;
end;
basename = [fastif(typecomp,'Channel ', 'IC') int2str(numcompo) ];
% title([ basename fastif(typecomp, ' location', ' map')], 'fontsize', 14); 
title(basename, 'fontsize', 12); 

% plotting erpimage
% -----------------
hhh = axes('Units','Normalized', 'Position',[45 67 48 33].*s+q); %era height 38
eeglab_options; 
if EEG.trials > 1
    % put title at top of erpimage
    axis off
    hh = axes('Units','Normalized', 'Position',[45 67 48 33].*s+q);
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
      ei_smooth = 1;
    else
      ei_smooth = 3;
    end
    if typecomp == 1 % plot component
         offset = nan_mean(EEG.data(numcompo,:));
         erpimage( EEG.data(numcompo,:)-offset, ones(1,EEG.trials)*10000, EEG.times , ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp');   
    else % plot channel
          if option_computeica  
                  offset = nan_mean(EEG.icaact(numcompo,:));
                  erpimage( EEG.icaact(numcompo,:)-offset, ones(1,EEG.trials)*10000, EEG.times , ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '');   
          else
                
              icaacttmp = (EEG.icaweights(numcompo,:) * EEG.icasphere) ...
                                   * reshape(EEG.data(1:size(EEG.icaweights,1),:,:), EEG.nbchan, EEG.trials*EEG.pnts);
%                     icaacttmp = (EEG.icaweights(numcompo,:) * EEG.icasphere) ...
%                                    * EEG.data(EEG.nbchan, EEG.trials*EEG.pnts);
                  offset = nan_mean(icaacttmp);
                  erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times, ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '');   
          end;
    end;
    axes(hhh);
    title(sprintf('%s activity \\fontsize{10}(global offset %3.3f)', basename, offset), 'fontsize', 12);
else

    % put title at top of erpimage
    EI_TITLE = 'Continous data';
    axis off
    hh = axes('Units','Normalized', 'Position',[45 62 48 38].*s+q);
    ERPIMAGELINES = 200; % show 200-line erpimage
    while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
       ERPIMAGELINES = round(0.9 * ERPIMAGELINES);
    end
    if ERPIMAGELINES > 2   % give up if data too small
      if ERPIMAGELINES < 10
         ei_smooth == 1;
      else
        ei_smooth = 3;
      end
      erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
      erpimageframestot = erpimageframes*ERPIMAGELINES;
      eegtimes = linspace(0, erpimageframes-1, EEG.srate/1000);
      if typecomp == 1 % plot component
           offset = nan_mean(EEG.data(numcompo,:));
           erpimage( reshape(EEG.data(numcompo,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar');   
      else % plot channel
            if option_computeica  
                    offset = nan_mean(EEG.icaact(numcompo,:));
                    erpimage( ...
              reshape(EEG.icaact(numcompo,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ...
                   ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel', '');   
            else
%                     icaacttmp = reshape(EEG.icaweights(numcompo,:) * EEG.icasphere) ...
%                                      * reshape(EEG.data, erpimageframes, ERPIMAGELINES);
                        
                        icaacttmp = EEG.icaweights(numcompo,:) * EEG.icasphere ...
                                     *EEG.data(:,1:erpimageframes*ERPIMAGELINES);

                    offset = nan_mean(icaacttmp);
                    erpimage( icaacttmp-offset, ones(1,ERPIMAGELINES)*10000, eegtimes, ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar', 'yerplabel', '');   
            end;
      end
    else
            axis off;
            text(0.1, 0.3, [ 'No erpimage plotted' 10 'for small continuous data']);
    end;
    axes(hhh);
end;	

% plotting spectrum
% -----------------
if ~exist('winhandle')
    winhandle = NaN;
end;
if ~isnan(winhandle)
	h = axes('units','normalized', 'position',[10 25 85 25].*s+q);
    %h = axes('units','normalized', 'position',[5 10 95 35].*s+q); %%%
    %CHANGE!
else
	h = axes('units','normalized', 'position',[10 15 85 30].*s+q);
    %h = axes('units','normalized', 'position',[5 0 95 40].*s+q); %%%
    %CHANGE!
end;
%h = axes('units','normalized', 'position',[45 5 60 40].*s+q);
try
	eeglab_options;
    %next instr added for correct function! Andrea
    option_computeica=1;
	if typecomp == 1
		%[spectra freqs] = spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate, spec_opt{:},'freqrange', [0 45] );
        [spectra freqs] = spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate, 'freqrange', [0 45] );
	else 
		if option_computeica  
            
			%[spectra freqs] = spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), spec_opt{:}, 'freqrange', [0 45]);
            % CONTROL ADDED FOR CONTINUOUS DATA
            if size(EEG.data,3)==1 
                EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
            end
            [spectra freqs] = spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo),  'freqrange', [0 45]);
		else
			if exist('icaacttmp')~=1, 
                
				icaacttmp = (EEG.icaweights(numcompo,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts); 
			end;
            
			%[spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), spec_opt{:} ,'freqrange', [0 45]);
            [spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), 'freqrange', [0 45]);
		end;
	end;
    set(gca,'fontsize',8);
    % set up new limits
    % -----------------
    %freqslim = 50;
	%set(gca, 'xlim', [0 min(freqslim, EEG.srate/2)]);
    %spectra = spectra(find(freqs <= freqslim));
	%set(gca, 'ylim', [min(spectra) max(spectra)]);
    
	%tmpy = get(gca, 'ylim');
    %set(gca, 'ylim', [max(tmpy(1),-1) tmpy(2)]);
	set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 8); 
	set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 8); 
	title('Activity power spectrum', 'fontsize', 12); 
catch
	axis off;
	text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;	

% ----------------------------------------------------------------
% plotting IC properties
% -----------------
if ~exist('winhandle')
    winhandle = NaN;
end;
if ~isnan(winhandle)
	h = axes('units','normalized', 'position',[3 2 95 10].*s+q);
else
	h = axes('units','normalized', 'position',[3 0 95 10].*s+q);
end;

axis off
str='ADJUST - Detected as ';
if is_H
    str=[str 'HEM '];
end
if is_V
    str=[str 'VEM '];
end
if is_B
    str=[str 'EB '];
end
if is_D
    str=[str 'GD'];
end

if (is_H || is_V || is_B || is_D)==0
    str='ADJUST - Not detected';
end
% text(0,0,[str 10 'TK ' num2str(meanK) '(' num2str(soglia_K) '); SAD ' num2str(SAD) '(' num2str(soglia_SAD) ...
%     '); SVD ' num2str(diff_var) '(' num2str(soglia_DV) '); SED ' num2str(SED) '(' num2str(soglia_SED) ...
%     ')' 10 'MEDR ' num2str(maxdin) '(' num2str(soglia_D) '); MEV ' num2str(maxvar) '(' num2str(soglia_V) ...
%     '); SDR ' num2str(topog_DR) '(' num2str(soglia_TDR) ')'],'FontSize',8);

% compute bar graph entries
% E=[ SAD/soglia_SAD SED/soglia_SED GDSF/soglia_GDSF maxvar/soglia_V meanK/soglia_K];
E=[(SAD-med2_SAD)/(soglia_SAD-med2_SAD) (SED-med2_SED)/(soglia_SED-med2_SED) (GDSF-med2_GDSF)/(soglia_GDSF-med2_GDSF) (maxvar-med2_V)/(soglia_V-med2_V) (meanK-med2_K)/(soglia_K-med2_K)];

% set bar colors
C={[1 0 0],[.6 0 .2],[1 1 0],[0 1 0], [0 1 1]};
% horizontal line
l=ones([1, length(E)+2]);
% plot
plot(0:length(E)+1 , l , 'Linewidth',2,'Color','k');
hold on
for i=1:length(E)
    v=zeros(1,length(E));
    v(i)=E(i);
    bar(v,'facecolor',C{i});
    title(str);
    set(gca,'XTickLabel',{'';'SAD';'SED';'GDSF';'MEV';'TK';''},'YTickLabel',{'0';'Threshold';'2*Threshold'},'YLim',[0 2])
end




% -----------------------------------------------------------------


% display buttons
% ---------------

if ~isnan(winhandle)
	COLREJ = '[1 0.6 0.6]';
	COLACC = '[0.75 1 0.75]';
	% CANCEL button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

	% VALUE button
	% -------------
	hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

	% REJECT button
	% -------------
	status = EEG.reject.gcompreject(numcompo);
	hr = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
				'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
	command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
				'if get(gcbo, ''userdata''),' ...
				'     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
				'else ' ...
				'     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
				'end;' ];					
	set( hr, 'callback', command); 

	% HELP button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop_ADJ'');');

	% OK button
	% ---------
 	command = [ 'global EEG;' ...
 				'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
 				'EEG.reject.gcompreject(' num2str(numcompo) ') = tmpstatus;' ]; 
	if winhandle ~= 0
	 	command = [ command ...
	 				sprintf('if tmpstatus set(%3.15f, ''backgroundcolor'', %s); else set(%3.15f, ''backgroundcolor'', %s); end;', ...
					winhandle, COLREJ, winhandle, COLACC)];
	end;				
	command = [ command 'close(gcf); clear tmpstatus' ];
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

	% draw the figure for statistical values
	% --------------------------------------
	index = num2str( numcompo );
	command = [ ...
		'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
		'' ...
		'pos = get(gcf,''Position'');' ...
		'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
		'pos = get(gca,''position'');' ...
		'q = [pos(1) pos(2) 0 0];' ...
		's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
		'axis off;' ...
		''  ...
		'txt1 = sprintf(''(\n' ...
						'Entropy of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ' AND                 \t\t\t----\n\n' ...
					    'Kurtosis of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ') OR                  \t\t\t----\n\n' ...
					    'Kurtosis distibution \t\t\t%2.2f\n' ...
					    '> Rejection threhold\t\t\t%2.2f\n\n' ...
					    '\n' ...
					    'Current thesholds sujest to %s the component\n\n' ...
					    '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
						'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
						'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
		'' ...				
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
		'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
		'clear txt1 q s h pos;' ];
	set( hval, 'callback', command); 
	if isempty( EEG.stats.compenta )
		set(hval, 'enable', 'off');
	end;
	
    % MODIFICA
	%com = sprintf('pop_prop( %s, %d, %d, 0, %s);', inputname(1), typecomp, numcompo, vararg2str( { spec_opt } ) );
    	com = sprintf('pop_prop_ADJ( %s, %d, %d, 0, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s);',...
            inputname(1), typecomp, numcompo, is_H, is_V, is_B, is_D,...
            soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin );

else
	%com = sprintf('pop_prop( %s, %d, %d, NaN, %s);', inputname(1), typecomp, numcompo, vararg2str( { spec_opt } ) );
    	com = sprintf('pop_prop_ADJ( %s, %d, %d, NaN, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %s);',...
            inputname(1), typecomp, numcompo, is_H, is_V, is_B, is_D,...
            soglia_DV, diff_var, soglia_K, med2_K, meanK, soglia_SED, med2_SED, SED, soglia_SAD, med2_SAD, SAD, ...
        soglia_GDSF, med2_GDSF, GDSF, soglia_V, med2_V, maxvar, soglia_D, maxdin );

end;

function com = pop_prop_FST(EEG, typecomp, numcompo, winhandle,listprops) %,spec_opt)

% pop_prop_FST() - overloaded pop_prop() for FASTER plugin in autorejICA.
%              
%
% Usage:
%   >> com = pop_prop_FST(EEG, typecomp, numcompo, winhandle,listprops);
%
% Inputs:
%   EEG        - current dataset structure or structure array 
%   typecomp   - [0|1] compute electrode property (1) or component 
%                property (0). Default is 1.
%   numcompo   - channel or component number
%   winhandle  - if this parameter is present or non-NaN, buttons for the
%                rejection of the component are drawn. If 
%                non-zero, this parameter is used to backpropagate
%                the color of the rejection button.
%   listprops   - property listing as returned by FASTER
%               a ncomps x 5 matrix with 
%                     1 Median gradient value, for high frequency stuff
%                     2 Mean slope around the LPF band (spectral)
%                     3 Kurtosis of spatial map
%                     4 Hurst exponent
%                     5 Eyeblink correlations
%

% ORIGINAL HELP:
% pop_prop() - plot the properties of a channel or of an independent
%              component. 
% Usage:
%   >> pop_prop( EEG, typecomp); % pops up a query window 
%   >> pop_prop( EEG, typecomp, chan, winhandle);
%
% Inputs:
%   EEG        - dataset structure (see EEGGLOBAL)
%   typecomp   - [0|1] compute electrode property (1) or component 
%                property (0). Default is 1.
%   chan       - channel or component number
%   winhandle  - if this parameter is present or non-NaN, buttons for the
%                rejection of the component are drawn. If 
%                non-zero, this parameter is used to backpropagate
%                the color of the rejection button.
%   spectral_options - [cell array] cell arry of options for the spectopo()
%                      function.
% 

% Copyright (C) 2009 Andrea Mognon and Marco Buiatti, 
% Center for Mind/Brain Sciences, University of Trento, Italy
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


com = '';
if nargin < 1
	help pop_prop_FST;
	return;   
end;

if nargin == 1
	typecomp = 1;
end;
if typecomp == 0 & isempty(EEG.icaweights)
   error('No ICA weights recorded for this set, first run ICA');
end;   
if nargin == 2
	promptstr    = { fastif(typecomp,'Channel number to plot:','Component number to plot:') ...
                     'Spectral options (see spectopo help):' };
	inistr       = { '1' '''freqrange'', [2 50]' };
	result       = inputdlg2( promptstr, 'Component properties - pop_prop_FST()', 1,  inistr, 'pop_prop_FST');
	if size( result, 1 ) == 0 return; end;
   
	numcompo   = eval( [ '[' result{1} ']' ] );
    spec_opt   = eval( [ '{' result{2} '}' ] );
end;

% plotting several component properties - STILL TO CHANGE
% -------------------------------------
if length(numcompo) > 1
    for index = numcompo
        pop_prop(EEG, typecomp, index);
    end;
	com = sprintf('pop_prop( %s, %d, [%s]);', inputname(1), typecomp, int2str(numcompo));
    return;
end;

if numcompo < 1 | numcompo > EEG.nbchan
   error('Component index out of range');
end;   

% assumed input is numcompo
% -------------------------
try, icadefs; 
catch, 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end;
basename = [fastif(typecomp,'Channel ', 'Component ') int2str(numcompo) ];

fh = figure('name', ['pop_prop_FST() - ' basename ' properties'], 'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'off');
pos = get(gcf,'Position');
set(gcf,'Position', [pos(1) pos(2)-500+pos(4) 500 500], 'visible', 'on');
pos = get(gca,'position'); % plot relative to current axes
hh = gca;
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
axis off;

% plotting topoplot
% -----------------
h = axes('Units','Normalized', 'Position',[-10 65 40 35].*s+q); 

%topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs); axis square; 

if typecomp == 1 % plot single channel locations
	topoplot( numcompo, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
             'electrodes','off', 'style', 'blank', 'emarkersize1chan', 12); axis square;
else             % plot component map
	topoplot( EEG.icawinv(:,numcompo), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
             'shading', 'interp', 'numcontour', 3); axis square;
end;
basename = [fastif(typecomp,'Channel ', 'IC') int2str(numcompo) ];
% title([ basename fastif(typecomp, ' location', ' map')], 'fontsize', 14); 
title(basename, 'fontsize', 12); 

% plotting erpimage
% -----------------
hhh = axes('Units','Normalized', 'Position',[45 67 48 33].*s+q); %era height 38
eeglab_options; 
if EEG.trials > 1
    % put title at top of erpimage
    axis off
    hh = axes('Units','Normalized', 'Position',[45 67 48 33].*s+q);
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
      ei_smooth = 1;
    else
      ei_smooth = 3;
    end
    if typecomp == 1 % plot component
         offset = nan_mean(EEG.data(numcompo,:));
         erpimage( EEG.data(numcompo,:)-offset, ones(1,EEG.trials)*10000, EEG.times , ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp');   
    else % plot channel
          if option_computeica  
                  offset = nan_mean(EEG.icaact(numcompo,:));
                  erpimage( EEG.icaact(numcompo,:)-offset, ones(1,EEG.trials)*10000, EEG.times , ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '');   
          else
                
              icaacttmp = (EEG.icaweights(numcompo,:) * EEG.icasphere) ...
                                   * reshape(EEG.data(1:size(EEG.icaweights,1),:,:), EEG.nbchan, EEG.trials*EEG.pnts);
%                     icaacttmp = (EEG.icaweights(numcompo,:) * EEG.icasphere) ...
%                                    * EEG.data(EEG.nbchan, EEG.trials*EEG.pnts);
                  offset = nan_mean(icaacttmp);
                  erpimage( icaacttmp-offset, ones(1,EEG.trials)*10000, EEG.times, ...
                       '', ei_smooth, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '');   
          end;
    end;
    axes(hhh);
    title(sprintf('%s activity \\fontsize{10}(global offset %3.3f)', basename, offset), 'fontsize', 12);
else

    % put title at top of erpimage
    EI_TITLE = 'Continous data';
    axis off
    hh = axes('Units','Normalized', 'Position',[45 62 48 38].*s+q);
    ERPIMAGELINES = 200; % show 200-line erpimage
    while size(EEG.data,2) < ERPIMAGELINES*EEG.srate
       ERPIMAGELINES = round(0.9 * ERPIMAGESLINES);
    end
    if ERPIMAGELINES > 2   % give up if data too small
      if ERPIMAGELINES < 10
         ei_smooth == 1;
      else
        ei_smooth = 3;
      end
      erpimageframes = floor(size(EEG.data,2)/ERPIMAGELINES);
      erpimageframestot = erpimageframes*ERPIMAGELINES;
      eegtimes = linspace(0, erpimageframes-1, EEG.srate/1000);
      if typecomp == 1 % plot component
           offset = nan_mean(EEG.data(numcompo,:));
           erpimage( reshape(EEG.data(numcompo,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar');   
      else % plot channel
            if option_computeica  
                    offset = nan_mean(EEG.icaact(numcompo,:));
                    erpimage( ...
              reshape(EEG.icaact(numcompo,1:erpimageframestot),erpimageframes,ERPIMAGELINES)-offset, ...
                   ones(1,ERPIMAGELINES)*10000, eegtimes , ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar','yerplabel', '');   
            else
%                     icaacttmp = reshape(EEG.icaweights(numcompo,:) * EEG.icasphere) ...
%                                      * reshape(EEG.data, erpimageframes, ERPIMAGELINES);
                        
                        icaacttmp = EEG.icaweights(numcompo,:) * EEG.icasphere ...
                                     *EEG.data(:,1:erpimageframes*ERPIMAGELINES);

                    offset = nan_mean(icaacttmp);
                    erpimage( icaacttmp-offset, ones(1,ERPIMAGELINES)*10000, eegtimes, ...
                         EI_TITLE, ei_smooth, 1, 'caxis', 2/3, 'cbar', 'yerplabel', '');   
            end;
      end
    else
            axis off;
            text(0.1, 0.3, [ 'No erpimage plotted' 10 'for small continuous data']);
    end;
    axes(hhh);
end;	

% plotting spectrum
% -----------------
if ~exist('winhandle')
    winhandle = NaN;
end;
if ~isnan(winhandle)
	h = axes('units','normalized', 'position',[10 25 85 25].*s+q);
    %h = axes('units','normalized', 'position',[5 10 95 35].*s+q); %%%
    %CHANGE!
else
	h = axes('units','normalized', 'position',[10 15 85 30].*s+q);
    %h = axes('units','normalized', 'position',[5 0 95 40].*s+q); %%%
    %CHANGE!
end;
%h = axes('units','normalized', 'position',[45 5 60 40].*s+q);
try
	eeglab_options;
    %next instr added for correct function! Andrea
    option_computeica=1;
	if typecomp == 1
		%[spectra freqs] = spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate, spec_opt{:},'freqrange', [0 45] );
        [spectra freqs] = spectopo( EEG.data(numcompo,:), EEG.pnts, EEG.srate, 'freqrange', [0 45] );
	else 
		if option_computeica  
            
			%[spectra freqs] = spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), spec_opt{:}, 'freqrange', [0 45]);
            % CONTROL ADDED FOR CONTINUOUS DATA
            if size(EEG.data,3)==1 
                EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
            end
            [spectra freqs] = spectopo( EEG.icaact(numcompo,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo),  'freqrange', [0 45]);
		else
			if exist('icaacttmp')~=1, 
                
				icaacttmp = (EEG.icaweights(numcompo,:)*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts); 
			end;
            
			%[spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), spec_opt{:} ,'freqrange', [0 45]);
            [spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,numcompo), 'freqrange', [0 45]);
		end;
	end;
    set(gca,'fontsize',8);
    % set up new limits
    % -----------------
    %freqslim = 50;
	%set(gca, 'xlim', [0 min(freqslim, EEG.srate/2)]);
    %spectra = spectra(find(freqs <= freqslim));
	%set(gca, 'ylim', [min(spectra) max(spectra)]);
    
	%tmpy = get(gca, 'ylim');
    %set(gca, 'ylim', [max(tmpy(1),-1) tmpy(2)]);
	set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 8); 
	set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 8); 
	title('Activity power spectrum', 'fontsize', 12); 
catch
	axis off;
	text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 ' make sure you have the ' 10 'signal processing toolbox']);
end;	

% ----------------------------------------------------------------
% plotting IC properties
% -----------------
if ~exist('winhandle')
    winhandle = NaN;
end;
if ~isnan(winhandle)
	h = axes('units','normalized', 'position',[3 2 95 10].*s+q);
else
	h = axes('units','normalized', 'position',[3 0 95 10].*s+q);
end;

axis off
str='FASTER - Detected as ';
FASTER_reasons = {'High freq ' 'flat spectrum ' 'spatial Kurtosis ' 'Hurst exponent ' 'EOG correl '};
%                     1 Median gradient value, for high frequency stuff
%                     2 Mean slope around the LPF band (spectral)
%                     3 Kurtosis of spatial map
%                     4 Hurst exponent
%                     5 Eyeblink correlations
zlist = zscore(listprops);
for i = 1:size(listprops,2)
    fst(:,i) = min_z(listprops(:,i));
end
reasons = FASTER_reasons(fst(numcompo,:));
if isempty(reasons)
    str = 'FASTER says ok.';
else
    str = [str reasons{:}];
end
% set bar colors
C={[1 0 0],[.6 0 .2],[1 1 0],[0 1 0], [0 1 1]};
hold on
for i = 1:size(listprops,2)
    bar(i,zlist(numcompo,i),'facecolor',C{i});
end
xlim([0 size(listprops,2)+1])
hold on
hline([-3 3],'k');
set(gca, 'ytick',[-3 3],'yticklabel',[-3 3],'xtick',1:5,'xticklabel',FASTER_reasons,'fontsize',10,...
    'color','none');
h = legend (FASTER_reasons,'location','eastoutside');
set(h,'box','off','fontsize',10)
ylabel('zscore')
title(str)



% -----------------------------------------------------------------


% display buttons
% ---------------

if ~isnan(winhandle)
	COLREJ = '[1 0.6 0.6]';
	COLACC = '[0.75 1 0.75]';
	% CANCEL button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

	% VALUE button
	% -------------
	hval  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

	% REJECT button
	% -------------
	status = EEG.reject.gcompreject(numcompo);
	hr = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
				'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
	command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
				'if get(gcbo, ''userdata''),' ...
				'     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
				'else ' ...
				'     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
				'end;' ];					
	set( hr, 'callback', command); 

	% HELP button
	% -------------
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''pop_prop_FST'');');

	% OK button
	% ---------
 	command = [ 'global EEG;' ...
 				'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
 				'EEG.reject.gcompreject(' num2str(numcompo) ') = tmpstatus;' ]; 
	if winhandle ~= 0
	 	command = [ command ...
	 				sprintf('if tmpstatus set(%3.15f, ''backgroundcolor'', %s); else set(%3.15f, ''backgroundcolor'', %s); end;', ...
					winhandle, COLREJ, winhandle, COLACC)];
	end;				
	command = [ command 'close(gcf); clear tmpstatus' ];
	h  = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

	% draw the figure for statistical values
	% --------------------------------------
	index = num2str( numcompo );
	command = [ ...
		'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
		'' ...
		'pos = get(gcf,''Position'');' ...
		'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
		'pos = get(gca,''position'');' ...
		'q = [pos(1) pos(2) 0 0];' ...
		's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
		'axis off;' ...
		''  ...
		'txt1 = sprintf(''(\n' ...
						'Entropy of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ' AND                 \t\t\t----\n\n' ...
					    'Kurtosis of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ') OR                  \t\t\t----\n\n' ...
					    'Kurtosis distibution \t\t\t%2.2f\n' ...
					    '> Rejection threhold\t\t\t%2.2f\n\n' ...
					    '\n' ...
					    'Current thesholds sujest to %s the component\n\n' ...
					    '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
						'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
						'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
		'' ...				
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
		'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
		'clear txt1 q s h pos;' ];
	set( hval, 'callback', command); 
	if isempty( EEG.stats.compenta )
		set(hval, 'enable', 'off');
	end;
	
    % MODIFICA
	%com = sprintf('pop_prop( %s, %d, %d, 0, %s);', inputname(1), typecomp, numcompo, vararg2str( { spec_opt } ) );
    	com = 'idontcare';

else
	%com = sprintf('pop_prop( %s, %d, %d, NaN, %s);', inputname(1), typecomp, numcompo, vararg2str( { spec_opt } ) );
    	com = 'idontcare';

end;




% trim_and_max() - Computes maximum value from vector 'vettore'
% after removing the top 1% of the values
% (to be outlier resistant)
%
% Usage:
%   >> valore=trim_and_max(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result            
%
%
% Author: Andrea Mognon, Center for Mind/Brain Sciences, University of
% Trento, 2009

% Motivation taken from the following comment to our paper:
% "On page 11 the authors motivate the use of the max5 function when computing
% Maximum Epoch Variance because the simple maximum would be too sensitive
% to spurious outliers. This is a good concern, however the max5 function would
% still be sensitive to spurious outliers for very large data sets. In other words, if
% the data set is large enough, one will be very likely to record more than five
% outliers. The authors should use a trimmed max function that computes the
% simple maximum after the top say .1% of the values have been removed from
% consideration. This rejection criteria scales appropriately with the size of the data
% set."

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


function valore=trim_and_max(vettore)


dim=floor(.01*size(vettore,2)); % = 1% of vector length

tmp=sort(vettore);
valore= tmp(length(vettore)-dim);



% trim_and_mean() - Computes average value from vector 'vettore'
% after removing the top .1% of the values
% (to be outlier resistant)
%
% Usage:
%   >> valore=trim_and_mean(vettore);
%
% Inputs:
%   vettore    - row vector
%
% Outputs:
%   valore     - result            
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


function valore=trim_and_mean(vettore)


dim=floor(.01*size(vettore,2)); % = 1% of vector length

tmp=sort(vettore);
valore= mean (tmp(1:(length(vettore)-dim)));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% END  ADJUST   CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   BELOW IS FASTER CODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   END FASTER CODE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



