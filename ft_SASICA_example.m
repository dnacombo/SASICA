addpath('/home/maximilien.chaumon/owncloud/MATLAB/fieldtrip')
ft_defaults

cfg = [];
cfg.dataset = '/network/lustre/iss01/cenir/analyse/meeg/REFERENCE/laurent_l/150707/damier/damier02_tsss.fif';
cfg.trialdef.eventtype = 'STI101';
cfg.trialdef.eventvalue = 3841;
cfg.trialdef.prestim = .1;
cfg.trialdef.poststim = .4;
cfg = ft_definetrial(cfg);
cfg.channel = {'meg' 'BIO*'};
data_orig = ft_preprocessing(cfg);

%%
cfg = [];
cfg.megscale = 1e8;
cfg.gradscale = 1e-1;
[~, cfg.mychan] = chnb('BIO.*',data_orig.label);
cfg.mychanscale = repmat(2,1,numel(cfg.mychan));
cfg.viewmode = 'vertical';
cfg.preproc.demean = 'yes';
ft_databrowser(cfg,data_orig)
%%
cfg = [];
cfg.resamplefs = 100;
data_rs = ft_resampledata(cfg,data_orig);



cfg = [];
cfg.channel = 'megmag';
datamag = ft_selectdata(cfg,data_rs);
cfg.channel = 'meggrad';
datagrad = ft_selectdata(cfg,data_rs);
cfg.channel = 'BIO*';
data_bio = ft_selectdata(cfg,data_rs);

cfg = [];
cfg.method = 'runica';
cfg.numcomponent = 60;
compmag = ft_componentanalysis(cfg,datamag);
compgrad = ft_componentanalysis(cfg,datagrad);

nudatamag = ft_appenddata([],datamag,data_bio);
nudatagrad = ft_appenddata([],datamag,data_bio);

%%
cfg = [];
cfg.autocorr.enable = 1;

toreject = ft_SASICA_neuromag(cfg,compmag,data_rs);
