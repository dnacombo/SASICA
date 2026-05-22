function [heart_IC, meas, aaa_parameters_find_heart_IC, output_for_zscore_corMatrix_ROC, output_for_user] = CARACAS(cfg, comp)

%% DOCSTRING

%%%%%% AUTHOR %%%%%%
% Pierre Champetier (2024) Contact: pi.champetier@gmail.com


% %%%%%% AIM %%%%%%
% This function takes EEG independent components (IC) as inputs, and finds the cardiac IC. --> To do so:
%
% The script detects cardiac events (PQRST) in all IC.
% It computes several features (std/mean ratio, skewness, kurtosis...) for the distribution of RR intervals, R, Q and S amplitudes.
% It identifies the top 3 IC with the lowest values for each feature (or highest values depending on the feature), and gives a score (1 or 0) for the feature.
% It computes a global score (sum of score/nb of feature).
% After sanity check (physiological beat per minute (bpm) and regularity of IC timecourse), cardiac IC are identified based on the global score: -method 1: cardiac IC if global score >= threshold (ex: 0.5) -method 2: cardiac if global score >= mean(all_global scores) + n*std(all_global_scores)


% %%%%%% INPUTS %%%%%%
% 1) cfg --> all the parameters including:
%   -method_chosen --> 'absolute_threshold' (method 1) or 'mean_std' (method 2)
%   -plot_heart_IC --> 1 or 0 (to plot the IC labelled as cardiac)
%   -path_output --> path where you want to save i) the distribution of the scores (among all IC) and ii) timecourse of the identified cardiac IC (used only if plot_heart_IC == 1)
%   -file_info --> name of your recording (to add it in the name of your plot files) (used only if plot_heart_IC == 1)
%   -nb_IC_wanted --> number of IC selected for each metric (kurtosis, skewness...) [default: 3, to select the top 3 IC for each metric]
%   -bpm_min and bpm_max  --> expected heart beat per min, for sanity check [default: 45 and 90]
%   -threshold_cond_IC_method1 --> minimum proportion of conditions that must be met in order that an IC could be considered as a potential heart IC [default: 0.5, so if method_chosen == 'absolute_threshold', an IC must be in the top 3 for at least 50% of the metrics]
%   -threshold_std_method2 --> if method_chosen == 'mean_std', an IC will be considered as a potential heart IC if its proportion of conditions met (i.e., its score) is above mean(all_score) + threshold_std_method2 * std(all_score) [default: 2.5]
%   -min_recording_duration_sec --> minimum duration (in sec) of the IC timecourse (default: 20]
%   -mini_bouts_duration_for_SignalAmplRange --> for sanity check (avoids false positive): the time course of a potential heart IC must be ~regular. The timecourse will be divided into mini-segments of this duration, and we will check that the amplitude between these mini-bouts is ~similar. [default: 10]
% - threshold_regularity_signal_minmax --> For each mini-bout, the averaged signal amplitude is computed. The IC timecourse will be considered as irregular if: (max(Mean_Amp_minibout) - min(Mean_Amp_minibout)) / min(Mean_Amp_minibout) > threshold_regularity_signal_minmax [default: 1.5]
%
% 2) comp --> your IC in FieldTrip format


% %%%%%% USAGE %%%%%%
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);
%
% Example 1 (all default parameters):
% EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',round(dataRank/5)); %Run your ICA with EEGLAB
% comp = eeglab2fieldtrip(EEG, 'comp'); % Convert into FieldTrip format
% cfg = [];
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

% Example 2 (If you want to use other parameters, specify them in cfg):
% cfg = [];
% cfg.nb_IC_wanted = 5; % The heart IC must be in the top 5 of all IC for skewness, kurtosis... of Rampl, RRintervals...
% cfg.bpm_max = 120; % Max physiological bpm
% [rejected_heart_IC, table_heart_IC, method_reject_cardiac_IC] = A_fct_find_cardiac_IC(cfg, comp);

% Rq: If plot_heart_IC == 1, you must specify cfg.path_output and cfg.file_info (output path and the name of your recording to have it in the file name of your plot)


% %%%%%% TOOLBOX USED TO DETECT CARDIAC EVENTS %%%%%%
% Based on R. Sanghavi, F. Chheda, S. Kanchan and S. Kadge, "Detection Of Atrial Fibrillation in Electrocardiogram Signals using Machine Learning," 2021
% 2nd Global Conference for Advancement in Technology (GCAT), 2021, pp. 1-6, doi: 10.1109/GCAT52182.2021.9587664.
% Info : https://fr.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox
% Citation pour cette source: Rohan Sanghavi (2024). ECG SIGNAL PQRST PEAK DETECTION TOOLBOX
% (https://www.mathworks.com/matlabcentral/fileexchange/73850-ecg-signal-pqrst-peak-detection-toolbox), MATLAB Central File


% %%%%%% DEPENDENCIES %%%%%%
% A_fct_test_unif_NEW.m (To test the uniform distribution of cardiac events)
% ECG_PQRST_VERSION_3_PC (Toolbox used to detect cardiac events) --> Rq: I did small modifications of 'compute_fudicial_peaks_live17_c' and 'preprocess_window_ecg' functions to avoid errors with some IC (all modifications are preceeded by a comment "PIERRE").

% add ECG toolbox to the path
% addpath(fullfile(fileparts(which(mfilename)),'ECG_PQRST_VERSION_3_PC'))
addpath(fullfile(fileparts(which(mfilename)),'heart_functions'))
%% CONSTANT PARAMETERS

fs = comp.fsample;

if numel(comp.trial) > 1
    error('data should be continuous');
end

% Window for ERP
window = 0.2*fs; % 0.2s before and 0.2s after the R peak

%% Extract parameters

% method_chosen
if isfield(cfg, 'method_chosen') == 1
    method_chosen = cfg.method_chosen;
else
    method_chosen = 'absolute_threshold'; % absolute_threshold (method 1) or mean_std (method 2)
end

% plot_heart_IC
if isfield(cfg, 'plot_heart_IC')
    plot_heart_IC = cfg.plot_heart_IC;
else
    plot_heart_IC = 0;
end

% path_output and file_info
if plot_heart_IC == 1
    try
        path_output = cfg.path_output;
        file_info = cfg.file_info;
    catch
        error('When plot_heart_IC == 1, you must specified cfg.path_output and cfg.file_info')
    end
end

% nb_IC_wanted
if isfield(cfg, 'nb_IC_wanted') == 1
    nb_IC_wanted = cfg.nb_IC_wanted;
else
    nb_IC_wanted = 3; % default
end

% bpm_min
if isfield(cfg, 'bpm_min') == 1
    bpm_min = cfg.bpm_min;
else
    bpm_min = 35;
end

% bpm_max
if isfield(cfg, 'bpm_max') == 1
    bpm_max = cfg.bpm_max;
else
    bpm_max = 90;
end

% threshold_cond_IC_method1
if isfield(cfg, 'threshold_cond_IC_method1') == 1
    threshold_cond_IC_method1 = cfg.threshold_cond_IC_method1;
else
    threshold_cond_IC_method1 = 0.6;
end

% threshold_std_method2
if isfield(cfg, 'threshold_std_method2') == 1
    threshold_std_method2 = cfg.threshold_std_method2;
else
    threshold_std_method2 = 2.5;
end

% threshold_regularity_signal_minmax
if isfield(cfg, 'threshold_regularity_signal_minmax') == 1
    threshold_regularity_signal_minmax = cfg.threshold_regularity_signal_minmax;
else
    threshold_regularity_signal_minmax = 1.5;
end

% min_recording_duration_sec
if isfield(cfg, 'min_recording_duration_sec') == 1
    min_recording_duration_sec = cfg.min_recording_duration_sec;
else
    min_recording_duration_sec = 20;
end

% mini_bouts_duration_for_SignalAmplRange
if isfield(cfg, 'mini_bouts_duration_for_SignalAmplRange') == 1
    mini_bouts_duration_for_SignalAmplRange = cfg.mini_bouts_duration_for_SignalAmplRange;
else
    mini_bouts_duration_for_SignalAmplRange = 10;
end

% IC to not analyze (because user has already labelled it for instance)
if isfield(cfg, 'IC_to_not_analyze') == 1
    IC_to_not_analyze = cfg.IC_to_not_analyze;
    % Check user has given number(s) for IC to not analyze
    if isnumeric(IC_to_not_analyze) == 0 && ~isequal(IC_to_not_analyze, [])
        error('Error in cfg.IC_to_not_analyze: please enter a number (or use [] /do not use cfg.IC_to_not_analyze).')
    end
    % Check user has given an IC that exists
    if sum(ismember(IC_to_not_analyze, [1:length(comp.label)])) ~= length(IC_to_not_analyze)
        error('Error in cfg.IC_to_not_analyze: you want to remove IC that does not exist.')
    end
else
    IC_to_not_analyze = [];
end

%% Pre-defined col_names (as it might be needed in "Extract parameters section if recording_duration_sec < min_recording_duration_sec)

col_names = {'IC_bug', 'IC_bpm_all', 'IC_bpm_ok', 'SignalAmpl_range_all', 'SignalAmpl_range_ok',...
    'IC_ERPampl_median',  'IC_ERPampl_std', 'IC_ERPampl_skew', 'IC_ERPampl_kurt', 'IC_ERPampl_chi2stat', 'IC_ERPampl_std_median',...
    'IC_RR_std_mean', 'IC_RR_skew', 'IC_RR_kurt', 'IC_RR_chi2stat', ...
    'IC_Rampl_std_mean', 'IC_Rampl_skew', 'IC_Rampl_kurt', 'IC_Rampl_chi2stat', ...
    'IC_heart_method2', 'IC_heart_method1', 'Prop_cond_ok_heart_IC', 'Prop_cond_ok_all_IC', 'Prop_cond_ok_other_IC', 'Max_prop_cond_ok_other_IC', 'bpm_final_IC', 'Nbr_sec_ICA'};


%% Extract parameters

aaa_parameters_find_heart_IC = [];
aaa_parameters_find_heart_IC.nb_IC_wanted = nb_IC_wanted;
aaa_parameters_find_heart_IC.bpm_min = bpm_min;
aaa_parameters_find_heart_IC.bpm_max = bpm_max;
aaa_parameters_find_heart_IC.threshold_regularity_signal_minmax = threshold_regularity_signal_minmax;
aaa_parameters_find_heart_IC.threshold_cond_IC_method1 = threshold_cond_IC_method1;
aaa_parameters_find_heart_IC.threshold_std_method2 = threshold_std_method2;
aaa_parameters_find_heart_IC.min_recording_duration_sec = min_recording_duration_sec;
aaa_parameters_find_heart_IC.mini_bouts_duration_for_SignalAmplRange = mini_bouts_duration_for_SignalAmplRange;
aaa_parameters_find_heart_IC.IC_to_not_analyze = IC_to_not_analyze;

% Checking that recording_duration_sec > mini_bouts_duration_for_SignalAmplRange
sample_tot = 0;
for i = 1:length(comp.trial)
    blocs = isnan(comp.trial{i}(1,:));
    d = diff(blocs);
    s = find(d == 1);
    sample_tot = sample_tot + size(comp.trial{i},2) - round(0.5 * fs) * numel(s) - sum(blocs);
    % assuming HB detection isn't working properly on edges so we discard
    % 0.5s for each trial in this sample_tot.
    % In case the trial is continuous, but has nans in it (if trial data
    % has been turned to continuous prior to entering CARACAS), we discard
    % this 0.5s for each block.
    % This is used only for bpm computation at the end.
end
recording_duration_sec = sample_tot / fs;

if recording_duration_sec < mini_bouts_duration_for_SignalAmplRange
    error('Error: attempt to divide the recording into bouts with longer than the original recording (for the Signal Ampl Range checking). Change the value of mini_bouts_duration_for_SignalAmplRange variable in the function.')
end

% Checking that recording_duration_sec > min_recording_duration_sec
if recording_duration_sec < min_recording_duration_sec
    warning('recording_duration_sec < min_recording_duration_sec')
    table_cardiac_IC = table({[1:length(comp.label)]}, {[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]},{[]});
    table_cardiac_IC.Properties.VariableNames = col_names;
    heart_IC = [];
else



    %% OLD: Extract timecourse of each IC for all 1s-mini trials

    IC_timecourse = [];
    for i = 1:numel(comp.trial)
        IC_timecourse = [IC_timecourse, comp.trial{1,i}];
    end

    %% NEW: Recreate continuous comp data
    % cfgtmp = [];
    % if isfield(comp, 'sampleinfo')
    %     cfgtmp.trl = [1, comp.sampleinfo(end,2), 0];
    % else
    %     cfgtmp.trl = [1, sum(cellfun(@numel,comp.time)), 0]; % Start, End, Offset
    % end
    % comp_continu = ft_redefinetrial(cfgtmp, comp);
    %
    % Extract timecourse
    timecourse_all_IC = comp.trial{1};

    %% DETECT CARDIAC EVENTS

    plot_data = 0; % (see FUNCTION DOCUMENTATION)

    IC_bug = [];
    IC_with_less_than_10_R_peak = [];

    bpm_all_IC = [];

    % RR features
    RR_mean_all = [];
    RR_std_all = [];
    RR_skew_all = [];
    RR_kurt_all = [];
    RR_ratio_std_mean_all_IC = [];
    RR_chi2stat_all_IC = [];

    % Rampl features
    Rampl_kurt_all_IC = [];
    Rampl_skew_all_IC = [];
    Rampl_ratio_std_mean_all_IC = [];
    Rampl_chi2stat_all_IC = [];

    SignalAmpl_range_all_IC = [];

    % Qampl features
    Qampl_kurt_all_IC = [];
    Qampl_skew_all_IC = [];
    Qampl_ratio_std_mean_all_IC = [];
    Qampl_chi2stat_all_IC = [];

    Sampl_kurt_all_IC = [];
    Sampl_skew_all_IC = [];
    Sampl_ratio_std_mean_all_IC = [];
    Sampl_chi2stat_all_IC = [];

    ERP_average_for_each_IC = [];
    ERP_details_for_each_IC = [];
    ERPampl_median = [];
    ERPampl_std = [];
    ERPampl_std_median = [];
    ERPampl_skew = [];
    ERPampl_kurt = [];
    ERPampl_chi2stat = [];



    nbr_cardiac_event = [];
    IC_not_cardiac_bc_PQstd = [];
    IC_not_cardiac_bc_RRstd = [];
    IC_not_cardiac_bc_Ramplstd = [];
    IC_not_cardiac_bc_skewcorr = NaN([1,numel(comp.label)]);

    meas = [];
    for comp_iter = 1:length(comp.label)
        % Skip if user don't want to analyze this IC
        if ismember(comp_iter, IC_to_not_analyze)
            ERP_average_for_each_IC = [ERP_average_for_each_IC; repmat(NaN,1,window*2+1);];
            ERP_details_for_each_IC{1,comp_iter} = {};
            ERPampl_median = [ERPampl_median; comp_iter, NaN];

            bpm_all_IC = [bpm_all_IC; comp_iter, NaN];
            SignalAmpl_range_all_IC = [SignalAmpl_range_all_IC; comp_iter, NaN];

            RR_mean_all = [RR_mean_all; comp_iter, NaN];
            RR_std_all = [RR_std_all; comp_iter, NaN];

            RR_ratio_std_mean_all_IC = [RR_ratio_std_mean_all_IC; comp_iter, NaN];
            RR_skew_all = [RR_skew_all; comp_iter, NaN];
            RR_kurt_all = [RR_kurt_all; comp_iter, NaN];
            RR_chi2stat_all_IC = [RR_chi2stat_all_IC; comp_iter, NaN];
            Rampl_ratio_std_mean_all_IC = [Rampl_ratio_std_mean_all_IC; comp_iter, NaN];
            Rampl_skew_all_IC = [Rampl_skew_all_IC; comp_iter, NaN];
            Rampl_kurt_all_IC = [Rampl_kurt_all_IC; comp_iter, NaN];
            Rampl_chi2stat_all_IC = [Rampl_chi2stat_all_IC; comp_iter, NaN];
            Qampl_ratio_std_mean_all_IC = [Qampl_ratio_std_mean_all_IC; comp_iter, NaN];
            Qampl_skew_all_IC = [Qampl_skew_all_IC; comp_iter, NaN];
            Qampl_kurt_all_IC = [Qampl_kurt_all_IC; comp_iter, NaN];
            Qampl_chi2stat_all_IC = [Qampl_chi2stat_all_IC; comp_iter, NaN];
            Sampl_ratio_std_mean_all_IC = [Sampl_ratio_std_mean_all_IC; comp_iter, NaN];
            Sampl_skew_all_IC = [Sampl_skew_all_IC; comp_iter, NaN];
            Sampl_kurt_all_IC = [Sampl_kurt_all_IC; comp_iter, NaN];
            Sampl_chi2stat_all_IC = [Sampl_chi2stat_all_IC; comp_iter, NaN];
            continue
        end

        % Select IC iter
        % ecg_1 = IC_timecourse(comp_iter,:);

        % Detect cardiac events
        % try
        % ecg_f = preprocess_window_ecg(ecg_1, fs);
        % ecg_f = ecg_1;
        % [locs_P,locs_Q,locs_R,locs_S,locs_T] = compute_fudicial_peaks_live17_c(ecg_f, fs, plot_data); % Why the nb of samples in ecg_f is lower than in ecg_1?? For some IC it's super low....

        cfg_peak = [];
        % cfg_peak.plotall         = 1;
        % cfg_peak.plotthresh      = 0;
        % cfg_peak.plotbeat        = 1;
        % cfg_peak.plotcorr        = 0;
        % cfg_peak.plotfinal       = 1;
        cfg_peak.channel = comp.label{comp_iter};
        cfg_peak.corthresh = 0.2;
        cfg_peak.absPT = 0;

        % try
        % cfg = [];
        % cfg.viewmode  = 'component';
        % cfg.layout    = evalin('base','layout');
        % ft_databrowser(cfg, comp);

        [HeartBeats] = heart_peak_detect(cfg_peak,comp);


        locs_P = [HeartBeats.P_sample];
        locs_Q = [HeartBeats.Q_sample];
        locs_R = [HeartBeats.R_sample];
        locs_S = [HeartBeats.S_sample];
        locs_T = [HeartBeats.T_sample];





        meas(comp_iter).PQ = gimme_interval('PQ', HeartBeats);
        meas(comp_iter).QS = gimme_interval('QS', HeartBeats);
        meas(comp_iter).ST = gimme_interval('ST', HeartBeats);
        meas(comp_iter).PR = gimme_interval('PR', HeartBeats);
        meas(comp_iter).RT = gimme_interval('RT', HeartBeats);
        meas(comp_iter).PT = gimme_interval('PT', HeartBeats);

        meas(comp_iter).sk = HeartBeats.sk;
        if meas(comp_iter).sk < 2
            IC_not_cardiac_bc_skewcorr = [IC_not_cardiac_bc_skewcorr, comp_iter];
        end
        meas(comp_iter).ku = HeartBeats.ku;
        % todo mettre un critère pour ku

        if meas(comp_iter).PQ > 1 / 3
            IC_not_cardiac_bc_PQstd = [IC_not_cardiac_bc_PQstd, comp_iter];
        end




        % RR intervals
        %%%%%%%%%%%%%%
        RR_intervals_sec = diff([HeartBeats.R_time]);
        % figure;
        % hist(RR_intervals_sec);

        % Remove lowest 5% and highest 20% lowest values (that biase the std)
        low_threshold = prctile(RR_intervals_sec, 0);
        high_threshold = prctile(RR_intervals_sec, 70);

        filtered_RR_intervals = RR_intervals_sec(RR_intervals_sec >= low_threshold & RR_intervals_sec <= high_threshold);

        % figure;
        % hist(filtered_RR_intervals);
        meas(comp_iter).RR = std(filtered_RR_intervals) / mean(filtered_RR_intervals);
        if meas(comp_iter).RR > 1 / 3
            IC_not_cardiac_bc_RRstd = [IC_not_cardiac_bc_RRstd, comp_iter];
        end


        % Rampl
        %%%%%%%
        time_course_IC_iter = timecourse_all_IC(comp_iter,:);
        Rampl = time_course_IC_iter(locs_R);
        Rampl = abs(Rampl);

        % figure;
        % hist(Rampl);

        % Remove lowest 15% and highest 15% lowest values (that biase the std)
        low_threshold = prctile(Rampl, 15);  % 5th percentile
        high_threshold = prctile(Rampl, 85); % 95th percentile

        filtered_Rampl = Rampl(Rampl >= low_threshold & Rampl <= high_threshold);

        % figure;
        % hist(filtered_Rampl);
        meas(comp_iter).Rampl = std(filtered_Rampl) / abs(mean(filtered_Rampl));
        if meas(comp_iter).Rampl > 1 /3
            IC_not_cardiac_bc_Ramplstd = [IC_not_cardiac_bc_Ramplstd, comp_iter];
        end



        % catch ME
        %     warning('error in heart_peak_detect')
        %     % nbr_cardiac_event = [nbr_cardiac_event; comp_iter 0];
        %     locs_P = [];
        % end

        nbr_cardiac_event = [nbr_cardiac_event; comp_iter length(locs_P)];




        %%  Metric to check homogeneous SignalAmpl across the recording
        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

        % Remove NaN
        time_course_comp_iter = timecourse_all_IC(comp_iter,:);
        time_course_comp_iter(find(isnan(time_course_comp_iter))) = [];


        % Divide the recordings into mini segments and check that the abs(mean) is similar for each segment
        nbr_mini_segments = round((length(time_course_comp_iter) / fs) / mini_bouts_duration_for_SignalAmplRange);

        % Extract the absolute value of the signal
        array = abs(time_course_comp_iter);

        % Determine the size of each small array
        sizeOfSmallArray = floor(length(array)/ nbr_mini_segments);

        % Preallocate cell array for the small arrays
        smallArrays = cell(1, nbr_mini_segments);

        % Fill the small arrays
        for i = 1:nbr_mini_segments
            startIdx = (i-1) * sizeOfSmallArray + 1;
            if i < nbr_mini_segments
                endIdx = i * sizeOfSmallArray;
            else
                % Last array takes any remaining elements
                endIdx = length(array);
            end
            smallArrays{i} = array(startIdx:endIdx);
        end

        % Compute mean for each segment
        SignalAmpl_mean_segments = [];
        for i = 1:nbr_mini_segments
            SignalAmpl_mean_segments = [SignalAmpl_mean_segments, range(smallArrays{i})];
        end

        % Compute a metric to evaluate if the signal is regular over the recording (max-min)/mean
        SignalAmpl_range = (max(SignalAmpl_mean_segments) - min(SignalAmpl_mean_segments)) / min(SignalAmpl_mean_segments);
        SignalAmpl_range_all_IC = [SignalAmpl_range_all_IC; comp_iter, SignalAmpl_range];

        meas(comp_iter).Ampl_var = SignalAmpl_range;
    end



    %% Find heart IC


    % 1) bpm
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % Option 1: with data-driven threshold
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    threshold_IC_4std = median(nbr_cardiac_event(:,2)) +  4*std(nbr_cardiac_event(:,2));
    threshold_IC_3std = median(nbr_cardiac_event(:,2)) +  3*std(nbr_cardiac_event(:,2));


    threshold_IC = threshold_IC_3std;
    % heart_IC = nbr_cardiac_event(find(nbr_cardiac_event(:,2) > threshold_IC),1);


    % Option 2: with physiological bpm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nbr_cardiac_event(:,3) = nbr_cardiac_event(:,2)/(recording_duration_sec /60);
    heart_IC = nbr_cardiac_event(nbr_cardiac_event(:,3) >= bpm_min & nbr_cardiac_event(:,3) <= bpm_max, :);
    heart_IC = heart_IC(:,1);

    [meas.bpm] = rep2struct(nbr_cardiac_event(:,3));

    % 2) Remove potential cardiac if too high std for RR interval or Rampl
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

    % For QS interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_PQstd));
    % heart_IC(idx_IC_with_too_high_std) = [];

    % % For QR interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_QRstd));
    % heart_IC(idx_IC_with_too_high_std) = [];
    %
    % % For RS interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_RSstd));
    % heart_IC(idx_IC_with_too_high_std) = [];
    %
    % % For ST interval
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_STstd));
    % heart_IC(idx_IC_with_too_high_std) = [];

    % For RR interval
    idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_RRstd));
    heart_IC(idx_IC_with_too_high_std) = [];

    % For Rampl
    % idx_IC_with_too_high_std = find(ismember(heart_IC, IC_not_cardiac_bc_Ramplstd));
    % heart_IC(idx_IC_with_too_high_std) = [];

    % skewcorr
    idx_IC_with_too_low_skewcorr = find(ismember(heart_IC, IC_not_cardiac_bc_skewcorr));
    heart_IC(idx_IC_with_too_low_skewcorr) = [];


    % 3) Remove potential cardiac if too high SignalAmpl range
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % idx_SignalAmpl_range = find(SignalAmpl_range_all_IC(:,2) > threshold_regularity_signal_minmax);
    % IC_not_cardiac_bc_Ramplstd = SignalAmpl_range_all_IC(idx_SignalAmpl_range,1);
    % 
    % idx_IC_with_too_high_SignalAmpl_range = find(ismember(heart_IC, IC_not_cardiac_bc_Ramplstd));
    % heart_IC(idx_IC_with_too_high_SignalAmpl_range) = [];


    %% Save output

    output_for_zscore_corMatrix_ROC =[];

    output_for_user = [];
    output_for_user.threshold_IC_3std = threshold_IC_3std;
    output_for_user.threshold_IC_4std = threshold_IC_4std;

    output_for_user.threshold_IC = threshold_IC;
    output_for_user.IC_nbr_cardiac_event = nbr_cardiac_event(:,2);

    if isempty(heart_IC) == 0
        output_for_user.heart_IC = heart_IC;
        output_for_user.heart_IC_nbr_cardiac_event = nbr_cardiac_event(find(nbr_cardiac_event(:,2) > threshold_IC),2);
    else
        output_for_user.heart_IC = NaN;
        output_for_user.heart_IC_nbr_cardiac_event = NaN;
    end
    if plot_heart_IC

        figure(94480); clf;
        set(gcf,'WindowState', 'maximized'); % Open a maximized figure window

        %         %%%%%%%%%%%%%%% Several plots per window %%%%%%%%%%%%%
        nbr_plot_per_window = 9; % If it is changed, change also subplot parameters
        nbr_IC = length(comp.label);

        for i = 1 : nbr_plot_per_window : nbr_IC


            IC_start = i;
            IC_stop =  min(i + nbr_plot_per_window - 1, nbr_IC); % Make sure we don't go beyond n

            figure(94480); clf;
            for IC_to_plot = IC_start:IC_stop
                position_idx = mod(IC_to_plot, nbr_plot_per_window); % 10 --> 1, 11 --> 2...
                if position_idx == 0
                    position_idx = 9;
                end
                subplot(3, 3, position_idx);

                % Generate a plot
                timetoplot = [0 20]; % <<- parametrize?
                t = timepts(timetoplot, comp.time{1});
                plot(comp.time{1}(t), comp.trial{1}(IC_to_plot,t), ifelse(ismember(IC_to_plot, heart_IC),'r',''));

                % Add a title and grid
                title(['IC ' num2str(IC_to_plot)])
                grid on;

                % Optional: Label axes (only if necessary)
                xlabel('Time (s)');
                ylabel('Voltage');
            end
            % %%%%% SAVE fig as .png %%%%%%%%
            filename = strcat(cfg.path_output, '/IC_timecourse/', file_info, '_IC ', num2str(IC_start), '_to_IC_', num2str(IC_stop), '.png');
            mymkdir(fileparts(filename))
            saveas(gcf, filename);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
    end

end


end


function interval = gimme_interval(l, HeartBeats)

intervals_sec = [HeartBeats.([l(1) '_time'])] - [HeartBeats.([l(2) '_time'])];

% figure;
% hist(intervals_sec)

low_threshold = prctile(intervals_sec, 15);
high_threshold = prctile(intervals_sec, 85);

filtered_intervals = intervals_sec(intervals_sec >= low_threshold & intervals_sec <= high_threshold);

% figure;
% hist(filtered_intervals);
interval = std(filtered_intervals)/(abs(mean(filtered_intervals)));
end
