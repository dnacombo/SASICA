function features = extract_time_freq_features(EEG)
%                             - 1st row: Average Local Skewness
%                             - 2nd row: lambda
%                             - 3rd row: Band Power 8 - 13 Hz
%                             - 4rd row: Fit Error
%
data = EEG.data(EEG.icachansind,:,:);
fs = EEG.srate; %sampling frequency

% transform epoched data into continous data
data = data(:,:);

%downsample (to 100-200Hz)
factor = max(floor(EEG.srate/100),1);
data = data(:, 1:factor:end);
fs = round(fs/factor);

%compute icaactivation and standardise variance to 1
icacomps = (EEG.icaweights * EEG.icasphere * data)';
icacomps = icacomps./repmat(std(icacomps,0,1),length(icacomps(:,1)),1);
icacomps = icacomps';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate featues
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ic=1:size(icacomps,1)  %for each component
    fprintf('.');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Proc Spectrum for Channel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [pxx, freq] = pwelch(icacomps(ic,:), ones(1, fs), [], fs, fs);
    pxx = 10*log10(pxx * fs/2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The average log band power between 8 and 13 Hz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = 0;
    for i = 8:13
        p = p + pxx(find(freq == i,1));
    end
    Hz8_13 = p / (13-8+1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % lambda and FitError: deviation of a component's spectrum from
    % a protoptypical 1/frequency curve
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1.x = 2; %first point: value at 2 Hz
    p1.y = pxx(find(freq == p1.x,1));

    p2.x = 3; %second point: value at 3 Hz
    p2.y = pxx(find(freq == p2.x,1));

    %third point: local minimum in the band 5-13 Hz
    p3.y = min(pxx(find(freq == 5,1):find(freq == 13,1)));
    p3.x = freq(find(pxx == p3.y,1));

    %fourth point: min - 1 in band 5-13 Hz
    p4.x = p3.x - 1;
    p4.y = pxx(find(freq == p4.x,1));

    %fifth point: local minimum in the band 33-39 Hz
    p5.y = min(pxx(find(freq == 33,1):find(freq == 39,1)));
    p5.x = freq(find(pxx == p5.y,1));

    %sixth point: min + 1 in band 33-39 Hz
    p6.x = p5.x + 1;
    p6.y = pxx(find(freq == p6.x,1));

    pX = [p1.x; p2.x; p3.x; p4.x; p5.x; p6.x];
    pY = [p1.y; p2.y; p3.y; p4.y; p5.y; p6.y];

    myfun = @(x,xdata)(exp(x(1))./ xdata.^exp(x(2))) - x(3);
    xstart = [4, -2, 54];
    fittedmodel = lsqcurvefit(myfun,xstart,double(pX),double(pY), [], [], optimset('Display', 'off'));

    %FitError: mean squared error of the fit to the real spectrum in the band 2-40 Hz.
    ts_8to15 = freq(find(freq == 8) : find(freq == 15));
    fs_8to15 = pxx(find(freq == 8) : find(freq == 15));
    fiterror = log(norm(myfun(fittedmodel, ts_8to15)-fs_8to15)^2);

    %lambda: parameter of the fit
    lambda = fittedmodel(2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Averaged local skewness 15s
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    interval = 15;
    abs_local_scewness = [];
    for i=1:interval:length(icacomps(ic,:))/fs-interval
        abs_local_scewness = [abs_local_scewness, abs(skewness(icacomps(ic, i * fs:(i+interval) * fs)))];
    end

    if isempty(abs_local_scewness)
        error('MARA needs at least 15ms long ICs to compute its features.')
    else
        mean_abs_local_scewness_15 = log(mean(abs_local_scewness));
    end;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Append Features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    features(:,ic)= [mean_abs_local_scewness_15, lambda, Hz8_13, fiterror];
end
disp('.');


