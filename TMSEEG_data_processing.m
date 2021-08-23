
%% Import data 

%Load NeurOne EEG data
session_folder=; %insert path to session's NeuroOne folder
session_number=; %insert path to session's measurement 
EEG_data = module_read_neurone(session_folder, session_number); %(Toolbox: NeurOne Tools for Matlab) 

EEGchannels =  setxor(fieldnames(EEG_data.signal),{'APBr';'FDIr'}); %define EMG channels
EMGchannels =  setxor(fieldnames(EEG_data.signal),EEGchannels); %define EEG channels

%Create variables for EEG and EMG whole session time course
for ch = 1:length(EEGchannels)
    data_eeg(ch,:) = eval(['EEG_data.signal.' EEGchannels{ch} '.data']);
end
for ch = 1:length(EMGchannels)
    data_emg(ch,:) = eval(['EEG_data.signal.' EMGchannels{ch} '.data']);
end

markers={'A', 'B'}; %TMS pulse markers
idx_test_markers=[];
for mm = 1:length(markers)
    idx_test_markers = [idx_test_markers; find(ismember(EEG_data.markers.port,markers{mm}))];
end
%extract TMS markers' time poits
idx_test_markers  = sort(idx_test_markers);
sample_markers = EEG_data.markers.index(idx_test_markers);
fsample = EEG_data.properties.samplingRate;

%Create variables for EEG and EMG epoched time course (trials: -1000ms to 1500ms around TMS) 
epochs_eeg = zeros(size(data_eeg,1),...
    length(sample_markers(1)-round(1*fsample):sample_markers(1)+round(1.5*fsample)),length(sample_markers));
epochs_emg = zeros(size(data_emg,1),...
    length(sample_markers(1)-round(1*fsample):sample_markers(1)+round(1.5*fsample)),length(sample_markers));
for i = 1:length(sample_markers)
    epochs_eeg(:,:,i) = data_eeg(:,sample_markers(i)-round(1*fsample):sample_markers(i)+round(1.5*fsample));
    epochs_emg(:,:,i) = data_emg(:,sample_markers(i)-round(1*fsample):sample_markers(i)+round(1.5*fsample));   
end

%time couse axis (x-axis)
time_ax = linspace(-1,1.5,size(epochs_eeg,2));

%exclude TMS artifact [-5ms to 20ms]
epochs_eeg_nan = epochs_eeg;
epochs_eeg_nan(:,find(time_ax>=-0.005 & time_ax <=0.020),:) = NaN;

%Change markers: 0=SHAM, 1=REAL TMS
for i=1:length(EEG_data.markers.port)
    if EEG_data.markers.port{i}=='A'; EEG_data.markers.port{i}='1';
    elseif EEG_data.markers.port{i}=='B'; EEG_data.markers.port{i}='0';
    end
end

%Create FieldTrip structure for EEG
eeg_ft = [];
eeg_ft.time = mat2cell(repmat(time_ax,size(epochs_eeg,3),1),ones(size(epochs_eeg_nan,3),1))';
eeg_ft.trial = squeeze(mat2cell(epochs_eeg_nan,size(epochs_eeg_nan,1),size(epochs_eeg_nan,2),ones(size(epochs_eeg_nan,3),1)))';
eeg_ft.fsample = fsample;
eeg_ft.label = EEGchannels;
eeg_ft.sampleinfo(:,1) = sample_markers-repmat(round(0.5*fsample),length(sample_markers),1);
eeg_ft.sampleinfo(:,2) = sample_markers+repmat(round(1*fsample),length(sample_markers),1);
eeg_ft.trialinfo = cellfun(@(x) str2double(x),EEG_data.markers.port(idx_test_markers));
%Create FieldTrip structure for EMG
emg_ft=eeg_ft;
emg_ft.trial=squeeze(mat2cell(epochs_emg,size(epochs_emg,1),size(epochs_emg,2),ones(size(epochs_emg,3),1)))';
emg_ft.label={'APBr';'FDIr'};
emg_ft.time = mat2cell(repmat(time_ax,size(epochs_emg,3),1),ones(size(epochs_emg,3),1))';

%Interpolate excluded TMS in EEG
cfg = [];
cfg.method      ='pchip';
cfg.prewindow   = 1; 
cfg.postwindow  = 1; 
eeg_data_interpolated = ft_interpolatenan(cfg, eeg_ft);

%Downsample to 1000Hz
cfg = [];
cfg.resamplefs = 1000; 
eeg_data_downsampled = ft_resampledata(cfg,eeg_data_interpolated);
emg_data_downsampled = ft_resampledata(cfg,emg_ft);

%Baseline correction EEG
cfg=[];
cfg.demean     = 'yes'; 
cfg.baselinewindow = [-0.1 -0.05 ] ;
eeg_data_processed = ft_preprocessing(cfg, eeg_data_downsampled);
eeg_data_processed.sampleinfo=eeg_ft.sampleinfo;

%Baseline correction EMG
cfg.detrend    = 'yes'; 
emg_data_processed = ft_preprocessing(cfg, emg_data_downsampled);
cfg=[];
cfg.latency = [0.020 0.060];
emg_data_processed = ft_selectdata(cfg, emg_data_downsampled);
emg_data_processed.sampleinfo=emg_ft.sampleinfo;


% Prepare neighbours (for plotting and statistics)
EEG_cap_template=; %insert path to EEG cap standard 10-10 layout (FieldTrip folder)
EEGchannels = eeg_data_processed.label;
elec       = ft_convert_units(ft_read_sens(EEG_cap_template),'mm');
index = find(any(cell2mat(cellfun(@strcmpi,EEGchannels',num2cell(repmat(elec.label,1,numel(EEGchannels)),1),'uniformoutput',false)),2));
elec.chanpos    = elec.chanpos(index,:);
elec.elecpos    = elec.elecpos(index,:);
elec.label      = elec.label(index,:);
neighbour_layout.elec       = elec;
cfg  = [];
neighbour_layout.layout  = ft_prepare_layout(cfg,neighbour_layout);
cfg             = [];
cfg.method      = 'triangulation';
cfg.elec        = elec;
cfg.channels    = 'EEG';
neighbour_layout.neighbor   = ft_prepare_neighbours(cfg);
%% Visual rejection of MEP trials

%create variables forn info on rejected data
reject_info = [];
reject_info.trials = []; 
reject_info.MEP = []; 
reject_info.chans  = {};
reject_info.total_trials  = []; 

% Visual detection of MEPs
cfg                 = [];
cfg.method          = 'trial';
cfg.preproc.demean = 'yes';
[rejectMEP] = ft_rejectvisual(cfg, emg_data_processed);
reject_info.MEP   = find(ismember(emg_data_processed.sampleinfo(:,1), rejectMEP.cfg.artfctdef.trial.artifact(:,1)));

%% Visual rejection of EEG trials

cfg                 = [];
cfg.method          = 'trial';
cfg.trials=find(~ismember([1:length(eeg_data_processed.trial)], reject_info.MEP));
cfg.layout = neighbour_layout;
cfg.preproc.demean = 'yes';
[rejectEEGtrials] = ft_rejectvisual(cfg, eeg_data_processed);

reject_info.trials = find(ismember(eeg_data_processed.sampleinfo(:,1), rejectEEGtrials.cfg.artfctdef.trial.artifact(:,1)));
reject_info.chans = setxor(eeg_data_processed.label,rejectEEGtrials.label);
reject_info.total_trials=unique([reject_info.trials; reject_info.MEP])  ; 

% Remove detected bad trials and channels
cfg=[];
cfg.trials=find(~ismember([1:length(eeg_data_processed.trial)], reject_info.total_trials));
cfg.channel=neighbour_layout.elec.label(find(~ismember(neighbour_layout.elec.label, reject_info.chans)));
data_cleaned=ft_selectdata(cfg,eeg_data_processed);

%% 1st ROUND of ICA
cfg = [];
cfg.channels={'all'};
cfg.method = 'fastica';  
cfg.fastica.approach = 'symm'; 
cfg.fastica.g = 'gauss';
cfg.fastica.interactivePCA = 'off';
ICA_comp1 = ft_componentanalysis(cfg, data_cleaned); 

%Plot and Identify the bad components
close all
nfig = ceil(length(ICA_comp1.label)./8);
topoICA = [];
topoICA.time = 1;
topoICA.label = ICA_comp1.topolabel;
topoICA.fsample = 1000;
topoICA.dimord = 'chan_time';
cfg =[];
cfg.layout = neighbour_layout.layout;
cfg.parameter = 'topography';
cfg.comment = 'no'; 
comp = 0; 
for ff = 1 :nfig
    figure
    pp = 0; 
    while pp< 4 && comp<length(ICA_comp1.label)
        pp = pp+1; 
        comp = comp+1;
        topoICA.topography = ICA_comp1.topo(:, comp);        
        subplot(4, 5, 1+5*(pp-1))  
        ft_topoplotER(cfg,topoICA);
        colormap jet
        title(['ICA comp'  num2str(comp)], 'FontSize', 18)
        timecourse = cell2mat(cellfun(@(x) x(comp,:),ICA_comp1.trial,'UniformOutput',false)');
        subplot(4, 5, 2+5*(pp-1))
        plot(ICA_comp1.time{1},mean(timecourse,1),'LineWidth',2,'Color',[64 64 64]./255)
        xlim([-0.200 0.400])
        ylim([-1000 1000])
        title(['Average timecourse ICA comp'  num2str(comp)])     
        comp = comp+1;
        topoICA.topography = ICA_comp1.topo(:, comp);        
        subplot(4, 5, 4+5*(pp-1))
        ft_topoplotER(cfg,topoICA);
        colormap jet
        title(['ICA comp'  num2str(comp)], 'FontSize', 18)       
        subplot(4, 5, 5+5*(pp-1))
        timecourse = cell2mat(cellfun(@(x) x(comp,:),ICA_comp1.trial,'UniformOutput',false)');
        plot(ICA_comp1.time{1},mean(timecourse,1),'LineWidth',2,'Color',[64 64 64]./255)
        xlim([-0.200 0.400])
        ylim([-1000 1000])
        title(['Average timecourse ICA comp'  num2str(comp)])
         set(gcf,'Position', get(0, 'Screensize'))
    end    
end


  %% Signal bad components from 1st round
%Reject componments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reject_comp1= []; 
%%%%%%%%%%%%%%%%%% 49%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
cfg.component = reject_comp1;  
cfg.updatesens   = 'yes';
cfg.baselinewindow  = [-0.5 -0.05]; 
cfg.demean     = 'yes'; 
[pos_ICA1] = ft_rejectcomponent(cfg,ICA_comp1, data_cleaned);
ICA_comp1.excluded=reject_comp1;

%% 2nd ROUND of ICA

%%% pre-processing data pre ICA 2nd round
cfg=[];
cfg.bpfilter   = 'yes';
cfg.bpfreq     = [0.5 100]; %BAND PASS
cfg.bpfiltord  = 3; % Order of the filter
cfg.bpfilttype = 'but'; % Butterworth
pre_ICA2 = ft_preprocessing(cfg, pos_ICA1);    

%%% --> 2nd ROUND of ICA
cfg = [];
cfg.baselinewindow  = [-0.5 -0.05]; 
cfg.demean = 'yes';
cfg.method = 'fastica'; 
cfg.fastica.approach = 'symm'; 
cfg.fastica.interactivePCA = 'off';
cfg.fastica.g = 'gauss';
ICA_comp2 = ft_componentanalysis(cfg, pre_ICA2);
   
%Plot and Identify the bad components
nfig = ceil(length(ICA_comp2.label)./4);
topoICA = [];
topoICA.time = 1;
topoICA.label = ICA_comp2.topolabel;
topoICA.fsample = 1000;
topoICA.dimord = 'chan_time';
cfg =[];
cfg.layout = neighbour_layout.layout;
cfg.parameter = 'topography';
cfg.comment = 'no'; 
comp = 0; 
for ff = 1 :nfig
    figure
    pp = 0; 
    while pp< 4 && comp<length(ICA_comp2.label)
        pp = pp+1; 
        comp = comp+1;
        topoICA.topography = ICA_comp2.topo(:, comp);      
        subplot(4, 4, 1+4*(pp-1))  
        ft_topoplotER(cfg,topoICA);
        colormap jet
        title(['ICA comp'  num2str(comp)], 'FontSize', 18)
        timecourse = cell2mat(cellfun(@(x) x(comp,:),ICA_comp2.trial,'UniformOutput',false)');
        timecourse_pre= timecourse(:,1:490);
        fourierspctrm = fft(timecourse_pre')';
        fourierspctrm = fourierspctrm(:,1:size(timecourse_pre, 2)/2+1);
        powspctrm = (1/(topoICA.fsample*size(timecourse_pre, 2)))*abs(fourierspctrm).^2;
        powspctrm(:,2:end-1) = 2*powspctrm(:,2:end-1);
        freq =  0:topoICA.fsample/size(timecourse_pre, 2):topoICA.fsample/2;   
        subplot(4, 4, 2+4*(pp-1))
        plot(ICA_comp2.time{1},mean(timecourse,1),'LineWidth',1,'Color',[64 64 64]./255)
        xlim([-0.300 0.500])
        ylim([-250 250])
        title(['Average timecourse ICA comp'  num2str(comp)])    
        subplot(4, 4, 3+4*(pp-1))
        imagesc(ICA_comp2.time{1},1:size(timecourse,1),timecourse, [-5*mean(mean(abs(timecourse))) 5*mean(mean(abs(timecourse)))])
        colormap jet
        title(['Single trial timecourse ICA comp'  num2str(comp)])
        subplot(4, 4, 4+4*(pp-1))
        plot(freq(freq>=4 & freq<= 100),mean(powspctrm(:,freq>=4 & freq<= 100)),'LineWidth',1,'Color',[64 64 64]./255)
        title(['Power spectrum ICA comp'  num2str(comp)])
        set(gcf, 'Position', get(0, 'Screensize'));
    end
end

  %% Signal bad components from 2nd round
%Reject componments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reject_comp2= [1 4 22 26 27 40];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
cfg.component = reject_comp2;  
cfg.updatesens   = 'yes';
cfg.baselinewindow  = [-0.5 -0.1]; 
cfg.demean     = 'yes'; 
[pos_ICA2] = ft_rejectcomponent(cfg,ICA_comp2, pos_ICA1);
ICA_comp2.excluded=reject_comp2;
%% Final Processing

%Interpolate excluded channels 
missingChan = setxor(pos_ICA2.label,neighbour_layout.elec.label);
cfg = [];
cfg.method         = 'spline';
cfg.missingchannel = missingChan;
cfg.elec           = neighbour_layout.elec;
cfg.neighbours     = neighbour_layout.neighbor;
data_eeg            = ft_channelrepair(cfg,pos_ICA2);
% Put data in the original order
[~, original_order] = ismember(neighbour_layout.elec.label,data_eeg.label);
data_eeg.label = data_eeg.label(original_order);
data_eeg.trial = cellfun(@(x) x(original_order,:),data_eeg.trial,'UniformOutput',false);
data_eeg.fsample = pos_ICA2.fsample;
final_eeg_data = data_eeg;

% Average referencing
cfg = [];
cfg.reref         = 'yes';
cfg.refchannel    = {'all'};
cfg.refmethod     = 'avg';
cfg.baselinewindow  = [-0.5 -0.1]; 
cfg.demean     = 'yes'; % Remove the mean  
EEG_interp = ft_preprocessing(cfg,final_eeg_data);
    

%  PLOT TEPs 
cfg=[];
cfg.lpfilter = 'yes';
cfg.lpfreq = 45;
cfg.covariance       = 'yes';
cfg.covariancewindow = 'all';
cfg.keeptrials = 'no';
cfg.trials = find(EEG_interp.trialinfo(:,1)==0);
plot_sham = ft_timelockanalysis(cfg, EEG_interp);
cfg.trials = find(EEG_interp.trialinfo(:,1)==1);
plot_real = ft_timelockanalysis(cfg, EEG_interp);
cfg=[];
cfg.operation =  'subtract';
cfg.parameter = 'avg';
plot_difference=ft_math(cfg, plot_real, plot_sham);

subplot(3,1,1)
hold on
plot(plot_sham.time, plot_sham.avg,'Linewidth',1,'Color',[144,144,144]./255)
plot(plot_sham.time, plot_sham.avg(30,:),'Linewidth',4,'Color',[53,94,182]./255)
xlim([-0.2 0.4])
ylim([-20 20])
plot([0 0],ylim,'-k','Linewidth',2)

subplot(3,1,2)
hold on
plot(plot_real.time, plot_real.avg,'Linewidth',1,'Color',[144,144,144]./255)
plot(plot_real.time, plot_real.avg(30,:),'Linewidth',4,'Color',[53,94,182]./255)
xlim([-0.2 0.4])
ylim([-20 20])
plot([0 0],ylim,'-k','Linewidth',2)

subplot(3,1,3)
hold on
plot(plot_difference.time, plot_difference.avg,'Linewidth',1,'Color',[144,144,144]./255)
plot(plot_difference.time, plot_difference.avg(30,:),'Linewidth',4,'Color',[53,94,182]./255)
xlim([-0.2 0.4])
ylim([-20 20])
plot([0 0],ylim,'-k','Linewidth',2)


figure
j = [0.025 0.050 0.080 0.120 0.180 0.250 0.300];    
for k = 1:7
    subplot(3,7,k);
    cfg = [];   
    cfg.xlim=[j(k)-0.01 j(k)+-0.01];  
    cfg.layout = neighbour_layout;
    cfg.zlim =[-5 5];
    cfg.interactive = 'no';
    cfg.comment= 'no';
    ft_topoplotER(cfg, plot_sham );   
    colormap jet
    title([num2str(j(k)*1000) 'ms']); 
    subplot(3,7,k+7);
    ft_topoplotER(cfg, plot_real );   
    colormap jet
    subplot(3,7,k+14);
    ft_topoplotER(cfg, plot_difference );   
    colormap jet
end


%% TIME FREQUENCY RESPONSE

cfg=[];
cfg.trials=find(EEG_interp.trialinfo(:,1)==1);
active_trials=ft_selectdata(cfg,EEG_interp);
cfg=[];
cfg.trials=find(EEG_interp.trialinfo(:,1)==0);
sham_trials=ft_selectdata(cfg,EEG_interp);

cfg = [];
cfg.preproc.demean = 'yes';
cfg.keeptrials = 'yes';
ERPactive = ft_timelockanalysis(cfg, active_trials );
ERPsham = ft_timelockanalysis(cfg, sham_trials );
INDactive=ERPactive;
INDsham=ERPsham;

%subtract the time-locked evoked response from the signal prior to
%time-frequency decomposition (Premoli et al. 2017)
for tr = 1:size(ERPactive.trial,1)
    INDactive.trial(tr,:,:) = squeeze(ERPactive.trial(tr,:,:))-squeeze(mean(ERPactive.trial));
end
for tr = 1:size(ERPsham.trial,1)
    INDsham.trial(tr,:,:) = squeeze(ERPsham.trial(tr,:,:))-squeeze(mean(ERPsham.trial));
end

freq = 3:1:100;
for fr = 1:length(freq)
ncycles = 9/73*(freq(fr))+156/73;  
        cfg = [];
        cfg.polyremoval = 1;
        cfg.method     = 'wavelet';
        cfg.width      = ncycles;
        cfg.output     = 'pow';
        cfg.foi        = freq(fr);
        cfg.toi        = -1:0.01:1.5;
        cfg.keeptrials = 'yes';
        induced_active(fr) = ft_freqanalysis(cfg, INDactive);
        induced_sham(fr) = ft_freqanalysis(cfg, INDsham);
end       

% Concatenate along frequency dimension
INDUCEDactive = induced_active(1);
INDUCEDactive.powspctrm = cat(3,induced_active.powspctrm);
INDUCEDactive.freq = cat(2,induced_active.freq);
INDUCEDsham = induced_sham(1);
INDUCEDsham.powspctrm = cat(3,induced_sham.powspctrm);
INDUCEDsham.freq = cat(2,induced_sham.freq);

%Signal Standardization
ERSPfullTBz_4D = INDUCEDactive;
avgPerFreq = nanmean(INDUCEDactive.powspctrm,4); % average over time bins
avgPerFreqMat = repmat(avgPerFreq,[1 1 1 size(INDUCEDactive.powspctrm,4)]); % set average value for all time bins
sdPerFreq = nanstd(INDUCEDactive.powspctrm,0,4); % sd over time bins
sdPerFreqMat =  repmat(sdPerFreq,[1 1 1 size(INDUCEDactive.powspctrm,4)]); % set sd value for all time bins
ERSPfullTBz_4D.powspctrm = (INDUCEDactive.powspctrm - avgPerFreqMat) ./ sdPerFreqMat; % z-transform per frequency relative to full epoch
ERSPfullTBz_INDUCEDactive = ft_freqdescriptives([],ERSPfullTBz_4D); % average ERSPfullTBz over trials (4D-->3D)
clear ERSPfullTBz_4D

ERSPfullTBz_4D = INDUCEDsham;
avgPerFreq = nanmean(INDUCEDsham.powspctrm,4); % average over time bins
avgPerFreqMat = repmat(avgPerFreq,[1 1 1 size(INDUCEDsham.powspctrm,4)]); % set average value for all time bins
sdPerFreq = nanstd(INDUCEDsham.powspctrm,0,4); % sd over time bins
sdPerFreqMat =  repmat(sdPerFreq,[1 1 1 size(INDUCEDsham.powspctrm,4)]); % set sd value for all time bins
ERSPfullTBz_4D.powspctrm = (INDUCEDsham.powspctrm - avgPerFreqMat) ./ sdPerFreqMat; % z-transform per frequency relative to full epoch
ERSPfullTBz_INDUCEDsham = ft_freqdescriptives([],ERSPfullTBz_4D); % average ERSPfullTBz over trials (4D-->3D)
clear ERSPfullTBz_4D

% Remove baselines
cfg = [];
cfg.baselinetype = 'absolute';
cfg.baseline = [-0.5 -0.1];
INDUCED_active = ft_freqbaseline(cfg, ERSPfullTBz_INDUCEDactive);
INDUCED_sham = ft_freqbaseline(cfg, ERSPfullTBz_INDUCEDsham);

%%
