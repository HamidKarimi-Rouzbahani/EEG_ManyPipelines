clc;
clear all;
close all;

for Subject=[1:33]
    pathToFile=['Z:\projects\Hamid\Projects\EEGManyPipelines\EMP_data\EMP_data\eeg_eeglab\'];
    Loaded_set = pop_loadset(['EMP',num2str(Subject,'%02.f'),'.set'], pathToFile);
    %% EEGLAB Band-Pass filtering
    lowband=0.05;
    highband=200;
    Preprocessed_data = pop_eegfiltnew(Loaded_set, lowband, highband, [],0,[],0,0,0);
    %% EEGLAB Band-Stop (notch) filtering
    lowband=45;
    highband=55;
    Preprocessed_data = pop_eegfiltnew(Preprocessed_data, lowband, highband, [],1,[],0,0,0);
    %% EEGLAB Theta band filtering
    lowband=4;
    highband=8;
    Preprocessed_data_theta = pop_eegfiltnew(Loaded_set, lowband, highband, [],0,[],0,0,0);
    %% EEGLAB Alpha band filtering
    lowband=8;
    highband=12;
    Preprocessed_data_alpha = pop_eegfiltnew(Loaded_set, lowband, highband, [],0,[],0,0,0); 
    clearvars Loaded_set
    %% Epoching
    trial=0;
    epoch_wind=[-256:512];
    EEG_channels_to_keep=[1:29 31:70];
    Preprocessed_data.data=Preprocessed_data.data(EEG_channels_to_keep,:);
    Preprocessed_data_theta.data=Preprocessed_data_theta.data(EEG_channels_to_keep,:);
    Preprocessed_data_alpha.data=Preprocessed_data_alpha.data(EEG_channels_to_keep,:);
    Data_matrix=nan(length(EEG_channels_to_keep),length(epoch_wind),1200);
    Data_matrix_theta=nan(length(EEG_channels_to_keep),length(epoch_wind),1200);
    Data_matrix_alpha=nan(length(EEG_channels_to_keep),length(epoch_wind),1200);
    for event=1:length(Preprocessed_data.event)
        trial=trial+1;
        trial_data=Preprocessed_data.data(:,Preprocessed_data.event(1,event).latency+epoch_wind(1):Preprocessed_data.event(1,event).latency+epoch_wind(end));
        Data_matrix(:,1:length(epoch_wind),trial)=trial_data-repmat(nanmean(trial_data(:,1:256),2),[1 length(epoch_wind)]);
        
        trial_data_theta=Preprocessed_data_theta.data(:,Preprocessed_data_theta.event(1,event).latency+epoch_wind(1):Preprocessed_data_theta.event(1,event).latency+epoch_wind(end));
        Data_matrix_theta(:,1:length(epoch_wind),trial)=trial_data_theta-repmat(nanmean(trial_data_theta(:,1:256),2),[1 length(epoch_wind)]);

        trial_data_alpha=Preprocessed_data_alpha.data(:,Preprocessed_data_alpha.event(1,event).latency+epoch_wind(1):Preprocessed_data_alpha.event(1,event).latency+epoch_wind(end));
        Data_matrix_alpha(:,1:length(epoch_wind),trial)=trial_data_alpha-repmat(nanmean(trial_data_alpha(:,1:256),2),[1 length(epoch_wind)]);
    end
    events=Preprocessed_data.event;
    channels=Preprocessed_data.chanlocs(EEG_channels_to_keep);
    clearvars Preprocessed_data Preprocessed_data_theta Preprocessed_data_alpha
    %% Saving
    mkdir(['Subject',num2str(Subject,'%02.f')]);
    file_name=['Z:\projects\Hamid\Projects\EEGManyLabs\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_preprocessed_data.mat'];
    save(file_name,'Data_matrix','Data_matrix_theta','Data_matrix_alpha','events','channels');
    clearvars -except Subject
    [Subject]
    pause(120);
end