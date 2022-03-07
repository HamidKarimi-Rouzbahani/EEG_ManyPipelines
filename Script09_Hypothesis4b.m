clc;
clear all;
close all;
%% Hypothesis 4b: subsequent memory effect in power spectrum
for Subject=[1:33]
    All_channels=[1:65];
    Fs=512; % Sampling freq
    file_name=['Z:\projects\Hamid\Projects\EEGManyPipelines\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_preprocessed_data.mat'];
    load(file_name,'Data_matrix','events','channels')
    % Decoding
    tr1=0;
    tr2=0;
    for trial=1:length(events)
        if strcmp(events(trial).subsequent_memory,'subsequent_remembered')
            tr1=tr1+1;
            ClassA(:,:,tr1)=Data_matrix(All_channels,:,trial);
        elseif strcmp(events(trial).subsequent_memory,'subsequent_forgotten')
            tr2=tr2+1;
            ClassB(:,:,tr2)=Data_matrix(All_channels,:,trial);
        end
    end
    % Extracting time-frequency power (Event-related spectral purturbation
    % as defined in EEGLAB): for every trial for decoding
    epoch_wind=[-256:512];
    erspA=nan(65,7500,size(ClassA,3));
    erspB=nan(65,7500,size(ClassB,3));
    for ch=1:65
        for trial=1:size(ClassA,3)
            [tmpA,~,~,times,freqs,~,~]=newtimef(ClassA(ch,:,trial),size(Data_matrix,2),...
                [epoch_wind(1) epoch_wind(end)]*1000/Fs,Fs, 0,...
                'timesout',[150],'freqs',[0 200],'nfreqs',50,...
                'plotersp','off','plotitc','off');
            erspA(ch,:,trial)=reshape(tmpA,1,[]);
        end
        for trial=1:size(ClassB,3)
            [tmpB,~,~,times,freqs,~,~]=newtimef(ClassB(ch,:,trial),size(Data_matrix,2),...
                [epoch_wind(1) epoch_wind(end)]*1000/Fs,Fs, 0,...
                'timesout',[150],'freqs',[0 200],'nfreqs',50,...
                'plotersp','off','plotitc','off');
            erspB(ch,:,trial)=reshape(tmpB,1,[]);
        end
    end
    % Equalizing number of trials between classes by undersampling the
    % dominant class and repeating the procedure untill all data is used:
    % this is necessary to avoid bias in decoding
    for time_freq=1:size(erspA,2)
        ClassA_time=squeeze(erspA(:,time_freq,:));
        ClassB_time=squeeze(erspB(:,time_freq,:));
        if size(ClassA_time,2)>size(ClassB_time,2)
            inds_tmp=1:size(ClassA_time,2);
            inds_second=size(ClassB_time,2);
            folds=floor(size(ClassA_time,2)./size(ClassB_time,2));
        else
            inds_tmp=1:size(ClassB_time,2);
            inds_second=size(ClassA_time,2);
            folds=floor(size(ClassB_time,2)./size(ClassA_time,2));
        end
        for fld=1:folds
            if size(ClassA_time,2)~=size(ClassB_time,2)
                inds=randsample(inds_tmp,inds_second);
                for j=1:length(inds)
                    inds_tmp((inds_tmp==inds(j)))=[];
                end
            else
                inds=inds_tmp;
            end
            if size(ClassA_time,2)>size(ClassB_time,2)
                Xready=[ClassA_time(:,inds)';ClassB_time'];
                Yready=[ones(size(ClassA_time(:,inds),2),1);zeros(size(ClassB_time,2),1)]';
            else
                Xready=[ClassB_time(:,inds)';ClassA_time'];
                Yready=[ones(size(ClassB_time(:,inds),2),1);zeros(size(ClassA_time,2),1)]';
            end
            Classifier_Model = fitcdiscr(Xready,Yready);
            decoding(time_freq,fld)=1-kfoldLoss(crossval(Classifier_Model));
        end
        [Subject time_freq]
    end
    decoding_accuracy=nanmean(decoding,2);
    
    file_name=['Z:\projects\Hamid\Projects\EEGManyPipelines\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_Results_Hypothesis4b.mat'];
    save(file_name,'decoding_accuracy','times','freqs');
    clearvars -except Subject All_channels
end
%% Plotting  and statistical testing
clc;
clear all;
close all;
% loading decoding results from all participants
Accuracies=nan(33,7500);
s=0;
for Subject=[1:33]
    s=s+1;
    file_name=['Z:\projects\Hamid\Projects\EEGManyPipelines\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_Results_Hypothesis4b.mat'];
    load(file_name);
    Accuracies(s,:)=decoding_accuracy;
    AccuraciesT(:,:,s)=reshape(Accuracies(s,:),[50 150]);
end

% Evaluation of significance of effect against decoding in the baseline
% period
for freq=1:size(AccuraciesT,1)
    for time_freq=1:size(AccuraciesT,2)
        [p(freq,time_freq),h(freq,time_freq)]=signrank(squeeze(AccuraciesT(freq,time_freq,:)),squeeze(nanmean(AccuraciesT(freq,1:45,:),2)),'tail','right');
    end
end


% Corrrecting p values for multiple comparisons using Combined Probability of Fisher
threshold=0.05;
[pCorrected,~, ~, ~] = mt_fisher(reshape(p,1,[]), threshold);
pnew=zeros(size(pCorrected));
pnew(pCorrected<threshold)=1;
time_span_of_clustering= 5; % 5 = 10 ms: only the time and frequency spans
% in which all points are sginificant will be marked as significant
pCorrected=bwareaopen(pnew,time_span_of_clustering);

% plotting decoding curves and indicating time points of significant effect
AccuraciesT=squeeze(nanmean(AccuraciesT,3));
limits_scale=round([min(min(AccuraciesT)) max(max(AccuraciesT))]*100);
imagesc(AccuraciesT);
AccuraciesT=uint8(round(rescale(AccuraciesT).*64));
ERSP=flipud(AccuraciesT);
image(ERSP)
hold on
pCorrected=reshape(pCorrected,50,150);
imcontour(flipud(pCorrected),1,'k')


line([46 46],[1 size(ERSP,1)],'color','k');

yticks=flip(freqs,2);

set(gca,'FontSize',10,'LineWidth',1,'XTick',...
    [16 46 75 105 135],'XTickLabel',...
    [-250:250:750],'YTick',...
    [1:5:size(ERSP,1)],'YTickLabel',{downsample( yticks , 5 )},...
    'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k','box','off');
xlim([1 size(ERSP,2)]);
ylim([1 size(ERSP,1)]);
xlabel('Time relative to stimulus onset (ms)');
ylabel('Frequency (Hz)');
title ('H4b: subsequent memory effect in spectral power');

c=colorbar('Ticks',[linspace(0,64,5)],...
    'TickLabels',{[linspace(limits_scale(1),limits_scale(end),5)]});
c.Label.String = 'Decoding accuracy (%)';
set(c,'FontSize',10,'LineWidth',1)