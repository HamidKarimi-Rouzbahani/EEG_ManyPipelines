clc;
clear all;
close all;
%% Hypothesis 3a: recognition effect in voltage amplitudes of all channels
for Subject=[1:33]
    All_channels=[1:65];
    file_name=['Z:\projects\Hamid\Projects\EEGManyPipelines\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_preprocessed_data.mat'];
    load(file_name,'Data_matrix','events','channels')
    % Separating trials into classes and performing decoding
    tr1=0;
    tr2=0;
    for trial=1:length(events)
        if strcmp(events(trial).behavior,'hit')
            tr1=tr1+1;
            ClassA(:,:,tr1)=Data_matrix(All_channels,:,trial);
        elseif strcmp(events(trial).behavior,'miss')
            tr2=tr2+1;
            ClassB(:,:,tr2)=Data_matrix(All_channels,:,trial);
        end
    end
    
    % Equalizing number of trials between classes by undersampling the
    % dominant class and repeating the procedure untill all data is used:
    % this is necessary to avoid bias in decoding
    for time=1:size(Data_matrix,2)
        ClassA_time=squeeze(ClassA(:,time,:));
        ClassB_time=squeeze(ClassB(:,time,:));
        
        folds=floor(size(ClassA_time,2)./size(ClassB_time,2));
        if size(ClassA_time,2)>size(ClassB_time,2)
            inds_tmp=1:size(ClassA_time,2);
            inds_second=size(ClassB_time,2);
        else
            inds_tmp=1:length(ClassB_time,2);
            inds_second=size(ClassA_time,2);
        end
        for fld=1:folds
            inds=randsample(inds_tmp,inds_second);
            for j=1:length(inds)
                inds_tmp((inds_tmp==inds(j)))=[];
            end
            if size(ClassA_time,2)>size(ClassB_time,2)
                Xready=[ClassA_time(:,inds)';ClassB_time'];
                Yready=[ones(size(ClassA_time(:,inds),2),1);zeros(size(ClassB_time,2),1)]';
            else
                Xready=[ClassB_time(:,inds)';ClassA_time'];
                Yready=[ones(size(ClassB_time(:,inds),2),1);zeros(size(ClassA_time,2),1)]';
            end
            Classifier_Model = fitcdiscr(Xready,Yready);
            decoding(time,fld)=1-kfoldLoss(crossval(Classifier_Model));
        end
        [Subject time]
    end
    decoding_accuracy=nanmean(decoding,2);
    file_name=['Z:\projects\Hamid\Projects\EEGManyPipelines\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_Results_Hypothesis3a.mat'];
    save(file_name,'decoding_accuracy');
    clearvars -except Subject All_channels
end

%% Plotting  and statistical testing
clc;
clear all;
close all;
% loading decoding results from all participants
Accuracies=nan(33,769);
for Subject=1:33
    file_name=['Z:\projects\Hamid\Projects\EEGManyPipelines\Analyses\Subject',num2str(Subject,'%02.f'),...
        '\','Subject',num2str(Subject,'%02.f'),'_Results_Hypothesis3a.mat'];
    load(file_name,'decoding_accuracy');
    Accuracies(Subject,:)=decoding_accuracy;
end

% Evaluation of significance of effect against decoding in the baseline
% period
for time=1:size(Accuracies,2)
    [p(time),h(time)]=signrank(Accuracies(:,time),nanmean(Accuracies(:,1:256),2),'tail','right');
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
xticks=[-256:512]*1000./512;
yticks=[45:5:60];
smoothing=5; % smoothing the decoding curve of clearer visualtion

plott=shadedErrorBar(xticks,nanmean(Accuracies*100),nanstd(Accuracies*100)*1.96./sqrt(33),{'color','b','LineWidth',2},1);
hold on;
plot([xticks(1) xticks(end)],[50 50],'k');
plot([0 0],[yticks(1) yticks(end)],'k');
Sig_plot=plot(xticks,pCorrected*47.5,'*');
set(gca,'FontSize',10,'LineWidth',1,'XTick',...
    [xticks(1):250:xticks(end)],'XTickLabel',...
    [xticks(1):250:xticks(end)],'YTick',...
    [yticks],'YTickLabel',{[yticks]},...
    'XMinorTick','on','YMinorTick','off','ycolor','k','tickdir','out','xcolor','k','box','off');
xlim([xticks(1) xticks(end)]);
ylim([yticks(1) yticks(end)]);
legend([plott.mainLine,plott.patch Sig_plot],{'Mean: hit vs. miss';'95% CI';'Significant'},'location','northwest')
xlabel('Time relative to stimulus onset (ms)');
ylabel('Decoding accuracy (%)');
title ('H3a: recognition effect in all channels');