%rename AD 'data' as 'temp'


no_trials=size(temp.raw_data,3);

fileDescription=cell(1,no_trials);
data=cell(1,no_trials);
for trial_id=1:no_trials;
    data{trial_id}=temp.raw_data(1:20000,:,trial_id)';
    fileDescription{trial_id}=strcat('trial_',num2str(trial_id));
end

BinWidth=1/20000;
chanInfo=[1:64];
plotOption = struct('raw',       1,...
                    'denoise',   0,...
                    'detection', 0,...
                    'stimulus',  0,...
                    'LFP',       0,...
                    'spiketrain',0,...
                    'Trials',    0....
                    );
                
clear temp trial_id 

%% blank stim artifacts

for trial_id=1:no_trials
    data{1,trial_id}(:,[2000:2020,2400:2420,2800:2820,3200:3220,3600:3620])=0; %blank 5x50Hz
end
    
    
clear no_trials trial_id

%% plot spike times from loaded spiketimes
no_trials=length(LFP);
figure; hold on
for neuron_id=1:size(data,2)
    no_events=size(data{1,neuron_id},2);
%     scatter(data{1,neuron_id}*BinWidth,repmat(neuron_id,1,no_events),'xk')
scatter(data{1,neuron_id},repmat(neuron_id,1,no_events),'xk')
end
%% align all trials to stim train
% time look up table
clear data_aligned
time_win=zeros(no_trials,2);time_win(:,1)=1:no_trials;
time_win(2:no_trials,1)=(time_win(2:no_trials,1)-1)*(lengthData/no_trials)+1;
time_win(:,2)=(time_win(:,1)+lengthData/no_trials-1);

figure; hold on
for neuron_id=1:length(data)
    for trial_id=1:length(LFP);
        temp=data{neuron_id};
        temp(temp<time_win(trial_id,1))=[]; temp(temp>time_win(trial_id,2))=[];
        temp=temp-time_win(trial_id,1)+1;
        temp=temp*BinWidth;
        data_aligned{trial_id,neuron_id}=temp;
        
        scatter(temp,repmat(neuron_id,1,numel(temp))) 
%         clear temp
    end
end
% clear time_win neuron_id trial_ids
%% rejiggle, convert to stamps....make a PDF
spikes.time_window=[0 1000];
spikes.binwidth=1;
for neuron_id=1:size(data_aligned,2)
    max_no_spikes=max(max(cell2mat(cellfun(@size,data_aligned(:,neuron_id),'UniformOutput',0))));
    spikes.spiketimes{neuron_id}=zeros(size(data_aligned,1),max_no_spikes);
     for trial_id=1:size(data_aligned,1);
         spikes.spiketimes{neuron_id}(trial_id,1:numel(data_aligned{trial_id,neuron_id}))=data_aligned{trial_id,neuron_id}*1000;
         spikes.spiketimes{neuron_id}(spikes.spiketimes{neuron_id}==0)=NaN;
     end
     spikes.timestamps{neuron_id}=spike_counts(spikes.spiketimes{neuron_id},...
                                               spikes.time_window,...
                                               spikes.binwidth);
end
% evaluate the kernel
spikes.kernel_sigma = .005; 
spikes.kernel_shoulder=3;
spikes.edges=-spikes.kernel_shoulder:0.001:spikes.kernel_shoulder;
spikes.kernel=normpdf(spikes.edges,0,spikes.kernel_sigma);

% Multiply by bin width so the probabilities sum to 1
spikes.kernel=spikes.kernel*spikes.binwidth*1E-3; 
% Find the index of the kernel center
spikes.kernel_center =ceil(length(spikes.edges)/2); 

figure; 
%Convolve time-stamped spike data with the kernel    
spikes.PDF_trimmed=[];
for neuron_id=1:size(data_aligned,2)
    for trial_id=1:size(data_aligned,1);
        spikes.PDF_trimmed{neuron_id}(trial_id,:)=...
            conv(spikes.timestamps{neuron_id}(trial_id,:),spikes.kernel);
        % Trim out the relevant portion of the spike density result        
    end
     % trim outliers
         spikes.PDF_trimmed{neuron_id}=...
            spikes.PDF_trimmed{neuron_id}(:,...
                spikes.kernel_center:spikes.time_window(2)+spikes.kernel_center-1); 
    spikes.PDF_trimmed{neuron_id}(spikes.PDF_trimmed{neuron_id}==0)=NaN;  %<--- Un-comment this line to discount silent regions

    spikes.mean_PDF{neuron_id}=nanmean(spikes.PDF_trimmed{neuron_id});
    spikes.std_PDF{neuron_id}=nanstd(spikes.PDF_trimmed{neuron_id});
    plot(spikes.mean_PDF{neuron_id}); hold on
end
spikes.mean_PDF=cell2mat(spikes.mean_PDF');
spikes.std_PDF=cell2mat(spikes.std_PDF');
%% Grand mean
spikes.grand_mean_PDF=nanmean(spikes.mean_PDF,1);
spikes.grand_std_PDF=nanstd(spikes.mean_PDF,1);


time_win=zeros(no_trials,2);time_win(:,1)=1:no_trials;
time_win(2:no_trials,1)=(time_win(2:no_trials,1)-1)*(lengthData/no_trials)+1;
time_win(:,2)=(time_win(:,1)+lengthData/no_trials-1);

h.f1=figure; 
h.f1_a=subplot(3,1,1);    
    for neuron_id=1:length(data)
        for trial_id=1:length(LFP);
            temp=data_aligned{trial_id,neuron_id};
            hold on
            scatter(data_aligned{trial_id,neuron_id}*1000,repmat(neuron_id,1,numel(temp)),'.k')
        end
    end
    clear time_win neuron_id trial_id temp
    set(h.f1_a,'Xtick',[])
    ylabel('Neuron number')
h.f1_b=subplot(3,1,2);
    plot(spikes.mean_PDF');
    set(h.f1_b,'Xtick',[])
    ylabel('spike pobability density')

h.f1_c=subplot(3,1,3);
    hold on
    ciplot(spikes.grand_mean_PDF-spikes.grand_std_PDF,...
           spikes.grand_mean_PDF+spikes.grand_std_PDF,...
           (spikes.time_window(1)+1:spikes.time_window(2))-100,'b')
    plot((spikes.time_window(1)+1:spikes.time_window(2))-100,spikes.grand_mean_PDF,'b','LineWidth',1.5)
    axis([-100 spikes.time_window(2)-100 0 Inf])
    xlabel('post-stimulus time (ms)')
    ylabel('spike pobability density')