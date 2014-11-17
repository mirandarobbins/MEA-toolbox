function spikes=MUA64plot(data,params,spikes)
%% plot all spikes//all reps//all channels
spikes.no_trials=size(spikes.spike_locs,2);
spikes.no_spikes=cellfun(@numel,spikes.spiketimes);
spikes.total_spikes=sum(sum(spikes.no_spikes));
params.flags.plot_online=0;
if params.flags.plot_online
figure; 
for channel_id=1:size(spikes.spiketimes,1)
    clear temp;hold on
%     temp=cell2mat(spikes.spiketimes(channel_id,:));
%     scatter(data{1,channel_id}*BinWidth,repmat(channel_id,1,no_events),'x
%     k')
% scatter(data{1,channel_id},repmat(channel_id,1,no_events),'xk')
temp=horzcat(spikes.spiketimes{channel_id,1:spikes.no_trials});temp=temp*10;
scatter(temp,repmat(channel_id,1,numel(temp)),'. k')
end; clear channel_id temp
axis([0 1000 0 size(spikes.spiketimes,1) ])
title('Plotting all spikes//all reps//all channels'); xlabel('Time (ms)'); ylabel('Channel no.')
else
end
%% rejiggle, convert to stamps....make a PDF
spikes.PDF=[];
spikes.PDF.data_aligned=spikes.spiketimes';
for c_idx=1:size(spikes.PDF.data_aligned,1);
    for r_idx=1:size(spikes.PDF.data_aligned,2);
        if isempty(spikes.PDF.data_aligned{c_idx,r_idx})
            spikes.PDF.data_aligned{c_idx,r_idx}(1)=NaN;
        end
    end
end
spikes.PDF.time_window=[0 1000];
spikes.PDF.binwidth=1;
% spikes.PDF.spiketimes=cell(1,size(spikes.PDF.data_aligned,2));
spikes.PDF.spiketimes=cell(1,params.last_sweep);
for channel_id=1:size(spikes.PDF.data_aligned,2)
    spikes.PDF.spiketimes{channel_id}=zeros(size(spikes.PDF.data_aligned,1),1); % need to force first with  AT LEAST 1 zero
    max_no_spikes=max(max(cell2mat(cellfun(@size,spikes.PDF.data_aligned(:,channel_id),'UniformOutput',0))));
    spikes.PDF.spiketimes{channel_id}=zeros(size(spikes.PDF.data_aligned,1),max_no_spikes); % +1 to catch cases with no spikes
    
    spikes.PDF.timestamps{channel_id}=zeros(params.last_sweep,spikes.PDF.time_window(2));
    for trial_id=data.sweep_sort.successful_sweeps
         spikes.PDF.spiketimes{channel_id}(trial_id,1:numel(spikes.PDF.data_aligned{trial_id,channel_id}))=...
            spikes.PDF.data_aligned{trial_id,channel_id}*10;
         spikes.PDF.spiketimes{channel_id}(spikes.PDF.spiketimes{channel_id}==0)=NaN;
        spikes.PDF.timestamps{channel_id}(trial_id,:)=spike_counts(spikes.PDF.spiketimes{channel_id}(trial_id,:)',...
                                               spikes.PDF.time_window,...
                                               spikes.PDF.binwidth);
    end
     
end
%% evaluate the kernel - NB this normalises only for trials that show significant LFP activity in search win
spikes.PDF.kernel_sigma =0.001;% .01; % finer detail =0.005
spikes.PDF.kernel_shoulder=3;
spikes.PDF.edges=-spikes.PDF.kernel_shoulder:0.001:spikes.PDF.kernel_shoulder;
spikes.PDF.kernel=normpdf(spikes.PDF.edges,0,spikes.PDF.kernel_sigma);

% Multiply by bin width so the probabilities sum to 1
spikes.PDF.kernel=spikes.PDF.kernel*spikes.PDF.binwidth*1E-3; 
% Find the index of the kernel center
spikes.PDF.kernel_center =ceil(length(spikes.PDF.edges)/2); 

figure; 
%Convolve time-stamped spike data with the kernel    
spikes.PDF.PDF_trimmed=[];
for channel_id=1:size(spikes.PDF.data_aligned,2)
%     for trial_id=1:size(spikes.PDF.data_aligned,1); % to include all trials
    for trial_id=(data.sweep_sort.successful_sweeps); % to include successful trials only
        spikes.PDF.PDF_trimmed{channel_id}(trial_id,:)=...
            conv(spikes.PDF.timestamps{channel_id}(trial_id,:),spikes.PDF.kernel);
        % Trim out the relevant portion of the spike density result        
    end
     % trim outliers
         spikes.PDF.PDF_trimmed{channel_id}=...
            spikes.PDF.PDF_trimmed{channel_id}(:,...
                spikes.PDF.kernel_center:spikes.PDF.time_window(2)+spikes.PDF.kernel_center-1); 
    spikes.PDF.PDF_trimmed{channel_id}(spikes.PDF.PDF_trimmed{channel_id}==0)=NaN;  %<--- Un-comment this line to discount silent regions

    spikes.PDF.mean_PDF{channel_id}=nanmean(spikes.PDF.PDF_trimmed{channel_id});
    spikes.PDF.std_PDF{channel_id}=nanstd(spikes.PDF.PDF_trimmed{channel_id});
    plot(spikes.PDF.mean_PDF{channel_id}); hold on
end
spikes.PDF.mean_PDF=cell2mat(spikes.PDF.mean_PDF');
spikes.PDF.std_PDF=cell2mat(spikes.PDF.std_PDF');
% Grand mean
spikes.PDF.grand_mean_PDF=nanmean(spikes.PDF.mean_PDF,1);
spikes.PDF.grand_std_PDF=nanstd(spikes.PDF.mean_PDF,1);
%% plot rasters and SDFs
channels_to_plot=1:size(spikes.PDF.mean_PDF,1);
% channels_to_plot=[3 11 19];
c_map=hot(size(spikes.PDF.mean_PDF,1));
c_map(1:8,:)  =repmat([1 0 0],8,1);
c_map(9:16,:) =repmat([1 0 0],8,1);
c_map(17:24,:)=repmat([0.5 1 0],8,1);
c_map(25:32,:)=repmat([0.5 1 0],8,1);
c_map(33:40,:)=repmat([0 0.5 1],8,1);
c_map(41:48,:)=repmat([0 0.5 1],8,1);
c_map(49:56,:)=repmat([0 0 1],8,1);
c_map(57:64,:)=repmat([0 0 1],8,1);

h.f1=figure; 
h.f1_a=subplot(3,1,1);    
hold on
    for channel_id=channels_to_plot             
        
        temp=horzcat(spikes.spiketimes{channel_id,1:spikes.no_trials})*10;
        scatter(temp,repmat(channel_id,1,numel(temp)),15,c_map(channel_id,:),'.')
    end
    clear time_win channel_id trial_id temp
    set(h.f1_a,'Xtick',[])
    ylabel('Channel no.')
    axis([0 1000 0 size(spikes.PDF.mean_PDF,1)])
h.f1_b=subplot(3,1,2);   
%     plot(spikes.PDF.mean_PDF');
for channel_id=channels_to_plot   
    

 hold on
%     ciplot(spikes.PDF.mean_PDF(channel_id,:)-spikes.PDF.std_PDF(channel_id,:),...
%            spikes.PDF.mean_PDF(channel_id,:)+spikes.PDF.std_PDF(channel_id,:),...
%            (spikes.PDF.time_window(1)+1:spikes.PDF.time_window(2)),c_map(channel_id,:))
       
    plot((spikes.PDF.time_window(1)+1:spikes.PDF.time_window(2)),spikes.PDF.mean_PDF(channel_id,:),'color',c_map(channel_id,:),'LineWidth',1.5)
    axis([0 spikes.PDF.time_window(2) 0 Inf])
    
    
    set(h.f1_b,'Xtick',[])
    ylabel('spike pobability density')
end
h.f1_c=subplot(3,1,3);
    hold on
    ciplot(spikes.PDF.grand_mean_PDF-spikes.PDF.grand_std_PDF,...
           spikes.PDF.grand_mean_PDF+spikes.PDF.grand_std_PDF,...
           (spikes.PDF.time_window(1)+1:spikes.PDF.time_window(2)),'b')
    plot((spikes.PDF.time_window(1)+1:spikes.PDF.time_window(2)),spikes.PDF.grand_mean_PDF,'b','LineWidth',1.5)
    axis([0 spikes.PDF.time_window(2) 0 Inf])
    xlabel('post-stimulus time (ms)')
    ylabel('spike pobability density')

%% Sort by max activity
% 
% for chan_id=1:size(spikes.PDF.mean_PDF,1)
%     first_spike(chan_id)=find(spikes.PDF.mean_PDF((chan_id),:),1); % find first response
% end
% 
% [R spikes.PDF.chan_order]= sort(max(spikes.PDF.mean_PDF(:,100:200),[],2)); %Sort by most activity
% [R spikes.PDF.time_order]= sort(first_spike); %sort by first response
% clear R first spike
% spikes.PDF.activity_sorted=zeros(size(spikes.PDF.mean_PDF));
% spikes.PDF.time_sorted=zeros(size(spikes.PDF.mean_PDF));
%     for i=1:size(spikes.PDF.mean_PDF,1)
%           spikes.PDF.activity_sorted(i,:)   =spikes.PDF.mean_PDF(spikes.PDF.chan_order==i,:); %Sort by most activity
%           spikes.PDF.time_sorted(i,:)       =spikes.PDF.mean_PDF(spikes.PDF.time_order==i,:); %sort by first response
%     end
%     figure;
%     subplot(1,2,1);imagesc(spikes.PDF.activity_sorted); colormap(flipud(gray))
%     subplot(1,2,2);imagesc(spikes.PDF.time_sorted); colormap(flipud(gray))
%     
    
%% plot mean PDF per channel
figure;
 for plot_id=1:params.no_channels

        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
%          hold on; box off; axis off
         plot(spikes.PDF.mean_PDF(plot_id,:),'LineWidth',1.2)
         axis([0 1000 0 1.1*max(max(spikes.PDF.mean_PDF))])
         
 end
 
 
 
 
 
 
 %% plot one column or row of electrodes
spikes.PDF.mean_PDF(isnan(spikes.PDF.mean_PDF))=0;
column_to_plot=4;
channels_to_plot=params.channel_index(:,column_to_plot);
trial_to_plot=data.sweep_sort.successful_sweeps;
% trial_to_plot=params.selected_rep;

figure; 
for plot_id=1:8        
    subplot(8,1,plot_id)%, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.1);
    hold on
    plot(data.tb,squeeze(data.filtered_spikes(:,channels_to_plot(plot_id),trial_to_plot)),'k')
        plot(spikes.PDF.mean_PDF(channels_to_plot(plot_id),:),'r','LineWidth',1.5)
%         plot(spikes.PDF.PDF_trimmed{channels_to_plot(plot_id)}(trial_to_plot,:)','r','LineWidth',0.5) % for including failures
        plot(spikes.PDF.PDF_trimmed{channels_to_plot(plot_id)}','r','LineWidth',0.5) % for trimmed failures

        text(800,0.07,strcat('Chan. ',num2str(channels_to_plot(plot_id))),'FontWeight','bold')

    axis([0 1000 -0.1 0.1])
    set(gca,'xtick',[])
    set(gca,'ytick',[-.1 .1])
    ylabel('Spike Density')
end
    set(gca,'xtick',[0:100:1000])

xlabel({'Time (ms)' ;''; strcat('Multi-unit activity on MEA column ', num2str(column_to_plot))})

%% plot MUA and spke desnity functions for whole array
spikes.PDF.mean_PDF(isnan(spikes.PDF.mean_PDF))=0;

figure; 
for plot_id=1:64        
    channels_to_plot=params.channel_index(plot_id);
    
        subaxis(8,8,plot_id, 'Spacing', 0.01, 'Padding', 0, 'Margin', 0.01);
    hold on
    plot(data.tb,squeeze(data.filtered_spikes(:,plot_id,trial_to_plot))*3,'k')
        plot(spikes.PDF.mean_PDF(plot_id,:)./2,'Color',[0.3 0.3 0.3],'LineWidth',1.5)

%         text(800,0.07,strcat('Chan. ',num2str(channels_to_plot(plot_id))),'FontWeight','bold')

    axis([200 400 -0.01 0.3])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
%     ylabel('Spike Density')
end
%     set(gca,'xtick',[0:100:1000])

% xlabel({'Time (ms)' ;''; strcat('Multi-unit activity on MEA column ', num2str(column_to_plot))})

assignin('base', 'spikes', spikes)



