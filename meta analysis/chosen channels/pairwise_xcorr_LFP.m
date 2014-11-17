function x_corr=pairwise_xcorr_LFP(L4,L23)
    x_corr.params.windowsize=100;
    x_corr.params.windowsize_ms=x_corr.params.windowsize*20000^-1*1000;
    x_corr.params.noverlap=20;
    x_corr.params.noverlap_ms=x_corr.params.noverlap*20000^-1*1000;
    x_corr.params.maxlags=500;                                                % samples
    x_corr.params.lags=(-x_corr.params.maxlags:x_corr.params.maxlags)';       % samples
    x_corr.params.maxlags_ms=x_corr.params.maxlags*20000^-1*1000;     % ms
    x_corr.params.lags_ms=(x_corr.params.lags)*20000^-1*1000;         % ms

t1=tic;

for sweep_id=1:size(L4,2)
    %%%%% Do X-corr
    signal_1=L4(:,sweep_id); 
    signal_2=L23(:,sweep_id); 
        signal_1([1994:2040, 2394:2440,  2794:2840, 3194:3240, 3594:3640]) = 0;
        signal_2([1994:2040, 2394:2440,  2794:2840, 3194:3240, 3594:3640]) = 0;
        signal_1=zscore(signal_1);signal_2=zscore(signal_2);
        [x_corr.c{sweep_id} ...
         x_corr.L...
         x_corr.T ]     =   corrgram(signal_1,signal_2,...
                                     x_corr.params.maxlags,...
                                     x_corr.params.windowsize,...
                                     x_corr.params.noverlap);        
        x_corr.c{sweep_id}=smooth2a(x_corr.c{sweep_id},20,0); %2D smoothing ...,lags , time
        x_corr.max_c{sweep_id}=zeros(1,size(x_corr.c{sweep_id},2));
        x_corr.max_c_lagindex{sweep_id}=zeros(2,size(x_corr.c{sweep_id},2));
        x_corr.actual_max_lag{sweep_id}=zeros(1,size(x_corr.c{sweep_id},2));
  
    
end
clear signal_2 signal_1 progbar channel_id channel_index chosen_channel t1 this_sweep sweep_id
%% Mean x-correlogram
x_corr.mean_c=[];
temp=[]; temp2=[];
% mean by trial
for sweep_id=1:size(L4,2)
    
    temp=x_corr.c(:,sweep_id);
    temp(cellfun(@numel,temp)==0)=[];
    for idx=1:numel(temp)
        test=temp{idx}; test(isnan(test))=0;
        if max(var(test))==0        
            test=zeros(size(test)); test(test==0)=NaN;
        end
            temp2(:,:,idx)=temp{idx};
    end; clear idx
    x_corr.mean_c.mean_c_by_sweep(:,:,sweep_id)=nanmean(temp2,3);
end
x_corr.mean_c.mean_c_by_sweep(isnan(x_corr.mean_c.mean_c_by_sweep))=0;
x_corr.mean_c.mean_c_grand = nanmean(x_corr.mean_c.mean_c_by_sweep,3);
x_corr.mean_c.SEM_c_grand = nansem(x_corr.mean_c.mean_c_by_sweep,3);
% x_corr.mean_c.mean_c_grand=smooth2a(x_corr.mean_c.mean_c_grand,0,10); %2D smoothing ...,lags , time
clear temp temp2 sweep_id plotYN chosen_channel this_sweep