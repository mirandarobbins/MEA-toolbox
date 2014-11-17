fifty.meta=meta;
twenty.meta=meta;

%% Plot WT vs KO twenty
figure
subplot(2,2,1); hold on
% plot(twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,'b')
% plot(twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,'r')
plot(twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean,'b')
plot(twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean,'r')
ciplot(twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean+twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean-twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       1:1000,'b')
ciplot(twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean+twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean-twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       1:1000,'r')
title('Mean LFP activity at burst centre')
xlabel('Time (ms)')   
ylabel('LFP amplitude (mV)')   
axis([0 1000 -0.1 0.1])
subplot(2,2,2); hold on
% plot(fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,'b')
% plot(fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_pop_mean,'r')
plot(fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean,'b')
plot(fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean,'r')
ciplot(fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean+fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean-fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       1:1000,'b')
ciplot(fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean+fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_mean-fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_LFP_grand_SEM,...
       1:1000,'r')
title('Mean LFP activity at burst centre')
xlabel('Time (ms)')   
ylabel('LFP amplitude (mV)')  
axis([0 1000 -0.1 0.1])


subplot(2,2,3); hold on
% plot(twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,'b')
% plot(twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,'r')
plot(twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean,'b')
plot(twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean,'r')
upper=twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean+...
      twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower=twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean-...
      twenty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower(isnan(lower))=0; upper(isnan(upper))=0;

ciplot(lower,upper,1:1000,'b')

upper=twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean+...
      twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
  
lower=twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean-...
      twenty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower(isnan(lower))=0; upper(isnan(upper))=0;
ciplot(lower,upper,1:1000,'r')
title('Mean spike density at burst centre')
xlabel('Time (ms)')   
ylabel('Spike density')  
axis([0 1000 0 1])

subplot(2,2,4); hold on
% plot(fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,'b')
% plot(fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_pop_mean,'r')
plot(fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean,'b')
plot(fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean,'r')
upper=fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean+...
      fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower=fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean-...
      fifty.meta.WT_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower(isnan(lower))=0; upper(isnan(upper))=0;

ciplot(lower,upper,1:1000,'b')

upper=fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean+...
      fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
  
lower=fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_mean-...
      fifty.meta.KO_exported.mean_responses.strongest_channel.ByStrongest_PDF_grand_SEM;
lower(isnan(lower))=0; upper(isnan(upper))=0;
ciplot(lower,upper,1:1000,'r')
title('Mean spike density at burst centre')
xlabel('Time (ms)')   
ylabel('Spike density')   
axis([0 1000 0 1])

