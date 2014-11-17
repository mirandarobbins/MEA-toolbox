% temp=data.filtered_lfp(:,:,1);
temp=CSD.csd_array(:,:,1);
temp=downsample(temp,20)';
temp=smooth(temp);
temp2=princomp(temp);
temp3=temp*temp2;
plot(temp')
figure;
plot(temp2(1,:),temp2(2,:)) 

%%
input=downsample(data.filtered_lfp(:,:,1),20);
f=input;%reshape(input(:,:,1),size(input,1)*size(input,2),1);
% mu=mean(f);
% f_hat=f-mu;
% sigma=f_hat*f';
% [EigValue EigVectors] = eig(sigma);
%%
figure; hold on;
for trial_id=data.sweep_sort.successful_sweeps
input=downsample(data.filtered_lfp(200:1000,:,trial_id),20);
f=input;
f_hat=input-mean(mean(f));
[COEFF,SCORE]=princomp(f);
VT=SCORE(:,1:3)';
alpha=VT*f_hat;
plot3(alpha(1,:),alpha(2,:),alpha(3,:))
end
% plot3(SCORE(:,1),SCORE(:,2),SCORE(:,3))

%%
time_lims=200:1000;
n_time_lims=numel(time_lims);
% temp=downsample(data.filtered_lfp(:,:,data.sweep_sort.successful_sweeps(1)),2);
temp=downsample(CSD.csd_array(:,:,data.sweep_sort.successful_sweeps(1)),20); input=temp(time_lims,:);

for trial_id=data.sweep_sort.successful_sweeps(2:numel(data.sweep_sort.successful_sweeps))
temp=downsample(CSD.csd_array(:,:,trial_id),20);
input2=temp(time_lims,:);
input=vertcat(input,input2);
end; clear temp input2
[COEFF,SCORE]=princomp(input);
% SCORE=smooth2a(SCORE,1,1);
figure; hold on
for trial_id=1:size(input,1)/numel(time_lims)
% view(3)
% shading interp
% camlight
% grid on
plot3(SCORE(1:n_time_lims,1),...
      SCORE(1:n_time_lims,2),...
      SCORE(1:n_time_lims,3),'b','LineWidth',1.5)
SCORE(1:n_time_lims,:)=[];
end;

grid on;  

%%
time_lims=1:1000;
n_time_lims=numel(time_lims);



temp=spikes.PDF.mean_PDF';
input=temp(time_lims,:);

clear temp input2
[COEFF,SCORE]=princomp(input);
% SCORE=smooth2a(SCORE,1,1);
figure; hold on
for trial_id=1
% view(3)
% shading interp
% camlight
% grid on
plot3(SCORE(1:n_time_lims,1),...
      SCORE(1:n_time_lims,2),...
      SCORE(1:n_time_lims,3),'b','LineWidth',1.5)
SCORE(1:n_time_lims,:)=[];
end;

grid on;  

