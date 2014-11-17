function CSD=AD_iCSD(data,params)
%% Setup

old_dir=pwd;
cd /Users/aleksanderdomanski/Desktop;
start_dir=pwd;

if isdir([start_dir,'/data'])
    rmdir ('data','s')
end
mkdir('data');
clear input_data

nx = 8;
ny = 8;
dx = 0.15;
dy = 0.15;
h  = 0.01;
zprofile='step';

VX = 1:8;%1:1:nx;
VY = 1:8;%1:1:ny;


CSD.params.nx = nx;
CSD.params.ny = ny;
CSD.params.dx = dx;
CSD.params.dy = dy;
CSD.params.h  = h;
CSD.params.zprofile  = zprofile;
fname=['temp: ',data.this_file];
name='8x8 map';
%% Prepare F- matrix in temp ('data') directory
initspline2d(name,nx+2,ny+2,dx,dy,h,zprofile) % use 2d spline interpolation
% initlin2d(name, nx, ny, dx, dy, h); % use linear 2d interpolation
%% iteratively calculate CSD

for trial_id=data.sweep_sort.successful_sweeps
    ['Analysing trial no.: ',num2str(trial_id),' of ',num2str(numel(data.sweep_sort.successful_sweeps))]
    cd ([start_dir,'/data'])
        input_data=data.filtered_lfp(:,:,trial_id);
%         input_data=data.mean_channels; % use grand mean
        input_data=reshape(input_data,20000,8,8);
        save(fname, 'input_data', 'name', 'dx', 'dy', 'h', ...
        'nx', 'ny','zprofile'); 
    
    cd(start_dir)
%     CSD_grid=load(['data/',name], 'F');  CSD_grid=CSD_grid.F; % for linear
    CSD_grid=load(['data/',name], 'Fn');  CSD_grid=CSD_grid.Fn; % For splines
    csd1=icsd2d(input_data,CSD_grid,'D');
%     csd2 = interp2d(csd1, VX, VY, 'splinem');

    CSD.csd_array(:,:,trial_id) = reshape(squeeze(csd1),20000,64,1); % linear
%     CSD.csd_array(:,:,trial_id) = reshape(squeeze(csd2),20000,size(csd2,2)^2,1); % spline
['Finished trial no.: ',num2str(trial_id),' of ',num2str(numel(data.sweep_sort.successful_sweeps))]
end
%% clean up
clear ans nx ny dx dy h zprofile fname name trial_id input_data  CSD_grid start_dir csd1 csd2 VX VY
rmdir ('data','s')
cd (old_dir)
%% mean CSD for monosynaptic and window
CSD.mean_monosynaptic=zeros(8,8,size(CSD.csd_array,3));
for sweep_id=1:size(CSD.csd_array,3)
    for channel_id=1:size(CSD.csd_array,2)
        temp(channel_id,sweep_id)=mean(CSD.csd_array(params.monosynaptic_win_samples(1):params.monosynaptic_win_samples(2),...
                                                        channel_id,sweep_id));
    end 
    CSD.mean_monosynaptic(:,:,sweep_id)=reshape(temp(:,sweep_id),8,8);
end
CSD.grand_mean_monosynaptic=mean(CSD.mean_monosynaptic,3);
CSD.grand_mean_monosynaptic=reshape(CSD.grand_mean_monosynaptic,8,8);
%% max STD in window
CSD.sweep_sort.detection_threshold=20;
    temp_win      =  CSD.csd_array(params.search_win_samples(1):params.search_win_samples(2),:,:);
    temp_baseline =  CSD.csd_array(params.baseline_win_samples(1):params.baseline_win_samples(2),:,:);
    
for channel_id=1:params.no_channels
        for sweep_id=data.sweep_sort.successful_sweeps
            CSD.sweep_sort.window_std(channel_id,sweep_id) =nanstd(temp_win(:,channel_id,sweep_id))./...
                                                             nanstd(temp_baseline(:,channel_id,sweep_id));
        end
end
    CSD.sweep_sort.window_std_max=nanmax(CSD.sweep_sort.window_std);

clear channel_id sweep_id
assignin('base', 'CSD', CSD);