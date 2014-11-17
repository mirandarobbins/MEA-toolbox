function AD_MED64_LFP_movie_Miranda(data)
clear input temp M
% params.selected_rep=1;
for trial_id=data.sweep_sort.successful_sweeps
%     close all
% figure('Name',strcat('LFP movie for trial ',num2str(trial_id)),'NumberTitle','off')
input=squeeze(data.filtered_lfp(:,:,trial_id)');
% input=squeeze(data.filtered_lfp(:,:,params.selected_rep)');
% input=data.mean_channels';
% input=CalCSD_AD(temp);

movie.cfg_idxbad_pads=[];
movie.cfg_idxbad_nonip=[];
movie.cfg_ipDepth=1;
movie.cfg_ipMethod='cubic';
movie.cfg_title='measurfplot';
movie.cfg_ystring='voltage [ \muV]' ;
movie.cfg_zlim=[-100 +100];%[min(input(:)) max(input(:))]
movie.cfg_zlim_blank=0;
movie.cfg_inv_cmap=0;
movie.cfg_cldata=0;

% [data_nanidx, surfobjh] = measurfplot(input, movie);

% CSD

% clear M M_smooth; close all
% writerObj = VideoWriter([data.this_file,'_trial_',num2str(trial_id),'.avi'])
% writerObj.FrameRate = 10;
% open(writerObj);
M=sweep2movie_LFPInterp3D_Miranda(input,[1 990],...
              'MicrosecondsPerTick',500,...
              'clim',[-0.0005 0.0005],...  %0.9*[min(min(temp)) max(max(temp))]
              'ipf',1,...
              'skip',1,...
              'framesize',[500 500],...
              'memlim',512e6,...
              'showTime','no',...
              'startendScale','ms',...
              'displayTimeStart',1,...
              'angle',[30 30]);
% %                     movieview(M);
% frame = getframe;
%    writeVideo(writerObj,frame);

%            close all;

% %    Write as .avi
% last_dir=pwd;
% if  ~isdir('movies')
%     mkdir('movies')
% end
% target_dir=strcat(last_dir,'\movies');
% cd movies
% movie2avi(M, ['LFP_',data.this_file,'_trial_',num2str(trial_id),'.avi'], 'compression', 'None')
% cd(last_dir)
close(gcf)

end
clear M movie trial_id input last_dir target_dir
disp( 'finished coding AVI files!')
%           mpgwrite(M, cmap, [pwd,'\',data.this_file])
% for I=1:size(M,2),
% %                plot_command
%                a(I) = getindexedframe;
%             end;
%             aviwrite(strcat(pwd,'\data.this_file.avi'),M,'menu')
%%
% cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat');
% cmap=cmap.cmap;
% 
% % target_dir	 = [pwd,'/' data.this_file];% 'c:\temp\jpeg_sequence';
% 
% for frameidx = 1:size(M,2)
%    [framedata, cmap2] = frame2im(M(:,frameidx));
%    frameName = ['0000' num2str(frameidx)];
%    imwrite(double(framedata),cmap, fullfile( target_dir, [data.this_file '_' frameName(end-4:end) '.jpg']), 'jpeg');
% end;