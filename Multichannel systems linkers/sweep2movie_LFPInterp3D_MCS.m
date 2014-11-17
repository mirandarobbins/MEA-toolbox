function M = sweep2movie_LFPInterp3D(DuV, startend, varargin);
% sweep2movie - generate matlab movie from a MEA data set
% M = sweep2movie(DuV, startend, varargin), were <DuV> is 
% a 64xN or 68xN matrix and <startend> is a 2-element vector
% delimiting a range of columns in <DuV>, generates a
% matlab-type movie matrix that can be played with movie(M).
% 
% startend must be in column indices or ms. If given in ms
% it must be combined with the parameter/value pair 
% ..., 'startendScale', 'ms', ...
%
% Further parameter/value options, default in parentheses
%
% 'MicrosecondsPerTick'    sampling resolution in µs
% 'clim'                   color scale, (0.9*[min(min(DuV)) max(max(DuV))])
% 'ipf'                    interpolation depth  (2)
% 'skip'                   no. of samples to skip between frames (1)
% 'framesize'              2-element for image size given in pixels ([200 200])
% 'memlim'                 limit before a memory warning is issued (256 Mb)
% 'showTime'               'yes' to show the frame time in ms ('yes')
% 'displayTimeStart'       start time to display in the movie (startframe *(MicrosecondsPerTick*1e-3)) 
%
% (c) U. Egert
%
% See also MPEGBUILDER MPGWRITE
% cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blueblackyellow.mat');blue_white_red
% cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat');
% cmap=flipud(cmap.cmap); %red is sink
% cmap=cmap.cmap; %red is source
cmap=jet
MicrosecondsPerTick = 50;
clim = 'auto';
ipf = 4;
skip = 1;
framesize = [200 200];
memlim = 256e6;
startendScale = 'frames';
showTime = 'yes';
angle=[30 30];
%--------------------------------------------------------------------
% evaluate input parameters
%--------------------------------------------------------------------
pvpmod(varargin);

if ~isvar('displayTimeStart')
   startTime = startframe *(MicrosecondsPerTick*1e-3);
end;

switch lower(startendScale)
case 'ms'
   startframe = floor(startend(1)/(MicrosecondsPerTick*1e-3));
   endframe = floor(startend(2)/(MicrosecondsPerTick*1e-3));
case 'frames'
   startframe = startend(1);
   endframe = startend(1);
end;

frames = 1:skip:diff([startframe, endframe]);
nframes = length(frames);

%--------------------------------------------------------------------
% prepare the figure
%--------------------------------------------------------------------
fh = setgcf('make movie', 'replace');
set(fh, 'renderer', 'zbuffer', 'pos',  [3 36 1184 943], 'menubar', 'none');
set(fh, 'units', 'normal');
p=get(fh, 'pos'); 
set(fh, 'pos', [0.01 .99-p(4) p(3:4)]);
%set(fh, 'renderer', 'OpenGL', 'units', 'normal', 'pos', [0 1-framesize(2) framesize ], 'menubar', 'none');

if strcmp(clim , 'auto')
%    clim  = [min(DuV(:)) max(DuV(:))]*.9;
end;
set(gca, 'visible', 'off', 'clim', 0.3*clim, 'nextplot', 'replacechildren', ...
    'plotboxaspectratio', [1 1 1], 'ydir', 'reverse', 'units', 'normal', 'pos', [0 0 1 1]);
axh = gca;
% view(2)

%--------------------------------------------------------------------
% generate the first frame and tests
%--------------------------------------------------------------------

% w = waitbar(0, 'processing frame ...');
% set(w, 'name', 'progress', 'units', 'normal');%, 'pos', [0.64 0.08 0.35 0.1]);
imgdata = interp2(reshape(DuV(:,startframe+frames(1)), 8,8), ipf, 'cubic');
ih = surf(imgdata'); 
view(angle)
   set(ih, 'cdata', imgdata',...
           'facecolor', 'interp',...
           'edgecolor', 'none',...
           'FaceLighting','phong ',...
           'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit');
    lightangle(30,30)
    colormap(cmap)
    
    axis off; box off; grid off ;
%     set(gca,'XTick',1:17.5:500,'XTickLabel','0|150|300|450|600|750|900|1050|1200')
%     xlabel('Electrode pitch(\mum)')
   % set(gca,'YTick',1:15:120,'YTickLabel','0|150|300|450|600|750|900|1050|1200')
   % ylabel('Electrode pitch(\mum)')
%     set(gca,'YTick',1:17.5:500,'YTickLabel','Pia | Layer II|III | Layer IV |  Layer IV | Layer Va |  Layer Vb')
%     ylabel('Cortical Lamina')
%     zlabel('LFP amplitude (mV)')
 set(gca,'zlim',clim);
 zoom(0.6)

% drawnow

S = length(moviein(1))*nframes;
if S*8 > memlim
   btn=questdlg(['generating the movie with the current settings will require ' num2str(S*8/1e6) ' Mbyte, do you want to continue?'], 'sweep2movie', 'no');
   switch lower(btn)
   case 'yes'
      clear S
   otherwise
      return
   end;
end

%--------------------------------------------------------------------
% generate the movie
%--------------------------------------------------------------------
M = moviein(nframes);
M(:,1) =   getframe(axh);
setfield(M(1),'colormap',cmap);

for i = 2:nframes-1;
   imgdata = interp2(reshape(DuV(:,startframe+frames(i)), 8,8), ipf, 'cubic');
   ih=surf(imgdata');
   
   view(angle)
   set(ih, 'cdata', imgdata',...
           'facecolor', 'interp',...
           'edgecolor', 'none',...
           'FaceLighting','phong ',...
           'AmbientStrength',.3,'DiffuseStrength',.8,...
    'SpecularStrength',.9,'SpecularExponent',25,...
    'BackFaceLighting','unlit');
    lightangle(30,30)
    colormap(cmap)
    
    axis off; box off; grid off 
%     set(gca,'XTick',1:17.5:500,'XTickLabel','0|150|300|450|600|750|900|1050|1200')
%     xlabel('Electrode pitch(\mum)','Rotation',-14);
   % set(gca,'YTick',1:15:120,'YTickLabel','0|150|300|450|600|750|900|1050|1200')
   % ylabel('Electrode pitch(\mum)')
%     set(gca,'YTick',1:25:500,'YTickLabel','Pia | Layer II | Layer III |  Layer IV | Layer Va |  Layer Vb')
%     ylabel('Cortical Lamina','Rotation',43);
%     zlabel('LFP amplitude (mV)')

    if strcmp(showTime, 'yes')
       th =text(1,1, [num2str((displayTimeStart+frames(i))*(MicrosecondsPerTick*1e-3)-0.1) ' ms']);
       set(th, 'fontsize', 14, 'color', 'k', 'units', 'normal', 'pos', [.89 .9 clim(2)]);
    else
       th =text(1,1, '');
       set(th, 'fontsize', 14, 'color', 'k', 'pos', [1.5 25 clim(1)-1]);
    end

%     waitbar(i/nframes);

    M(:,i) =   getframe(axh);
end;

for i = nframes-1
setfield(M(i),'colormap',cmap);
end
% delete(w);
