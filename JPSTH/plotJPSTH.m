function plotJPSTH(temp,chan_1,chan_2,tb,bw)
if nargin==1
chan_1=1;
chan_2=2;
tb=1:499;
bw=1;
end
fhandle=figure();
mataxes=gca;
ph=fhandle;

cmap=load('/Users/aleksanderdomanski/Documents/Aleks Domanski/MATLAB/colourmaps/blue_white_red.mat');
cmap=cmap.cmap;


% Problems with OpenGL - user painters
set(fhandle, 'Renderer', 'painters');
% Orient
orient(fhandle, 'landscape');

% plot JPSTH matrix
matrix=temp.normalized_jpsth(tb,tb);
hold on
pl=contourf(mataxes,matrix, 'LineColor', 'none');
pl_2=plot(mataxes,tb,':k');
 set(mataxes, 'XLim', [tb(1) tb(end)],...
        'YLim', [tb(1) tb(end)],...
        'ZLim', [min(min(matrix)) max(max(matrix))]);
%                 set(pl, 'LineColor', 'none');
        set(mataxes, 'YDir', 'normal',...
    'Units', 'normalized',...
    'XTickLabelMode', 'manual',...
    'XTickLabel',[],...
    'YTickLabelMode', 'manual',...
    'YTickLabel',[]);
colormap(cmap)
title(mataxes, 'PSTH');
set(mataxes, 'Position', [0.2 0.3 0.25 0.25]);
set(mataxes, 'Units', 'pixels');
p=get(mataxes, 'Position');


% Set data limits for colormap: this ignores NaNs
mn=min(min(matrix));
mx=max(max(matrix));
if mn>=0
    csc=[0 mx];
elseif mn<0 && mx<0
    csc=[mn 0];
else
    csc=max(abs([mn,mx]));
    csc=[-csc csc];
end

% Set colormap
set(mataxes, 'CLimMode', 'manual',...
    'Clim', csc);
t=max([p(3) p(4)]);
set(mataxes, 'Position', [p(1) p(2) t t]);
p=get(mataxes, 'Position');
set(mataxes, 'Units', 'normalized');


%PETH1
h=axes('Parent', ph, 'Position', [0.2 0.1 0.25 0.15]);
set(h, 'Units', 'pixels');
p1=get(h, 'Position');
p1(3)=t;
set(h, 'Position', p1);
x=temp.psth_1(tb);
h1=bar(tb,x,'k');
set(h1, 'Hittest', 'off');
set(h, 'XLim', [min(tb) max(tb)]);
set(h, 'Color', 'none', 'GridLineStyle', 'none', 'Box', 'off',...
    'YAxisLocation', 'left');
set(h, 'Units', 'normalized');
x=get(h, 'XLim');
y=get(h,'YLim');
title(h, 'PETH 1', 'Position', [x(1) y(1)-0.25],...
    'HorizontalAlignment', 'left');


% PETH2
h=axes('Parent', ph, 'Position', [0.01 0.3 0.15 0.25]);    
set(h, 'Units', 'pixels');
p2=get(h, 'Position');
p2(4)=t;
set(h, 'Position', p2);
x=temp.psth_2(tb);
h2=bar(tb, x,'k');
set(h, 'XLim', [min(tb) max(tb)]);
view(-90, 90);
p2(1)=p2(1)+p2(3)-p1(4)+10;
p2(3)=p1(4);
set(h, 'Position', p2);
set(h, 'Color', 'none', 'GridLineStyle', 'none', 'Box', 'off');
set(h, 'Units', 'normalized',...
    'ButtonDownFcn', {@LocalCallback tb x});
x=get(h, 'XLim');
y=get(h,'YLim');
title(h, 'PETH 2', 'Position', [x(1)-220 y(2)],...
    'HorizontalAlignment', 'left');

% Coincidence histogram
h=axes('Parent', ph, 'Position', [0.5075   0.3967    0.4882    0.6]);
x=temp.pstch(tb); x(lt(x,0))=0;
h3=bar3(tb, x, 'histc');
set(h3,'LineWidth',0.5)
set(h, 'Color', 'none', 'GridLineStyle', 'none', 'Box', 'off');
view(-145,43);
set(h, 'XLim', [min(tb) max(tb)]);
set(h, 'YLim', [min(tb) max(tb)]);
set(h, 'ZLim', [0 1.5*max(x)]);
z=get(h,'ZLim');
title(h, 'Coincidence', 'Position', [0 0 z(1)-3],...
    'HorizontalAlignment', 'center',...
    'VerticalAlignment', 'bottom');


% XCorr
h=axes('Parent', ph, 'Position', [0.2321    0.6    0.5258    0.3]);
tb= -50:50;x=temp.xcorr_hist;
h4=bar3(h,tb, x, 'histc');
set(h4,'LineWidth',0.5)
% hold on
% h4=bar(tb, x, 'histc');
% h5=stairs(tb, temp.sig_low,'r')
% h6=stairs(tb, temp.sig_high,'r')
view(-145,-47);
set(h, 'Color', 'none', 'GridLineStyle', 'none', 'Box', 'off');
set(h, 'XLim', [min(tb) max(tb)]);
set(h, 'YLim', [min(tb) max(tb)]);
% set(h, 'Units', 'normalized', 'ButtonDownFcn', {@LocalCallback tb*1e3 x});
z=get(h,'ZLim');
% set(h,[z(1)-2 z(2)]);
title(h, 'Correlation', 'Position', [0 0 z(2)],...
    'HorizontalAlignment', 'left');

% Colorbar
drawnow();
h=colorbar('peer', mataxes);
set(h, 'Position', [0.9 0.05 0.025 .3]);

% Set small font for clarity
h=findall(ph, 'Type', 'text');
h=[h; findall(ph, 'Type', 'axes')];
set(h, 'FontUnits', 'points', 'FontSize', 10);
set(h, 'FontUnits', 'normalized');

% Add title or name
set(ph, 'Name',['JPSTH for MUA on channel ',num2str(chan_1),' vs. channel ', num2str(chan_2)]);
