% multiplot64: plot MED graph 8x8 windows
%   data=loaded data
%   Yvalue=scale limit of Y axis
%   Xmin, Xmax=range of X axis
% mulliplot(data, Yvalue, Xmin, Xmax)

function multiplot64(data, Yvalue, Xmin, Xmax)
for i=1:64
        subplot(8,8,i);
        plot(data(:,1), data(:,i+1));
        ylim([-Yvalue Yvalue]);
        if nargin >= 3
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
end
