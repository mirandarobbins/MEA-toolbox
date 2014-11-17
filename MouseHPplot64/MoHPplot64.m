
function MoHPPlot64(inFilename, DESIRED_TRACES, MEM_LIMIT, Yvalue, Xmin, Xmax)

%MoHPPlot64: Load data from a Conductor file and Plot as 002 probe order (Mouse Hippocampus). 
%	MoPlot64(inFilename, DESIRED_TRACES [1], MEM_LIMIT [100 MB], Yvalue, Xmin, Xmax)
%       inFilename: 
%       DESIRED_TRACES: target number of sweep (default; 1)
%       MEM_LIMIT: Meory limit of this import rutin (default; 100MB)
%       Yvalue: scale limit of Y axis (mV)
%       Xmin, Xmax: range of X axis (msec); If there is no info of X axis, X axis has AUTO mode


%--- Parameters --

if( ~exist( 'MEM_LIMIT', 'var' ) |  isempty( MEM_LIMIT ) )
   MEM_LIMIT = 100;  % in MB. limit of memory size available to load data in. 
end

if( ~exist( 'DESIRED_TRACES', 'var' ) |  isempty( DESIRED_TRACES ) )
   DESIRED_TRACES = 1;  % trace to load 
end


%--- Load data ---
data = medload(inFilename, DESIRED_TRACES, MEM_LIMIT);

%--- Creat figure window ----
scrsz = get(0,'ScreenSize');
h_ca1 = figure('Position',[1 scrsz(4)*2/3-75 scrsz(3)/2 scrsz(4)/3]);
h_ca3 = figure('Position',[scrsz(3)/2 scrsz(4)/3 scrsz(3)/2 scrsz(4)/3]);
h_dg = figure('Position',[1 25 scrsz(3)/2 scrsz(4)/3*1.25]);
        
%--- Creat plotting axis (CA1) ----
    figure(h_ca1);
    for axis_box_index =1:24
        if (axis_box_index <= 5)
            h(axis_box_index)=axes('position', [.05 + axis_box_index * .14, .8, .13, .2])    
        elseif (axis_box_index <= 11)
            h(axis_box_index)=axes('position', [.12 + (axis_box_index - 6) * .14, .55, .13, .2])
        elseif (axis_box_index <= 18)
            h(axis_box_index)=axes('position', [.05 + (axis_box_index - 12) * .14, .3, .13, .2])
        else
            h(axis_box_index)=axes('position', [.12 + (axis_box_index - 19) * .14, .05, .13, .2])
        end
    end    
    
    figure(h_ca3);
    for axis_box_index = 25:43
        if (axis_box_index <= 28)
            h(axis_box_index)=axes('position', [.29 + (axis_box_index - 25) * .14, .8, .13, .2])
        elseif (axis_box_index <= 33)
            h(axis_box_index)=axes('position', [.22 + (axis_box_index - 29) * .14, .55, .13, .2])
        elseif (axis_box_index <= 38)
            h(axis_box_index)=axes('position', [.15 + (axis_box_index - 34) * .14, .3, .13, .2])
        else
            h(axis_box_index)=axes('position', [.08 + (axis_box_index - 39) * .14, .05, .13, .2])
        end
    end

    figure(h_dg);
    for axis_box_index = 44:64       
        if (axis_box_index <= 46)
            h(axis_box_index)=axes('position', [.31 + (axis_box_index - 44) * .14, .84, .13, .16])
        elseif (axis_box_index <= 50)
            h(axis_box_index)=axes('position', [.24 + (axis_box_index - 47) * .14, .64, .13, .16])
        elseif (axis_box_index <= 55)
            h(axis_box_index)=axes('position', [.17 + (axis_box_index - 51) * .14, .44, .13, .16])
        elseif (axis_box_index <= 60)
            h(axis_box_index)=axes('position', [.1 + (axis_box_index - 56) * .14, .24, .13, .16])    
        else
            h(axis_box_index)=axes('position', [.17 + (axis_box_index - 61) * .14, .04, .13, .16])
        end
    end

%--- re-plott ---
    axes(h(1)) % ch 11
        plot(data(:,1), data(:,12), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(2)) % ch 20
        plot(data(:,1), data(:,21), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(3)) % ch 4
        plot(data(:,1), data(:,5), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(4)) %ch 5
        plot(data(:,1), data(:,6), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(5)); % ch 21
        plot(data(:,1), data(:,22), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(6)); % ch 10
        plot(data(:,1), data(:,11), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(7)); % ch 2
        plot(data(:,1), data(:,3), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(8)); % ch 12
        plot(data(:,1), data(:,13), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(9)); % ch 13
        plot(data(:,1), data(:,14), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(10)); % ch 6
        plot(data(:,1), data(:,7), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(11)); % ch 14
        plot(data(:,1), data(:,15), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(12)); % ch 9
        plot(data(:,1), data(:,10), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end    
    axes(h(13)); % ch 1
        plot(data(:,1), data(:,2), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end    
    axes(h(14)); % ch 3
        plot(data(:,1), data(:,4), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(15)); % ch 29
        plot(data(:,1), data(:,30), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(16)); % ch 8
        plot(data(:,1), data(:,9), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(17)); % ch 22
        plot(data(:,1), data(:,23), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(18)); % ch 7
        plot(data(:,1), data(:,8), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(19)); % ch 18
        plot(data(:,1), data(:,19), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(20)); % ch 17
        plot(data(:,1), data(:,18), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(21)); % ch 19
        plot(data(:,1), data(:,20), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(22)); % ch 27
        plot(data(:,1), data(:,28), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(23)); % ch 26
        plot(data(:,1), data(:,27), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(24)); % ch 15
        plot(data(:,1), data(:,16), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(25)); % ch 16
        plot(data(:,1), data(:,17), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(26)); % ch 23
        plot(data(:,1), data(:,24), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(27)); % ch 24
        plot(data(:,1), data(:,25), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(28)); % ch 30
        plot(data(:,1), data(:,31), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(29)); % ch 62
        plot(data(:,1), data(:,63), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(30)); % ch 46
        plot(data(:,1), data(:,47), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(31)); % ch 37
        plot(data(:,1), data(:,38), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(32)); % ch 31
        plot(data(:,1), data(:,32), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(33)); % ch 32
        plot(data(:,1), data(:,33), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(34)); % ch 63
        plot(data(:,1), data(:,64), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(35)); % ch 64
        plot(data(:,1), data(:,65), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(36)); % ch 48
        plot(data(:,1), data(:,49), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(37)); % ch 39
        plot(data(:,1), data(:,40), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(38)); % ch 40
        plot(data(:,1), data(:,41), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(39)); % ch 54
        plot(data(:,1), data(:,55), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(40)); % ch 55
        plot(data(:,1), data(:,56), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(41)); % ch 56
        plot(data(:,1), data(:,57), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(42)); % ch 47
        plot(data(:,1), data(:,48), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(43)); % ch 38
        plot(data(:,1), data(:,39), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(44)); % ch 34
        plot(data(:,1), data(:,35), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(45)); % ch 33
        plot(data(:,1), data(:,34), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(46)); % ch 25
        plot(data(:,1), data(:,26), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(47)); % ch 41
        plot(data(:,1), data(:,42), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(48)); % ch 35
        plot(data(:,1), data(:,36), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(49)); % ch 28
        plot(data(:,1), data(:,29), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(50)); % ch 45
        plot(data(:,1), data(:,46), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(51)); % ch 49
        plot(data(:,1), data(:,50), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(52)); % ch 42
        plot(data(:,1), data(:,43), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(53)); % ch 50
        plot(data(:,1), data(:,51), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(54)); % ch 61
        plot(data(:,1), data(:,62), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(55)); % ch 53
        plot(data(:,1), data(:,54), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(56)); % ch 57
        plot(data(:,1), data(:,58), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(57)); % ch 43
        plot(data(:,1), data(:,44), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(58)); % ch 59
        plot(data(:,1), data(:,60), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(59)); % ch 52
        plot(data(:,1), data(:,53), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(60)); % ch 36
        plot(data(:,1), data(:,37), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(61)); % ch 58
        plot(data(:,1), data(:,59), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(62)); % ch 51
        plot(data(:,1), data(:,52), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(63)); % ch 44
        plot(data(:,1), data(:,45), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
    axes(h(64)); % ch 60
        plot(data(:,1), data(:,61), 'black');
        ylim([-Yvalue Yvalue]);
         if nargin >= 5
            xlim([Xmin Xmax]);
        else
            xlim('auto')
        end
        