%txt2mat.m
%
%	convert Medac data file from ASCII text format to matlab binary format
%	source file will be deleted after conversion
%	a faulty channel can be replaced by interpolation of neighbors
%	'label' appears on plot titles that use the data file
%	Fs: sampling frequency in KHz
%	txt2mat( filename, outFilename, label, Fs, nTraces, badChannel [default = 0 = none] ) 

function txt2mat( filename, outFilename, label, Fs, nTraces, badChannel )


%--- default parameters ---

if( nargin == 4 ), badChannel = 0;, end


%--- Load data file ---

disp( '[txt2mat] Loading data...' )
eval( ['data = load(''', filename, ''');'] )


%--- Fix bad channel ---

switch badChannel   
case 0
   disp( '[txt2mat] No invalid data channel' )
case 1
   disp( ['[txt2mat] Fixing invalid data of channel: ', num2str( badChannel )] )   
   data(:,2) = ( data(:,3) + data(:,10)) / 2;
case 64
   disp( ['[txt2mat] Fixing invalid data of channel: ', num2str( badChannel )] )   
   data(:,65) = ( data(:,64) + data(:,57)) / 2;
otherwise
   disp( ['[txt2mat] Fixing invalid data of channel: ', num2str( badChannel )] )
   data(:,badChannel+1) =  ( data(:,badChannel+1-1) + data(:,badChannel+1+1) ) / 2;   
end


%--- If Fs is not 10 KHz, resample data to make it so ---

if( Fs ~= 10 )
   disp( '[txt2mat] Resampling to 10KHz...' )
   desiredPeriod = 0.1;  % in [ms]
   
   tempData = [];
   sizeTrace = size(data,1) / nTraces;
   endTime = data(sizeTrace,1);
   timeInt = 0:desiredPeriod:endTime;
   
   for trace = 1:nTraces
      time = data(1 + (trace - 1) * sizeTrace:sizeTrace + (trace - 1) * sizeTrace, 1); 
      traceData = data(1 + (trace - 1) * sizeTrace:sizeTrace + (trace - 1) * sizeTrace, 2:65); 
      tempData = [tempData; timeInt' interp1( time , traceData, timeInt )];
   end
   
   data = tempData;
end


%--- Save data in matlab format ---

disp( '[txt2mat] Saving data...' )
eval( ['save ', outFilename, ' data label'] )
delete( filename )
clear

