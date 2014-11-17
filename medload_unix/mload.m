
function [data] = mload( inFilename, DESIRED_TRACES, MEM_LIMIT )

%MLOAD	Loads data from a medac data file (Spontaneous Data). 
%	[data] = medload( inFilename, DESIRED_TRACES [1], MEM_LIMIT [20 MB] ) 
%	loads the desired trace (one is default) from the data file defined 
%	by 'inFilename', respecting the desired memory usage limit defined 
%	by MEM_LIMIT (default 10 MB). If DESIRED_TRACES is 'all' or 'avg' or
%	'average' all the traces are loaded and averaged.
%	If the data file is larger than the defined memory limit it will be
%	loaded until the limit is reached. The resulting loaded data is 
%	truncated in time.
%	This function requieres both the '.med' and '.dat' files generated 
%	by a medac session. The input filename can be the name of any of them
%	of just their common name without the extension.
%	The first column of 'data' contains the time stamp, and the remainder 
%	columns contain the recording channels in order, with column 2 of 
%	'data' corresponding to channel one. The first row of 'data' corresponds
%	to time = 0. 
%
%   Bug fixed by T.H. on 2004/12/21


%--- Parameters --

if( ~exist( 'MEM_LIMIT', 'var' ) |  isempty( MEM_LIMIT ) )
   MEM_LIMIT = 20;  % in MB. limit of memory size available to load data in. 
end

if( ~exist( 'DESIRED_TRACES', 'var' ) |  isempty( DESIRED_TRACES ) )
   DESIRED_TRACES = 1;  % trace to load 
end

VERSION_NUMBER = 2.1;
VERSION_LABEL  = 'Version';
ONE_SHOT_SIZE_LABEL = 'One Shot Size/CH';
DAQ_CHANNELS_LABEL = 'DAQ Channels/BD';
GAIN_LABEL = 'Gain';
AMPLITUDE_LABEL = 'Amplitude(mV/V)';
%ONE_RECORD_SIZE_LABEL = 'One Record Size/CH';
SAMPLE_RATE_LABEL = 'Sample Rate(kHz)/CH';
FRAGMENTS_LABEL = 'Fragment(s)';
FRAG_SIZE_LABEL = 'One Fragment Size/CH';


%--- Input parameter must be a char string ---

if( ~ischar( inFilename ) ), error( 'Input parameter must be a char string' ), end


%--- Test for valid filename extension (if any) and strip it ---

[inFilename, extension] = strtok( inFilename, '.' );

if( ~isempty( extension ))
   if( ~strcmp( extension, '.med' ) & ~strcmp( extension, '.dat' ) )
      error( 'Input filename extension is not ''.dat'' or ''.med''. Use a valid extension or omit it.' )
   end
end


%--- Create the filenames for the header and data files associated to the input filename ---

inFilenameHead = [inFilename,'.med'];
inFilenameData = [inFilename,'.dat'];


%--- Test if the input file exists ---

if( ~existFile( inFilenameHead ) ), error( ['Unable to find file named: ', inFilenameHead] ), end
if( ~existFile( inFilenameData ) ), error( ['Unable to find file named: ', inFilenameData] ), end


%===== Read medac header file information =====

fh = fopen( inFilenameHead );  % open header file


%--- Check that format version of data file is correct ---

%versionStr = readLabeledString( fh, VERSION_LABEL );
%versionNum = str2num( versionStr(1:3) );

%if( versionNum ~= VERSION_NUMBER )
%   fclose( fh );   
%   error( ['Data file version number is not: ', num2str( VERSION_NUMBER ), '. Feature not supported.'] )
%end


%--- Extract number of channels per board, and number of boards ---

chPerCard = readLabeledValue( fh, DAQ_CHANNELS_LABEL );
nDaqCards = length( chPerCard );

if( max(chPerCard) ~= min(chPerCard) )
   fclose( fh );   
   error( 'DAQ cards have different number of channels. Feature not supported.' )
end

chPerCard = chPerCard(1);
nChannels = nDaqCards * chPerCard;


%--- Extract various data file parameters ---

oneShotSize = readLabeledValue( fh, ONE_SHOT_SIZE_LABEL );  % block size
gain = readLabeledValue( fh, GAIN_LABEL );  % gain
amplitude = readLabeledValue( fh, AMPLITUDE_LABEL );  % amplitude
%oneRecordSize = readLabeledValue( fh, ONE_RECORD_SIZE_LABEL );  % record size
samplingFreq = readLabeledValue( fh, SAMPLE_RATE_LABEL );  % sampling rate, in kHz
nFragments = readLabeledValue( fh, FRAGMENTS_LABEL );  % number of traces (fragments)
fragmentSize = readLabeledValue( fh, FRAG_SIZE_LABEL );  % size per channel of a fragment 


%--- Close header file ---

fclose( fh );


%===== Read medac data file contents =====

%fd = fopen( inFilenameData );  
fd = fopen( inFilenameData, 'r', 'l' ); % <== open data file in 'l' format. (T.H. 2004/12/21)


%--- Compute and constrain memory size of loaded data ---

%nBlocks = floor( oneRecordSize / oneShotSize );
nBlocks = floor( fragmentSize / oneShotSize);

blockDataMemSize = 8 * nChannels * oneShotSize;
blockTimeStampMemSize = 8 * oneShotSize;
blockMemorySize = blockTimeStampMemSize + blockDataMemSize;

dataMemorySize = nBlocks * blockMemorySize;

if( dataMemorySize / 1e6 > MEM_LIMIT )
   warning( ['Data file cannot be loaded in the defined memory limit of: ', num2str(MEM_LIMIT), ' MB. The loaded data is truncated in time.'] )
   nBlocks = floor( MEM_LIMIT * 1e6 / blockMemorySize );
end


%=== Multiple trace read loop (one trace is just boring special case of this) ===

if( strcmpi( DESIRED_TRACES, 'all' ) | strcmpi( DESIRED_TRACES, 'avg' ) | strcmpi( DESIRED_TRACES, 'average' ) )  % in this case all traces are averaged
   
   DESIRED_TRACES = [1:nFragments];
   
else
   
   if( min( DESIRED_TRACES ) < 1 | max( DESIRED_TRACES ) > nFragments )
      fclose( fd );
      error( ['Desired trace: ', num2str( DESIRED_TRACES ), ...
            ', does not exist in the input data file. Choose from the range: 1 to ', num2str( nFragments),'.'] ) 
   end
   
   firstTrace = DESIRED_TRACES;
   lastTrace = DESIRED_TRACES;
   
end

oneTraceSize = fragmentSize * 2 * nChannels;
data = zeros( nBlocks * oneShotSize, nChannels + 1 );  % first column for the time stamp
%scaleFactor = amplitude * gain * 10 / 2048 / 50;
scaleFactor = amplitude / 1024 * 2.5 / gain;     %<=== (T.H. 2004/12/21)

for traceNum = DESIRED_TRACES
   
   %--- Seek to the begining of the desired fragment (trace) ---
   
   fseek( fd, oneTraceSize * (traceNum - 1), 'bof' );
   
   
   %--- Data load loop ---
   
   for block = 1:nBlocks
      
      switch nDaqCards 
         
      case 1
         dataBlock = [fread( fd, [chPerCard, oneShotSize], 'int16' )'];
         
      case 2
         dataBlock = [fread( fd, [chPerCard, oneShotSize], 'int16' )', ...
               fread( fd, [chPerCard, oneShotSize], 'int16' )'];   
         
      case 4   
         dataBlock = [fread( fd, [chPerCard, oneShotSize], 'int16' )', ...
               fread( fd, [chPerCard, oneShotSize], 'int16' )', ...
               fread( fd, [chPerCard, oneShotSize], 'int16' )', ...
               fread( fd, [chPerCard, oneShotSize], 'int16' )'];
         
      otherwise
         fclose( fd );
         error( [ num2str( nDaqCards ), ' number of DAQ cards not supported.'] )
         
      end
      
      data(1 + (block - 1) * oneShotSize:oneShotSize + (block - 1) * oneShotSize, 2:nChannels + 1) = ...
         data(1 + (block - 1) * oneShotSize:oneShotSize + (block - 1) * oneShotSize, 2:nChannels + 1) + dataBlock;
      
   end
   
end

nTracesRead = length( DESIRED_TRACES );
data = data / nTracesRead;
data = data * scaleFactor;  %in units of mV

data(:, 1) = [0:length( data ) - 1]';	  % generate the time stamp in ms.
data(:, 1) = data(:, 1) / samplingFreq;  % generate the time stamp in ms.




%=================================================================================================

function [value] = readLabeledValue( f, label )


frewind( f )  % start search at the begining of the file
oneLine = fgetl( f );

while( ~strncmp( label, oneLine, length( label ) ) & ~feof( f ) )
   oneLine = fgetl( f ); 
end

if( feof( f ) )
   fclose( f );   
   error( ['Fail to find label: ', label, ', in provided file'] )
end

[labelTokenized, valueTokenized] = strtok( oneLine, '=' );
valueString = valueTokenized(2:length( valueTokenized ));
value = str2num( valueString );


%=================================================================================================

function [valueString] = readLabeledString( f, label )


frewind( f )  % start search at the begining of the file
oneLine = fgetl( f );

while( ~strncmp( label, oneLine, length( label ) ) & ~feof( f ) )
   oneLine = fgetl( f ); 
end

if( feof( f ) )
   fclose( f );   
   error( ['Fail to find label: ', label, ', in provided file'] )
end

[labelTokenized, valueTokenized] = strtok( oneLine, '=' );
valueString = valueTokenized(2:length( valueTokenized ));



%=================================================================================================

function [outResult] = existFile( filename )

outResult = exist( filename, 'file' ) == 2;  %  Test if file exists



