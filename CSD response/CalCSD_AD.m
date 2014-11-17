%calCSD.m
%
%	Show 8 x 8 response e of potential and currents
%	CalCSD( inFilename [MAT file] )

function [csd_data]=calCSD_nl(data)

disp( '[calCSD] Loading data...' )
% eval( ['load ',inFilename] )
% outFilename = [inFilename,'_csd'];

sideArray = 8;

%--- Animation loop ---
%---insert by ken
S = [0 1/8 0; 1/8 1/2 1/8; 0 1/8 0];
L = [0 1 0; 1 -4 1; 0 1 0];
%L = [1 4 1; 4 -20 4; 1 4 1]; %laplacian
%---insert end

disp( '[calCSD] Calc CSD...' )
for frame = 1:size(data, 1)
   
   %--- interpolate array ---
   
   Z1 = toMatrix( data(frame, 2:65), sideArray );  
   
   S1 = conv2( Z1, S, 'same'); 
   current = -conv2( grow( S1 ), L, 'valid');
   %current = conv2( grow( Z1 ), L, 'valid');
   
   csd_data(frame, 2:65) = reshape(current.', 1, 64);  
end

csd_data(:, 1) = data(:, 1); 

disp( '[calCSD] Saving data...' )
% eval(['save ', outFilename])
    

%----------------------------------
function [m] = toMatrix( vector, size )
m = zeros( size );

for lin = 1:size
   
   m(lin,1:size) = vector(1 + (lin-1)*size:size + (lin-1)*size);
   
end

%----------------------------------