%shMovieDRed.m
%
%	Show movie of potential and currents
%	shMovieDRed( filename, resolution [default = 4], frameStep [default = 10, at 10Khz that is 1 frame/ms] )


function shMovieDRed_nl( filename, resolution, frameStep )

%-- Please adjust maximum ---
maxV = 2.5;
%----------------------------

if( nargin <= 2), frameStep = 10;, end
if( nargin == 1), resolution= 4;, end

eval( ['load ',filename] )

status = mkdir('tempPot');
status = mkdir('tempCSD');

sideArray = 8;

steps = 8 * resolution;
[XI, YI] = meshgrid(1:(8-1)/(steps-1):8,1:(8-1)/(steps-1):8);
[X, Y] = meshgrid(1:8,1:8);

figure(1)
subplot(121)
h1 = image( rand( size( XI ) ) );
hTitle1 = title('*');
set( hTitle1, 'EraseMode', 'normal' )
set( h1, 'EraseMode', 'none' )
shading flat
axis('image')

subplot(122)
h2 = image( rand( size( XI ) ) );
hTitle2 = title('*');
set( hTitle2, 'EraseMode', 'normal' )
set( h2, 'EraseMode', 'none' )
colormap(hot);
shading flat
axis('image')


%--- Animation loop ---
L = [1 4 1; 4 -20 4; 1 4 1]; %laplacian

frameNum = 0;

for frame = 1:frameStep:size(data, 1)
   
   %--- interpolate array ---
   
   Z1 = toMatrix( data(frame, 2:65), sideArray );
   ZI1 = interp2( X, Y, Z1, XI, YI, 'cubic' );  
   
   current = conv2( Z1, L, 'same' ); 
   currentI = -interp2( X, Y, current, XI, YI, 'cubic' );
   
   
   %--- Update data in plot ---
   
   maxDV = maxV * 20;
   
   ZI1 = ( trunc11( ZI1/maxV ) + 1 ) / 2 * 64;
   currentI = ( trunc11( currentI/maxDV ) + 1 ) / 2 * 64;
   
   set( h1, 'cdata', ZI1 )
   set( h2, 'cdata', currentI )
   set( hTitle1, 'String', ['Potential. Time: ', num2str( data(frame, 1), '%04.1f' ), ' [ms]'] );
   set( hTitle2, 'String', 'CSD' );
   drawnow
   
   %--- Save frame ---
   frameNum = frameNum + 1;
   frameFilenameZ = ['tempPot\frame', int2str(frameNum), '.bmp'];
   frameFilenameI = ['tempCSD\frame', int2str(frameNum), '.bmp'];   
   eval( ['imwrite( ZI1, hot, ''', frameFilenameZ , ''', ''bmp'' ) '] )
   eval( ['imwrite( currentI, hot, ''', frameFilenameI , ''', ''bmp'' ) '] )
   
  
end

%----------------------------------
function [m] = toMatrix( vector, size )
m = zeros( size );

for lin = 1:size
   
   m(lin,1:size) = vector(1 + (lin-1)*size:size + (lin-1)*size);
   
end

%----------------------------------
function [A]= trunc11(A)

[N_A dim]= size(A);

for i=1:N_A
	for j=1:dim
		if (A(i,j)>1)
			A(i,j)= 1;
		else
			if (A(i,j)<-1)
				A(i,j)= -1;
			end
		end
	end
end