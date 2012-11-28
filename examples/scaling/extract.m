
% Remove any old data
!rm -f *.txt

% Loop over the input files
range = [ Inf , -Inf ];
avg = [ 0 , 0 , 0 ];
count = 0;
shift = [ 0 0 0 ];

% Get the file name
fname = 'snap_023_z000p503.hdf5';

% Get the coordinates
coord = double( h5read( fname , '/PartType0/Coordinates' )' );
coord = coord - repmat( shift , size( coord , 1 ) , 1 );

% Get the smoothing lengths
h = double( h5read( fname , '/PartType0/SmoothingLength' ) );

% Remove entries with too large smoothing lengths
% ind = (h < 150);
% coord = coord(ind,:);
% h = h(ind);

% Save the data
save Coordinates.txt -ascii -double -append coord
save SmoothingLength.txt -ascii -double -append h

% Get some statistics
count = size( coord , 1 );
avg = sum( coord , 1 ) / count;
range(1) = min( range(1) , min(min(coord)) );
range(2) = max( range(2) , max(max(coord)) );
    
% Display some statistics
disp( sprintf( 'read %i particles' , count ) );
disp( sprintf( 'range of coords is [ %e , %e ]' , range(1) , range(2) ) );
disp( sprintf( 'range of h is [ %e , %e ]' , min(h) , max(h) ) );
disp( sprintf( 'avg position is [ %e , %e , %e ]' , avg(1) , avg(2) , avg(3) ) );
