
% Remove any old data
!rm -f *.txt

% Loop over the input files
range = [ Inf , -Inf ];
avg = [ 0 , 0 , 0 ];
count = 0;
shift = [ 4.331456e+01 , 4.097030e+01 , 4.383432e+01 ];
for i=0:15

    % Get the file name
    fname = sprintf( 'snap_C09_200_000.%i.hdf5' , i );

    % Get the coordinates
    coord = double( h5read( fname , '/PartType0/Coordinates' )' );
    coord = coord - repmat( shift , size( coord , 1 ) , 1 );
    save Coordinates.txt -ascii -double -append coord
    
    % Adjust the range
    range(1) = min( range(1) , min(min(coord)) );
    range(2) = max( range(2) , max(max(coord)) );
    avg = avg + sum( coord , 1 );
    count = count + size( coord , 1 );
    
    % Get the smoothing lengths
    h = double( h5read( fname , '/PartType0/SmoothingLength' ) );
    save SmoothingLength.txt -ascii -double -append h
    
end

% Display some statistics
disp( sprintf( 'read %i particles' , count ) );
disp( sprintf( 'range of coords is [ %e , %e ]' , range(1) , range(2) ) );
disp( sprintf( 'avg position is [ %e , %e , %e ]' , avg(1)/count , avg(2)/count , avg(3)/count ) );
