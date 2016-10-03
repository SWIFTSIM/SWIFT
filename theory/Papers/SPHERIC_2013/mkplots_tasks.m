
%% Some global stuff
tpms = 1/2.4e6;
ms_max = 82;

%% Plot the task timelines for tasks allocation
% Load the data
tasks = importdata( 'SedovBlast_fixed_16.dump' );
tasks(:,6) = ( tasks(:,6) - tasks(:,5) ) * tpms;
start = min( tasks(:,5) );
tasks(:,5) = ( tasks(:,5) - start ) * tpms;
% for k=0:max(tasks(:,1))
%     temp = tasks( tasks(:,1)==k , : );
%     start = min( temp(:,5) );
%     temp(:,5) = ( temp(:,5) - start ) * tpms;
%     tasks( tasks(:,1)==k , :) = temp;
% end
tasks_sorts = tasks( tasks(:,2)==1 , : );
tasks_pairs = tasks( tasks(:,2)==3 | ( tasks(:,2)==4 & tasks(:,4)==0 ) , : );
tasks_selfs = tasks( tasks(:,2)==2 | ( tasks(:,2)==4 & tasks(:,4)==1 ) , : );
tasks_ghosts = tasks( tasks(:,2)==5 , : );
tasks_kick2s = tasks( tasks(:,2)==6 , : );
nr_cores = max( tasks(:,1) ) + 1;

% Init the plot
clf;
subplot('position',[ 0.05 , 0.1 , 0.9 , 0.8 ]);
hold on;

% Plot the pairs
for k=1:size(tasks_pairs,1)
    rectangle( 'Position' , [ tasks_pairs(k,5) , tasks_pairs(k,1)+0.5 , tasks_pairs(k,6) , 1 ] , ...
        'EdgeColor' , [ 0 0.8 0 ] , 'LineWidth' , 1 , 'FaceColor' , [ 0 1 0 ] );
end

% Plot the selfs
for k=1:size(tasks_selfs,1)
    rectangle( 'Position' , [ tasks_selfs(k,5) , tasks_selfs(k,1)+0.5 , tasks_selfs(k,6) , 1 ] , ...
        'EdgeColor' , [ 0 0 0.8 ] , 'LineWidth' , 1 , 'FaceColor' , [ 0 0 1 ] );
end

% Plot the sorts
for k=1:size(tasks_sorts,1)
    rectangle( 'Position' , [ tasks_sorts(k,5) , tasks_sorts(k,1)+0.5 , tasks_sorts(k,6) , 1 ] , ...
        'EdgeColor' , [ 0.8 0 0 ] , 'LineWidth' , 1 , 'FaceColor' , [ 1 0 0 ] );
end

% Plot the ghosts
for k=1:size(tasks_ghosts,1)
    rectangle( 'Position' , [ tasks_ghosts(k,5) , tasks_ghosts(k,1)+0.5 , tasks_ghosts(k,6) , 1 ] , ...
        'EdgeColor' , [ 0.9 0.9 0 ] , 'LineWidth' , 1 , 'FaceColor' , [ 1 1 0 ] );
end

% Plot the kick2s
for k=1:size(tasks_kick2s,1)
    rectangle( 'Position' , [ tasks_kick2s(k,5) , tasks_kick2s(k,1)+0.5 , tasks_kick2s(k,6) , 1 ] , ...
        'EdgeColor' , [ 1 0.5 0 ] , 'LineWidth' , 1 , 'FaceColor' , [ 1 0.5 0 ] );
end

% Set the axes and stuff.
hold off;
xlabel('time (ms)');
ylabel('core ID');
set(gca,'YTick',1:(max(tasks(:,1))+1))
title('SWIFT tasks');
axis([ 0 , ms_max , 0.5 , max(tasks(:,1))+1.5 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 16 4 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 16 4 ] );
print -depsc2 tasks_dynamic.eps
!epstopdf tasks_dynamic.eps 



%% Make a fake task ordering a la OpenMP.

% Init the core timers
timers = zeros( nr_cores , 1 );
for k=1:nr_cores
    timers(k) = min( tasks( tasks(:,1)==k , 5 ) );
end;
delta = 0.0035;

% Loop through the sort tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if ( tasks(k,2)==1 )
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end

% Re-set the timers
timers = ones( nr_cores , 1 ) * max( timers );

% Loop through the self density tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if tasks(k,2)==2 && tasks(k,3)==1
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end

% Re-set the timers
timers = ones( nr_cores , 1 ) * max( timers );

% Loop through the pair density tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if tasks(k,2)==3 && tasks(k,3)==1
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end

% Re-set the timers
timers = ones( nr_cores , 1 ) * max( timers );

% Loop through the ghost tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if tasks(k,2)==5
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end

% Re-set the timers
timers = ones( nr_cores , 1 ) * max( timers );

% Loop through the self force tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if tasks(k,2)==2 && tasks(k,3)==2
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end

% Re-set the timers
timers = ones( nr_cores , 1 ) * max( timers );

% Loop through the pair force tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if tasks(k,2)==3 && tasks(k,3)==2
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end

% Re-set the timers
timers = ones( nr_cores , 1 ) * max( timers );

% Loop through the kick tasks and put them at the front of the shortest queue
for k=1:size(tasks,1)
    if tasks(k,2)==6
        [ tic, ind ] = min( timers );
        tasks(k,5) = tic;
        tasks(k,1) = ind-1;
        timers(ind) = timers(ind) + tasks(k,6) + delta;
    end
end



% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 20 4 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 20 4 ] );
print -depsc2 tasks_static.eps
!epstopdf tasks_static.eps 

