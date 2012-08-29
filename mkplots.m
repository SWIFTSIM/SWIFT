
%% Scaling and efficiency plots for the first test
test = importdata( 'test_nosort.totals');
ncores = size( test , 1 );

clf
subplot('position',[ 0.05 , 0.1 , 0.4 , 0.8 ]);
plot( 1:ncores , test(1,10) ./ test(:,10) , '-k' , 'LineWidth' , 2 ); hold on;
xlabel('nr. cores');
plot( [1,ncores] , [1,ncores] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Speedup');
axis([ 1 , ncores , 0 , ncores ]);

subplot('position',[ 0.52 0.1 , 0.4 , 0.8 ]);
plot( 1:ncores , test(1,10) ./ test(:,10) ./ (1:ncores)' , '-k' , 'LineWidth' , 2 ); hold on;
plot( [1,ncores] , [1,1] , ':k' , 'LineWidth' , 1.4 );
text(4*ncores/5,1,sprintf('%.2f',min(test(:,10))),'BackgroundColor',[1,1,1],'FontSize',12);
xlabel('nr. cores');
hold off;
title('Efficiency');
axis([ 1 , ncores , 0 , 1.2 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 8.5 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 8 4 ] );
print -depsc2 test.eps
!epstopdf test.eps 


%% Scaling and efficiency plots for the second test
test2 = importdata( 'test2.totals');
ncores = size( test2 , 1 );

clf
subplot('position',[ 0.05 , 0.1 , 0.4 , 0.8 ]);
plot( 1:ncores , test2(1,10) ./ test2(:,10) , '-k' , 'LineWidth' , 2 ); hold on;
xlabel('nr. cores');
plot( [1,ncores] , [1,ncores] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Speedup');
axis([ 1 , ncores , 0 , ncores ]);

subplot('position',[ 0.52 0.1 , 0.4 , 0.8 ]);
plot( 1:ncores , test2(1,10) ./ test2(:,10) ./ (1:ncores)' , '-k' , 'LineWidth' , 2 ); hold on;
plot( [1,ncores] , [1,1] , ':k' , 'LineWidth' , 1.4 );
text(4*ncores/5,1,sprintf('%.2f',min(test2(:,10))),'BackgroundColor',[1,1,1],'FontSize',12);
xlabel('nr. cores');
hold off;
title('Efficiency');
axis([ 1 , ncores , 0 , 1.2 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 8.5 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 8 4 ] );
print -depsc2 test2.eps
!epstopdf test2.eps 


%% Components of the first test
test = importdata( 'test_nosort.totals');
ncores = size( test , 1 );
cols = [ 1 , 2 , 3 , 5 ];

clf
subplot('position',[ 0.1 , 0.1 , 0.8 , 0.8 ]);
plot( 1:ncores , test(:,cols) , 'LineWidth' , 2 ); hold on;
legend( 'sort' , 'self' , 'pair' , 'get task' , 'Location' , 'NorthWest' );
xlabel('nr. cores');
hold off;
title('ms per task type');
axis([ 1 , ncores , 0 , max(max(test(:,cols)))*1.1 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 4.5 4 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 4 4 ] );
print -depsc2 parts.eps
!epstopdf parts.eps 


%% Components of the second test
test2 = importdata( 'test2.totals');
ncores = size( test2 , 1 );
cols = [ 1 , 2 , 3 , 5 ];

clf
subplot('position',[ 0.1 , 0.1 , 0.8 , 0.8 ]);
plot( 1:ncores , test2(:,cols) , 'LineWidth' , 2 ); hold on;
legend( 'sort' , 'self' , 'pair' , 'get task' , 'Location' , 'NorthWest' );
xlabel('nr. cores');
hold off;
title('ms per task type');
axis([ 1 , ncores , 0 , max(max(test2(:,cols)))*1.1 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 4.5 4 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 4 4 ] );
print -depsc2 parts2.eps
!epstopdf parts2.eps 

