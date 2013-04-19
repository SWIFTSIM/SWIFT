
%% Make a combined plot for the CosmoVolume
gadget = importdata( 'data/CosmoVolume_fixed_Gadget2.totals' );
gadget = [ gadget(:,1) , gadget(:,6)/2.4e9*1000 ];
swift = importdata( 'data/CosmoVolume_fixed.totals');
swift = [ (1:size(swift,1))' , swift(:,end-1) ];
ncores = max( [ swift(:,1) ] );
nparts = 1841127; 

clf
subplot('position',[ 0.03 , 0.1 , 0.3 , 0.8 ]);
plot( swift(:,1) , swift(1,2) ./ swift(:,2) , '-k' , 'LineWidth' , 2 );  hold on;
plot( gadget(:,1) , gadget(1,2) ./ gadget(:,2) , '-.k' , 'LineWidth' , 2 );
% plot( gadget(:,1) , swift(1,2) ./ gadget(:,2) , '--k' , 'LineWidth' , 1 ); hold on;
text(ncores-0.2,gadget(1,2)/gadget(end,2)+0.2,sprintf('%.0f',min(gadget(:,2))),'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',12);
text(ncores-0.2,swift(1,2)/swift(end,2)+0.2,sprintf('%.0f',min(swift(:,2))),'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',12);
xlabel('nr. cores');
plot( [1,ncores] , [1,ncores] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Speedup Cosmological Volume');
axis([ 1 , ncores , 0 , ncores ]);

subplot('position',[ 0.38 0.1 , 0.6 , 0.8 ]);
plot( swift(:,1) , swift(1,2) ./ (swift(:,2).*swift(:,1)) , '-k' , 'LineWidth' , 2 ); hold on;
plot( gadget(:,1) , gadget(1,2) ./ (gadget(:,2).*gadget(:,1)) , '-.k' , 'LineWidth' , 2 );
xlabel('nr. cores');
legend('SWIFT','Gadget2','Location','SouthWest');
plot( [1,ncores] , [1,1] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Parallel Efficiency Cosmological Volume');
axis([ 1 , ncores , 0 , 1.1 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 12 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 12 4.4 ] );
print -depsc2 figures/CosmoVolume_scaling.eps
!epstopdf figures/CosmoVolume_scaling.eps 


%% Plot the timings for each part of the code.
swift = importdata( 'data/CosmoVolume_fixed.totals');
swift = [ (1:size(swift,1))' , swift(:,end-1).*(1:size(swift,1))' , sum(swift(:,[7,8]),2) , sum(swift(:,[5,6]),2) , swift(:,4) , swift(:,13) , sum(swift(:,[1:3]),2) ];
ncores = max( swift(:,1) );

clf
subplot( 'position' , [ 0.1 , 0.1 , 0.8 , 0.8 ] );
plot( swift(:,1) , swift(:,2:end)/1000 , 'LineWidth' , 1.4 );
legend( 'total' , 'pair' , 'self' , 'sort' , 'task' , 'int' , 'Location' , 'NorthWest' );
title('Total times');
axis([ 1 , ncores , 0 , 1.1*max(swift(:,2))/1000 ]);
xlabel('nr. cores');
ylabel('s');

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 8 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 8 4.4 ] );
print -depsc2 figures/CosmoVolume_times.eps
!epstopdf figures/CosmoVolume_times.eps 



%% Make a combined plot for the SodShock
swift = importdata( 'data/SodShock_fixed.totals');
swift = [ (1:size(swift,1))' , swift(:,end-1) ];
ncores = max( [ swift(:,1) ] );
nparts = 1841127; 

clf
subplot('position',[ 0.03 , 0.1 , 0.3 , 0.8 ]);
plot( swift(:,1) , swift(1,2) ./ swift(:,2) , '-k' , 'LineWidth' , 2 ); hold on;
text(ncores-0.2,swift(1,2)/swift(end,2)+0.2,sprintf('%.0f',min(swift(:,2))),'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',12);
xlabel('nr. cores');
plot( [1,ncores] , [1,ncores] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Speedup Sod-Shock');
axis([ 1 , ncores , 0 , ncores ]);

subplot('position',[ 0.38 0.1 , 0.6 , 0.8 ]);
plot( swift(:,1) , swift(1,2) ./ (swift(:,2).*swift(:,1)) , '-k' , 'LineWidth' , 2 ); hold on;
xlabel('nr. cores');
% legend('Sod-Shock','Location','SouthWest');
plot( [1,ncores] , [1,1] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Parallel Efficiency Sod-Shock');
axis([ 1 , ncores , 0 , 1.1 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 12 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 12 4.4 ] );
print -depsc2 figures/SodShock_scaling.eps
!epstopdf figures/SodShock_scaling.eps 



%% Make a combined plot for the SedovBlast
swift = importdata( 'data/SedovBlast_fixed.totals');
swift = [ (1:size(swift,1))' , swift(:,end-1) ];
ncores = max( [ swift(:,1) ] );
nparts = 1841127; 

clf
subplot('position',[ 0.03 , 0.1 , 0.3 , 0.8 ]);
plot( swift(:,1) , swift(1,2) ./ swift(:,2) , '-k' , 'LineWidth' , 2 ); hold on;
text(ncores-0.2,swift(1,2)/swift(end,2)+0.2,sprintf('%.0f',min(swift(:,2))),'HorizontalAlignment','Right','VerticalAlignment','Bottom','FontSize',12);
xlabel('nr. cores');
plot( [1,ncores] , [1,ncores] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Speedup Sedov Blast');
axis([ 1 , ncores , 0 , ncores ]);

subplot('position',[ 0.38 0.1 , 0.6 , 0.8 ]);
plot( swift(:,1) , swift(1,2) ./ (swift(:,2).*swift(:,1)) , '-k' , 'LineWidth' , 2 ); hold on;
xlabel('nr. cores');
% legend('Sedov Blast','Location','SouthWest');
plot( [1,ncores] , [1,1] , ':k' , 'LineWidth' , 1.4 );
hold off;
title('Parallel Efficiency Sedov Blast');
axis([ 1 , ncores , 0 , 1.1 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 12 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 12 4.4 ] );
print -depsc2 figures/SedovBlast_scaling.eps
!epstopdf figures/SedovBlast_scaling.eps 




%% Plot the Sod-shock density, pressure and velocity profiles
swift = importdata( 'data/SodShock_rhox.dat');
exact_rho = importdata( 'data/SodShock_exact_rho.dat' );
exact_P = importdata( 'data/SodShock_exact_P.dat' );
exact_v = importdata( 'data/SodShock_exact_v.dat' );

clf
subplot('position',[ 0.05 , 0.1 , 0.28 , 0.8 ]);
plot( 0.5*(swift(:,2)+swift(:,3)) - 0.5 , swift(:,9) , '-k' , 'LineWidth' , 2 ); hold on;
plot( exact_rho(:,1) , exact_rho(:,2) , ':k' , 'LineWidth' , 2 ); hold on;
xlabel('x'); ylabel('\rho');
hold off;
title('Sod-shock Density');
axis([ -0.25 0.25 0.2 1 ]);

subplot('position',[ 1/3+0.05 , 0.1 , 0.28 , 0.8 ]);
plot( 0.5*(swift(:,2)+swift(:,3)) - 0.5 , swift(:,10) , '-k' , 'LineWidth' , 2 ); hold on;
plot( exact_P(:,1) , exact_P(:,2) , ':k' , 'LineWidth' , 2 ); hold on;
xlabel('x'); ylabel('P');
hold off;
title('Sod-shock Pressure');
axis([ -0.25 0.25 0.1 1 ]);

subplot('position',[ 2/3+0.05 , 0.1 , 0.28 , 0.8 ]);
plot( 0.5*(swift(:,2)+swift(:,3)) - 0.5 , swift(:,5) , '-k' , 'LineWidth' , 2 ); hold on;
plot( exact_v(:,1) , exact_v(:,2) , ':k' , 'LineWidth' , 2 ); hold on;
xlabel('x'); ylabel('v_x');
hold off;
title('Sod-shock Velocity');
axis([ -0.25 0.25 -0.8 0.8 ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 12 4.2 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 12 4.2 ] );
print -depsc2 figures/SodShock_profile.eps
!epstopdf figures/SodShock_profile.eps 



%% Plot the Sedov blast density profile
swift_10 = importdata( 'data/SedovBlast_rdf_10.dat');
swift_15 = importdata( 'data/SedovBlast_rdf_15.dat');
swift_20 = importdata( 'data/SedovBlast_rdf_20.dat');
exact_10 = importdata( 'data/SedovBlast_exact_10.dat');
exact_15 = importdata( 'data/SedovBlast_exact_15.dat');
exact_20 = importdata( 'data/SedovBlast_exact_20.dat');

clf
subplot('position',[ 0.03 , 0.1 , 0.8 , 0.8 ]);
plot( exact_15(:,1) , exact_15(:,2) , ':k' , 'LineWidth' , 2 ); hold on;
plot( 0.5*(swift_15(:,2)+swift_15(:,3)) , swift_15(:,9) , 'ok' , 'LineWidth' , 2 );
plot( exact_10(:,1) , exact_10(:,2) , ':k' , 'LineWidth' , 2 ); hold on;
plot( 0.5*(swift_10(:,2)+swift_10(:,3)) , swift_10(:,9) , 'xk' , 'LineWidth' , 2 );
plot( exact_20(:,1) , exact_20(:,2) , ':k' , 'LineWidth' , 2 ); hold on;
plot( 0.5*(swift_20(:,2)+swift_20(:,3)) , swift_20(:,9) , '^k' , 'LineWidth' , 2 );
xlabel('r');
ylabel('\rho');
hold off;
title('Density Sedov blast');
axis([ 0 , 2.5 , 0 , max(exact_20(:,2)) ]);

% Print this plot
set( gcf , 'PaperSize' , 2.3*[ 8 4.5 ] );
set( gcf , 'PaperPosition' , 2.3*[ 0.25 0.25 8 4.4 ] );
print -depsc2 figures/SedovBlast_density.eps
!epstopdf figures/SedovBlast_density.eps 

