clearvars -except display;

%*******************General Domain Settings***********************
domain.lx = 1;    % length of the domain
domain.ly = 1;    % height of t
domain.on_obstacle = @(x,y) 0*x; %indicatorfunction of the Obstacle

domain.dirichlet_Boundary_u = @(x,y) 0*x+1; %indicatorfunctions of the dirichlet boundary
domain.dirichlet_Boundary_v = @(x,y) 0*x+1; 


%************************Mesh Settings****************************
mesh.h = 1/20;  % distance between closest datapoints of the same type
mesh.rounding_precision = 10; %precision in digits after , of on domain rounding (used for p)
mesh.bndry_iteration_number = 50; %therefore 2^-x error. Be carefull with numbers under 50 since this could lead to wrong boundry evaluation
mesh.number_normal_sample_vectors = 200;
mesh.relative_tangential_length = 1e-3;

%********************PDE Specific Settings************************
%Right hand side:
pde.fut = @(x,y,t) 0*x;     
pde.fvt = @(x,y,t) 0*x;
%Dirichlet Boundary Condition.
pde.u0t = @(x,y,t) domain.lx^3/6*x.*(domain.lx-x);     
pde.v0t = @(x,y,t) 0*x;     
%Neumann Boundary Condition
pde.unt = @(x,y,t) 0*x;
pde.vnt = @(x,y,t) 0*x;
%Dynamic settings
pde.dynamic.dynamic = false; %Dynamic Calculation
pde.dynamic.u_start = @(x,y) 0*x+1;   %Initial state
pde.dynamic.v_start = @(x,y) 0*x;   
pde.dynamic.time_Dependent = false;   %Toggles Reevaluation of the BC and RHS in each timestep.


integrator.T = 10; % Endtime of the simulation
integrator.nT = 1e4; % number of timesteps
integrator.plot_every = 1e2; % updates the figure every x time_steps
integrator.theta = 1; %Specifies which theta scheme used for time integration



%*********************General Display Settings********************
display.set_Figure_Position = true; %Opens the figure at specific location
display.figure_Position = [228 296.5000 1471 571];  %sets the location
display.nx = 40;    %sets the number of interpolationpoints for the amplitude and arrowplot
display.ny = 40;    
display.display_factor = 1; %sets the factor of the arrowlength for the arrowplot
display.mesh = true; %whether the mesh is plotted (time-intensive for small h)
display.full_HD_for_video = true; %sets the figuresize to the maximal size for mpeg4 lvl2
display.monitor_4k = true; %set this to true if you use a 4k monitor (since matlab always calculates positon of pixels for fullHD)


display.fix_Axis = false; %Fixes the Axes in the plot. (recomended for dynamic calculations)
display.range.u = [0,2];
display.range.v = [-2,2];
display.range.p = [-10,10];
display.range.amplitude = [0,4];
display.range.arrowplot = [0,domain.lx;0,domain.ly];



%************************Video Settings***************************
video.record = false;   %toggles the recording of the main figure in order to make a replayable movie
video.name = 'out'; %the name of the video
video.profile = 'MPEG-4';
video.FrameRate = 60;




%************************Debugging and Errorhandling**************
debugging.stop_video_on_rank_deficit = true;
debugging.rankcheck=false; %checks the rank of the assembled matrix (Carefull! This assembles the full matrix
                           % Therefore it needs a lot of memory and time if n is to high! feasable up to 1000 DOFs)
debugging.cleanup = true; %cleans up unnecessary/no longer needed variables





