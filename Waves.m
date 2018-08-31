Init_MAC

pde.dynamic.dynamic = 1;

integrator.nT = 1e4;%1000000;
integrator.T = 10;
integrator.plot_every = 10;%10000;


mesh.h = 1/20;

domain.dirichlet_Boundary_u = @(x,y) 0*x+1-near(x,domain.lx);

pde.dynamic.time_Dependent = 1;

f = 1;
pde.v0t = @(x,y,t) 6/sqrt(2)*sin(f*pi*t)*y.*(domain.ly-y).*near(x,0);
%u0t = @(x,y,t) 3*sqrt(4-sin(pi*t)^2)*y.*(ly-y).*(abs(x)<eps)  + 6*y.*(ly-y).*(abs(x-lx)<eps) ;%parabol BC
pde.u0t = @(x,y,t) 6*sqrt(1-1/2*sin(f*pi*t)^2)*y.*(domain.ly-y).*near(x,0);  %für neumann

display.fix_Axis = true;
display.range.u = [0,1.5];
display.range.v = [-1,1];
display.range.p = [-40,40];
display.range.amplitude = [0,2];

video.record = 1;

Solve_MAC