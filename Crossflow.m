Init_MAC

domain.lx = 6;
domain.ly = 2;

pde.dynamic.dynamic = 1;


integrator.nT = 1e3;%1000000;
integrator.T = 10;
integrator.plot_every = 1;%10000;


mesh.h = 1/20;


on_obstacle1 = @(x,y) x<2 & y<1;
on_obstacle2 = @(x,y) x>3 & y<1;



crossflow_inlet = @(t) 0.4*sin(pi*t);
domain.on_obstacle = @(x,y) on_obstacle1(x,y)|on_obstacle2(x,y);


pde.dynamic.time_Dependent = 1;
pde.u0t = @(x,y,t) -4*(y-1).*(y-2).*(abs(x)<eps)-4*(1+3/2*crossflow_inlet(t))*(y-1).*(y-2).*(abs(x-6)<eps);
pde.v0t = @(x,y,t) crossflow_inlet(t)*(abs(y)<eps);

pde.dynamic.u_start = @(x,y) 0*x;
pde.dynamicv_start = @(x,y) 0*x;

video.record = 0;

display.fix_Axis = 1;

Solve_MAC
