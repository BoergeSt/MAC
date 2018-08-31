Init_MAC

domain.lx = 3;
domain.ly = 1;

mesh.h = 1/40;

o1 = 1; 
o2 = 1;
o3 = 1;


on_obstacle1 = @(x,y) (x-.5).^2+(y-.5).^2 < .05^2; 


on_obstacle2 = @(x,y) (x-1.2).^2+(y-.5).^2 < .1^2; 

on_obstacle3 = @(x,y) (x-2).^2+(y-.5).^2 < .2^2; 

domain.on_obstacle = @(x,y) o1*on_obstacle1(x,y) | o2*on_obstacle2(x,y) | o3*on_obstacle3(x,y);

domain.dirichlet_Boundary_u = @(x,y) 1-near(x,domain.lx);


pde.u0t = @(x,y,t) 6*y.*(domain.ly-y).*near(x,0);%parabol DBC on the left side
pde.unt = @(x,y,t) 0*x; %homogen neumann on the left side

pde.v0t = @(x,y,t) 0*x;

Solve_MAC
