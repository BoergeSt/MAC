Init_MAC

mesh.h = 1/50;

o1 = 0;

on_obstacle_1 = @(x,y) (x-.3).^2+(y-.3).^2<.1^2;

domain.on_obstacle = @(x,y) (x+y>1)| o1*on_obstacle_1(x,y) ;



pde.u0t = @(x,y,t) 6*y.*(domain.ly-y).*near(x,0);%parabol BC
pde.v0t = @(x,y,t) -6*x.*(domain.lx-x).*near(y,0);

Solve_MAC
