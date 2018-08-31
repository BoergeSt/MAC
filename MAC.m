
Init_MAC;

mesh.h = 1/20;

o1 = 1; 
o2 = 0;
o3 = 0;
o4 = 0;

on_obstacle1 = @(x,y) 0*x+(x-.5).^2+(y-.5).^2 < .2^2; %circle

on_obstacle2 = @(x,y) 0*x+x>0.85 & y < 0.5; %step right,bottom

on_obstacle3 = @(x,y) 0*x+x<0.15 & y>0.5; %step left, top

on_obstacle4 = @(x,y) (x.^2+y.^2 < .2^2 )|((x-domain.lx).^2+(y-domain.ly).^2 < .2^2 ) ; %circle-corners

domain.on_obstacle = @(x,y) o1*on_obstacle1(x,y) | o2*on_obstacle2(x,y) | o3*on_obstacle3(x,y)| o4*on_obstacle4(x,y);


pde.u0t = @(x,y,t) 0*x+(1+o3)*near(x,0) + (1+o2)*near(x,domain.lx); %uniform Boundry Condition with no-slip
pde.v0t = @(x,y,t) 0*x;

Solve_MAC

