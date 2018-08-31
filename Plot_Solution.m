figure(display.calculated.figure)

% Velocity u
subplot(2,3,1)
surf(mesh.calculated.X_u',mesh.calculated.Y_u',U);
if display.fix_Axis
    zlim(display.range.u)
end
title('velocity u (in x - direction)')
gcaExpandable

% Velocity v
subplot(2,3,4)
surf(mesh.calculated.X_v',mesh.calculated.Y_v',V)
if display.fix_Axis
    zlim(display.range.v)
end
title('velocity v (in y - direction)')
gcaExpandable

% Pressure p
subplot(2,3,2)
surf(mesh.calculated.X_p',mesh.calculated.Y_p',P)
if display.fix_Axis
    zlim(display.range.p)
end
title('pressure p')
gcaExpandable


% Flow with arrows
    % interpolating with splines to achieve values on the gridpoints of 
    % [mesh.calculated.XView,mesh.calculated.YView]
subplot(2,3,3)
pu = interp2(mesh.calculated.X_u',mesh.calculated.Y_u',full(U),mesh.calculated.XView',mesh.calculated.YView','spline');
pv = interp2(mesh.calculated.X_v',mesh.calculated.Y_v',full(V),mesh.calculated.XView',mesh.calculated.YView','spline');
quiver(mesh.calculated.XView',mesh.calculated.YView',pu,pv,display.display_factor,'b');
title('flow')
if display.fix_Axis
    %todo
end
gcaExpandable

subplot(2,3,5)
surf(mesh.calculated.XView',mesh.calculated.YView',sqrt(pu.^2+pv.^2))
title('velocity magnitude |v|')
if display.fix_Axis
    zlim(display.range.amplitude)
end
gcaExpandable


if(pde.dynamic.dynamic)
    disp(['timestep ',int2str(time_step)])
end