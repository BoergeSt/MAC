function calculated = Generate_Mesh(domain,mesh,display)
calculated.nx = round(domain.lx/mesh.h);
calculated.ny = round(domain.ly/mesh.h);
% Attention: nx and ny are the number of blocks.
% i.e. ncalculated.X_u = nx+1, ncalculated.Y_u = ny+1

%Check, if grid size h fits the Rectancle:
if calculated.nx ~= domain.lx/mesh.h || calculated.ny ~= domain.ly/mesh.h
    error('Choose lx,ly,h such that the number of elements is integer');
end




calculated.nu = (calculated.nx+1)*calculated.ny;
calculated.nv = calculated.nx*(calculated.ny+1);
calculated.np = calculated.nx*calculated.ny;
calculated.DOFs = calculated.nu+calculated.nv+calculated.np;




[X,Y] = ndgrid(0:mesh.h:domain.lx,0:mesh.h:domain.ly);
[calculated.X_p,calculated.Y_p] = ndgrid(mesh.h/2:mesh.h:domain.lx-mesh.h/2,mesh.h/2:mesh.h:domain.ly-mesh.h/2);
[calculated.X_u,calculated.Y_u] = ndgrid(0:mesh.h:domain.lx,mesh.h/2:mesh.h:domain.ly-mesh.h/2);
[calculated.X_v,calculated.Y_v] = ndgrid(mesh.h/2:mesh.h:domain.lx-mesh.h/2,0:mesh.h:domain.ly);

[calculated.XView,calculated.YView] = ndgrid(0:display.calculated.hx:domain.lx, ...
                                                            0:display.calculated.hy:domain.ly);




domain_u = domain.calculated.on_domain(calculated.X_u,calculated.Y_u);
domain_v = domain.calculated.on_domain(calculated.X_v,calculated.Y_v);

% P RESTRICTION
% at least one neighbor
%domain_p = on_domain(calculated.X_p+h/2,calculated.Y_p) | on_domain(calculated.X_p-h/2,calculated.Y_p) |on_domain(calculated.X_p,calculated.Y_p+h/2) |on_domain(calculated.X_p,calculated.Y_p-h/2);
% all neighbors and rounded
domain_p =  domain.calculated.on_domain(round(calculated.X_p+mesh.h/2,mesh.rounding_precision),calculated.Y_p) & ...
            domain.calculated.on_domain(round(calculated.X_p-mesh.h/2,mesh.rounding_precision),calculated.Y_p)  & ...
            domain.calculated.on_domain(calculated.X_p,round(calculated.Y_p+mesh.h/2,mesh.rounding_precision))  & ...
            domain.calculated.on_domain(calculated.X_p,round(calculated.Y_p-mesh.h/2,mesh.rounding_precision));

% all neighbors and not rounded
%domain_p = on_domain(calculated.X_p+h/2,calculated.Y_p) & on_domain(calculated.X_p-h/2,calculated.Y_p)  & on_domain(calculated.X_p,calculated.Y_p+h/2)  & on_domain(calculated.X_p,calculated.Y_p-h/2); 



calculated.domain_vec = [domain_u(:);domain_v(:);domain_p(:)];

index_u = 1:numel(calculated.X_u);
index_v = 1:numel(calculated.X_v);
index_p = 1:numel(calculated.X_p);

obstacle_p = index_p(~domain_p);

on_border_u = domain_u & ~(domain.calculated.on_domain(calculated.X_u+mesh.h,calculated.Y_u) & ...
                           domain.calculated.on_domain(calculated.X_u-mesh.h,calculated.Y_u) & ...
                           domain.calculated.on_domain(calculated.X_u,calculated.Y_u+mesh.h) & ...
                           domain.calculated.on_domain(calculated.X_u,calculated.Y_u-mesh.h));
                       
on_border_v = domain_v & ~(domain.calculated.on_domain(calculated.X_v+mesh.h,calculated.Y_v) & ...
                           domain.calculated.on_domain(calculated.X_v-mesh.h,calculated.Y_v) & ...
                           domain.calculated.on_domain(calculated.X_v,calculated.Y_v+mesh.h) & ...
                           domain.calculated.on_domain(calculated.X_v,calculated.Y_v-mesh.h));

calculated.on_bndry_u = on_border_u & ~(domain.calculated.on_domain(calculated.X_u+10*eps,calculated.Y_u) & ...
                             domain.calculated.on_domain(calculated.X_u-10*eps,calculated.Y_u) & ...
                             domain.calculated.on_domain(calculated.X_u,calculated.Y_u+10*eps) & ...
                             domain.calculated.on_domain(calculated.X_u,calculated.Y_u-10*eps));
                         
calculated.on_bndry_v = on_border_v & ~(domain.calculated.on_domain(calculated.X_v+10*eps,calculated.Y_v) & ...
                             domain.calculated.on_domain(calculated.X_v-10*eps,calculated.Y_v) & ...
                             domain.calculated.on_domain(calculated.X_v,calculated.Y_v+10*eps) & ...
                             domain.calculated.on_domain(calculated.X_v,calculated.Y_v-10*eps));

on_border_u = on_border_u(:);
on_border_v = on_border_v(:);

calculated.on_bndry_u = calculated.on_bndry_u(:);
calculated.on_bndry_v = calculated.on_bndry_v(:);

calculated.on_ghost_u = on_border_u & ~calculated.on_bndry_u;
calculated.on_ghost_v = on_border_v & ~calculated.on_bndry_v;

calculated.on_bndry_u = index_u(calculated.on_bndry_u);
calculated.on_bndry_v = index_v(calculated.on_bndry_v);

calculated.on_ghost_u = index_u(calculated.on_ghost_u);
calculated.on_ghost_v = index_v(calculated.on_ghost_v);

calculated.ghost_dir_u = zeros(length(calculated.on_ghost_u),4); %[up right down left]
calculated.ghost_dir_v = zeros(length(calculated.on_ghost_v),4);

for i =1:length(calculated.on_ghost_u)
    x = calculated.X_u(calculated.on_ghost_u(i));
    y = calculated.Y_u(calculated.on_ghost_u(i));
    calculated.ghost_dir_u(i,1) = Find_Border_Distance(domain.calculated.on_domain,x,y,0,mesh.h,mesh.bndry_iteration_number);%up
    calculated.ghost_dir_u(i,2) = Find_Border_Distance(domain.calculated.on_domain,x,y,mesh.h,0,mesh.bndry_iteration_number);%right
    calculated.ghost_dir_u(i,3) = Find_Border_Distance(domain.calculated.on_domain,x,y,0,-mesh.h,mesh.bndry_iteration_number);%down
    calculated.ghost_dir_u(i,4) = Find_Border_Distance(domain.calculated.on_domain,x,y,-mesh.h,0,mesh.bndry_iteration_number);%left
end

delete_vector = [];
for i = 1:length(calculated.on_ghost_u)
    if 1-min(calculated.ghost_dir_u(i,:))<2^(-mesh.bndry_iteration_number+1)
        delete_vector = [delete_vector;i];
    end
end
calculated.ghost_dir_u(delete_vector,:)=[];
calculated.on_ghost_u(delete_vector)=[];

for i =1:length(calculated.on_ghost_v)
    x = calculated.X_v(calculated.on_ghost_v(i));
    y = calculated.Y_v(calculated.on_ghost_v(i));
    calculated.ghost_dir_v(i,1) = Find_Border_Distance(domain.calculated.on_domain,x,y,0,mesh.h,mesh.bndry_iteration_number);
    calculated.ghost_dir_v(i,2) = Find_Border_Distance(domain.calculated.on_domain,x,y,mesh.h,0,mesh.bndry_iteration_number);
    calculated.ghost_dir_v(i,3) = Find_Border_Distance(domain.calculated.on_domain,x,y,0,-mesh.h,mesh.bndry_iteration_number);
    calculated.ghost_dir_v(i,4) = Find_Border_Distance(domain.calculated.on_domain,x,y,-mesh.h,0,mesh.bndry_iteration_number);
end

delete_vector = [];
for i = 1:length(calculated.on_ghost_v)
    if 1-min(calculated.ghost_dir_v(i,:))<2^(-mesh.bndry_iteration_number+1)
        delete_vector = [delete_vector;i];
    end
end
calculated.ghost_dir_v(delete_vector,:)=[];
calculated.on_ghost_v(delete_vector)=[];


if display.mesh
    figure(display.calculated.figure)
    subplot(2,3,6)



    axis equal

    plot(X,Y,'k-',X',Y','k-');
    hold on
    XX_p = calculated.X_p(domain_p);
    YY_p = calculated.Y_p(domain_p);
    XX_u = calculated.X_u(domain_u);
    YY_u = calculated.Y_u(domain_u);
    Xcalculated.X_v = calculated.X_v(domain_v);
    YY_v = calculated.Y_v(domain_v);

    plot(XX_p,YY_p,'g +','DisplayName','pressure');
    plot(XX_u,YY_u,'r o','DisplayName','x-velocity');
    plot(Xcalculated.X_v,YY_v,'b o','DisplayName','y-velocity');
    plot(calculated.X_u(calculated.on_bndry_u),calculated.Y_u(calculated.on_bndry_u),'r *');
    plot(calculated.X_v(calculated.on_bndry_v),calculated.Y_v(calculated.on_bndry_v),'b *');
    plot(calculated.X_u(calculated.on_ghost_u),calculated.Y_u(calculated.on_ghost_u),'k *');
    plot(calculated.X_v(calculated.on_ghost_v),calculated.Y_v(calculated.on_ghost_v),'k *');

    hold off
    title('Mesh with velocity and pressure nodes')
    gcaExpandable
end
