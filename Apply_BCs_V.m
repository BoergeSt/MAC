function [ M,F ] = Apply_BCs_V( M,F,domain,mesh,pde)
I = speye(mesh.calculated.DOFs);

on_dirichlet_boundary_vec = domain.dirichlet_Boundary_v(mesh.calculated.X_v(mesh.calculated.on_bndry_v),...
    mesh.calculated.Y_v(mesh.calculated.on_bndry_v));
d_bndry_v = mesh.calculated.on_bndry_v(on_dirichlet_boundary_vec==1);
n_bndry_v = mesh.calculated.on_bndry_v(~on_dirichlet_boundary_vec);

M(mesh.calculated.on_bndry_v+mesh.calculated.nu,:) = I(mesh.calculated.on_bndry_v+mesh.calculated.nu,:);

F(d_bndry_v+mesh.calculated.nu) = pde.calculated.v0(mesh.calculated.X_v(d_bndry_v),mesh.calculated.Y_v(d_bndry_v));
F(n_bndry_v+mesh.calculated.nu) = pde.calculated.vn(mesh.calculated.X_v(n_bndry_v),mesh.calculated.Y_v(n_bndry_v));

for i = n_bndry_v
    normal = Get_Outer_Normal( on_domain, mesh.calculated.X_v(i),mesh.calculated.Y_v(i), mesh.h*1e-3,201);
    M(i+mesh.calculated.nu,i+mesh.calculated.nu) = mesh.h*sum(abs(normal));
    M(i+mesh.calculated.nu,i+mesh.calculated.nu-mesh.calculated.nx*sign(normal(2))) =M(i+mesh.calculated.nu,i+mesh.calculated.nu-mesh.calculated.nx*sign(normal(2))) -abs(normal(2))*mesh.h;
    M(i+mesh.calculated.nu,i+mesh.calculated.nu-sign(normal(1))) =M(i+mesh.calculated.nu,i+mesh.calculated.nu-sign(normal(1))) -abs(normal(1))*mesh.h;
end



M(mesh.calculated.on_ghost_v+mesh.calculated.nu,mesh.calculated.nu+1:mesh.calculated.nu+1+(mesh.calculated.ny+1)*mesh.calculated.nx) = 0;
%M(mesh.calculated.on_ghost_v,:)=0;
for i = 1:length(mesh.calculated.on_ghost_v)
   M(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu,mesh.calculated.on_ghost_v(i)+mesh.calculated.nu) = 1/mesh.h^2*(2/mesh.calculated.ghost_dir_v(i,1)/mesh.calculated.ghost_dir_v(i,3) + 2/mesh.calculated.ghost_dir_v(i,2)/mesh.calculated.ghost_dir_v(i,4));
   if abs(1-mesh.calculated.ghost_dir_v(i,1))<=eps%up
       M(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu,mesh.calculated.on_ghost_v(i)+mesh.calculated.nu+mesh.calculated.nx) = 1/mesh.h^2*(-2/(1+mesh.calculated.ghost_dir_v(i,3)));
   else
       F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu) =F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu)+ 2/mesh.h^2/mesh.calculated.ghost_dir_v(i,1)/(mesh.calculated.ghost_dir_v(i,1)+mesh.calculated.ghost_dir_v(i,3))*pde.calculated.v0(mesh.calculated.X_v(mesh.calculated.on_ghost_v(i)),mesh.calculated.Y_v(mesh.calculated.on_ghost_v(i))+mesh.h*mesh.calculated.ghost_dir_v(i,1));
       M(mesh.calculated.on_ghost_v(i),(mesh.calculated.nx+1)*mesh.calculated.ny+(mesh.calculated.ny+1)*mesh.calculated.nx+(1:mesh.calculated.nx*mesh.calculated.ny)) = 0;
   end
   if abs(1-mesh.calculated.ghost_dir_v(i,3))<=eps%down
       M(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu,mesh.calculated.on_ghost_v(i)+mesh.calculated.nu-mesh.calculated.nx) = 1/mesh.h^2*(-2/(1+mesh.calculated.ghost_dir_v(i,1)));
   else
       F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu) =F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu)+ 2/mesh.h^2/mesh.calculated.ghost_dir_v(i,3)/(mesh.calculated.ghost_dir_v(i,1)+mesh.calculated.ghost_dir_v(i,3))*pde.calculated.v0(mesh.calculated.X_v(mesh.calculated.on_ghost_v(i)),mesh.calculated.Y_v(mesh.calculated.on_ghost_v(i))-mesh.h*mesh.calculated.ghost_dir_v(i,3));
       M(mesh.calculated.on_ghost_v(i),(mesh.calculated.nx+1)*mesh.calculated.ny+(mesh.calculated.ny+1)*mesh.calculated.nx+(1:mesh.calculated.nx*mesh.calculated.ny)) = 0;
   end
   if abs(1-mesh.calculated.ghost_dir_v(i,2))<=eps%right
       M(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu,mesh.calculated.on_ghost_v(i)+1+mesh.calculated.nu) = 1/mesh.h^2*(-2/(1+mesh.calculated.ghost_dir_v(i,4)));
   else
       F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu) =F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu)+ 2/mesh.h^2/mesh.calculated.ghost_dir_v(i,2)/(mesh.calculated.ghost_dir_v(i,2)+mesh.calculated.ghost_dir_v(i,4))*pde.calculated.v0(mesh.calculated.X_v(mesh.calculated.on_ghost_v(i))+mesh.h*mesh.calculated.ghost_dir_v(i,2),mesh.calculated.Y_v(mesh.calculated.on_ghost_v(i)));
   end
   if abs(1-mesh.calculated.ghost_dir_v(i,4))<=eps%left
       M(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu,mesh.calculated.on_ghost_v(i)+mesh.calculated.nu-1) = -2/mesh.h^2/(1+mesh.calculated.ghost_dir_v(i,2));
   else
       F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu) =F(mesh.calculated.on_ghost_v(i)+mesh.calculated.nu)+ 2/mesh.h^2/mesh.calculated.ghost_dir_v(i,4)/(mesh.calculated.ghost_dir_v(i,2)+mesh.calculated.ghost_dir_v(i,4))*pde.calculated.v0(mesh.calculated.X_v(mesh.calculated.on_ghost_v(i))-mesh.h*mesh.calculated.ghost_dir_v(i,4),mesh.calculated.Y_v(mesh.calculated.on_ghost_v(i)));
   end
end







end

