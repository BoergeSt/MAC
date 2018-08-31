function [ M,F ] = Apply_BCs_U( M,F,domain,mesh,pde)
I = speye(mesh.calculated.DOFs);

on_dirichlet_boundary_vec = domain.dirichlet_Boundary_u(mesh.calculated.X_u(mesh.calculated.on_bndry_u),...
                                                        mesh.calculated.Y_u(mesh.calculated.on_bndry_u));
d_bndry_u = mesh.calculated.on_bndry_u(on_dirichlet_boundary_vec==1);
n_bndry_u = mesh.calculated.on_bndry_u(~on_dirichlet_boundary_vec);

M(mesh.calculated.on_bndry_u,:) = I(mesh.calculated.on_bndry_u,:);
F(d_bndry_u) = pde.calculated.u0(mesh.calculated.X_u(d_bndry_u),mesh.calculated.Y_u(d_bndry_u));
F(n_bndry_u) = pde.calculated.un(mesh.calculated.X_u(n_bndry_u),mesh.calculated.Y_u(n_bndry_u));

for i = n_bndry_u
    normal = Get_Outer_Normal( domain.calculated.on_domain, mesh.calculated.X_u(i),mesh.calculated.Y_u(i),...
                                mesh.h*mesh.relative_tangential_length,mesh.number_normal_sample_vectors);
    M(i,i) = mesh.h*(sum(abs(normal)));
    M(i,i+(mesh.calculated.nx+1)*sign(normal(2))) =M(i,i+(mesh.calculated.nx+1)*sign(normal(2)))  -abs(normal(2))*mesh.h;
    M(i,i-sign(normal(1))) =M(i,i-sign(normal(1))) -abs(normal(1))*mesh.h;
end


M(mesh.calculated.on_ghost_u,1:mesh.calculated.nu) = 0;
%M(on_ghost_u,:)=0;
for i = 1:length(mesh.calculated.on_ghost_u)
   M(mesh.calculated.on_ghost_u(i),mesh.calculated.on_ghost_u(i)) = ...
       2/mesh.h^2*(1/mesh.calculated.ghost_dir_u(i,1)/mesh.calculated.ghost_dir_u(i,3) +...
                   1/mesh.calculated.ghost_dir_u(i,2)/mesh.calculated.ghost_dir_u(i,4));
   if abs(1-mesh.calculated.ghost_dir_u(i,1))<=eps%up
       M(mesh.calculated.on_ghost_u(i),mesh.calculated.on_ghost_u(i)+mesh.calculated.nx+1) = ...
           1/mesh.h^2*(-2/(1+mesh.calculated.ghost_dir_u(i,3)));
   else
       F(mesh.calculated.on_ghost_u(i)) =F(mesh.calculated.on_ghost_u(i))+...
           2/mesh.h^2/mesh.calculated.ghost_dir_u(i,1)/(mesh.calculated.ghost_dir_u(i,1)+mesh.calculated.ghost_dir_u(i,3))*...
           pde.calculated.u0(mesh.calculated.X_u(mesh.calculated.on_ghost_u(i)),mesh.calculated.Y_u(mesh.calculated.on_ghost_u(i))+...
           mesh.h*mesh.calculated.ghost_dir_u(i,1));
   end
   if abs(1-mesh.calculated.ghost_dir_u(i,3))<=eps%down
       M(mesh.calculated.on_ghost_u(i),mesh.calculated.on_ghost_u(i)-mesh.calculated.nx-1) = ...
           1/mesh.h^2*(-2/(1+mesh.calculated.ghost_dir_u(i,1)));
   else
       F(mesh.calculated.on_ghost_u(i)) =F(mesh.calculated.on_ghost_u(i))+ 2/mesh.h^2/mesh.calculated.ghost_dir_u(i,3)/...
           (mesh.calculated.ghost_dir_u(i,1)+mesh.calculated.ghost_dir_u(i,3))*pde.calculated.u0(mesh.calculated.X_u(mesh.calculated.on_ghost_u(i)),...
           mesh.calculated.Y_u(mesh.calculated.on_ghost_u(i))-mesh.h*mesh.calculated.ghost_dir_u(i,3));
   end
   if abs(1-mesh.calculated.ghost_dir_u(i,2))<=eps%right
       M(mesh.calculated.on_ghost_u(i),mesh.calculated.on_ghost_u(i)+1) = 1/mesh.h^2*(-2/(1+mesh.calculated.ghost_dir_u(i,4)));
   else
       F(mesh.calculated.on_ghost_u(i)) =F(mesh.calculated.on_ghost_u(i))+ 2/mesh.h^2/mesh.calculated.ghost_dir_u(i,2)/...
           (mesh.calculated.ghost_dir_u(i,2)+mesh.calculated.ghost_dir_u(i,4))*pde.calculated.u0(...
           mesh.calculated.X_u(mesh.calculated.on_ghost_u(i))+mesh.h*mesh.calculated.ghost_dir_u(i,2),...
           mesh.calculated.Y_u(mesh.calculated.on_ghost_u(i)));
       M(mesh.calculated.on_ghost_u(i),mesh.calculated.nu+mesh.calculated.nv+(1:mesh.calculated.np)) = 0;
   end
   if abs(1-mesh.calculated.ghost_dir_u(i,4))<=eps%left
       M(mesh.calculated.on_ghost_u(i),mesh.calculated.on_ghost_u(i)-1) = -2/mesh.h^2/(1+mesh.calculated.ghost_dir_u(i,2));
   else
       F(mesh.calculated.on_ghost_u(i)) =F(mesh.calculated.on_ghost_u(i))+...
           2/mesh.h^2/mesh.calculated.ghost_dir_u(i,4)/(mesh.calculated.ghost_dir_u(i,2)+...
           mesh.calculated.ghost_dir_u(i,4))*pde.calculated.u0(mesh.calculated.X_u(mesh.calculated.on_ghost_u(i))-...
           mesh.h*mesh.calculated.ghost_dir_u(i,4),mesh.calculated.Y_u(mesh.calculated.on_ghost_u(i)));
       M(mesh.calculated.on_ghost_u(i),mesh.calculated.nu+mesh.calculated.nv+(1:mesh.calculated.np)) = 0;
   end
end







end

