function [ VxV ] = Create_VxV_Block( mesh )

e = ones(mesh.calculated.nx+1,1);
A_11 = spdiags([-e, 4.*e, -e], -1:1 ,mesh.calculated.nx,mesh.calculated.nx);
A_p1 = kron(speye(mesh.calculated.ny+1,mesh.calculated.ny+1),A_11);
A_12 = kron(spdiags([-ones(mesh.calculated.ny+1,1),-ones(mesh.calculated.ny+1,1)],[-1 1],...
                     mesh.calculated.ny+1,mesh.calculated.ny+1),speye(mesh.calculated.nx,mesh.calculated.nx));
VxV = 1/mesh.h^2 * (A_p1 + A_12);

end

