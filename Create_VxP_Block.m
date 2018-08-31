function [ UxP ] = Create_VxP_Block( mesh )

%         1   0   0   0 
%         0   1   0   0 
% B    =  0   0   1   0   nx 
%         0   0   0   1 
%               nx
%       
%
%         B  0  0  0  0
%        -B  B  0  0  0
% UxP  =  0 -B  B  0  0   ny +1
%         0  0 -B  B  0
%         0  0  0 -B  B
%         0  0  0  0 -B
%               ny

B = 1/mesh.h*speye(mesh.calculated.nx);
UxP =  kron(spdiags([-ones(mesh.calculated.ny+1,1),ones(mesh.calculated.ny+1,1)],[-1 0],...
                    mesh.calculated.ny+1,mesh.calculated.ny),B);
end