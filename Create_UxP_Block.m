function [ UxP ] = Create_UxP_Block( mesh )

%         1   0   0   0 
%        -1   1   0   0 
% B    =  0  -1   1   0   nx+1 
%         0   0  -1   1 
%         0   0   0  -1 
%               nx
%       
%
%         B  0  0  0  0
%         0  B  0  0  0
% UxP  =  0  0  B  0  0   ny
%         0  0  0  B  0
%         0  0  0  0  B
%               ny

e = 1/mesh.h*ones(mesh.calculated.nx+1,1);
B = spdiags([-e, e], -1:0 ,mesh.calculated.nx+1,mesh.calculated.nx);
UxP =  kron(speye(mesh.calculated.ny,mesh.calculated.ny),B);
end

