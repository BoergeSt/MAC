function [ UxU ] = Create_UxU_Block( mesh )
%% Assembly of the first equality in Matrixform
% 1/h^2 * (4u_(i,j)-u_(i-1,j)-u_(i+1,j)-u_(i,j-1)-u_(i,j+1)) + 1/h * ...
%   ... (p_(i,j)-p_(i-1,j)) = f_1(x_(i,j))

% This leads to the following matrices (no-size correspondence)
%         4  -1   0   0   0
%        -1   4  -1   0   0
% A_11 =  0  -1   4  -1   0  
%         0   0  -1   4  -1
%         0   0   0  -1   4
%       
%
%         A_11   -Id     0     0     0
%          -Id  A_11   -Id     0     0
% A_1  =     0   -Id  A_11   -Id     0
%            0     0   -Id  A_11   -Id
%            0     0     0   -Id  A_11

e = ones(mesh.calculated.nx+1,1);
A_11 = spdiags([-e, 4.*e, -e], -1:1 ,mesh.calculated.nx+1,mesh.calculated.nx+1);
A_p1 = kron(speye(mesh.calculated.ny,mesh.calculated.ny),A_11);
A_12 = kron(spdiags([-ones(mesh.calculated.ny,1),-ones(mesh.calculated.ny,1)],...
            [-1 1],mesh.calculated.ny,mesh.calculated.ny),...
            speye(mesh.calculated.nx+1,mesh.calculated.nx+1));
UxU =  1/mesh.h^2 * (A_p1 + A_12);
