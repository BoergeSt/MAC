
X = mesh.calculated.X_u;
Y = mesh.calculated.Y_u;
normalx=zeros(size(X));
normaly = normalx;


for i = 1:numel(X)
    
    normal = Get_Outer_Normal( domain.calculated.on_domain, X(i),Y(i),mesh.h*1e-1,500);
    normalx(i)=normal(1);
    normaly(i)=normal(2);
end

%normaly(abs(abs(normaly)-1)<1e-5)=0;


figure()
quiver(X,Y,normalx,normaly,1,'b');