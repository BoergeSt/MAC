function [ normal ] = Get_Outer_Normal( on_domain, x,y, h,n)

normal = [0;0];

theta = linspace(0,2*pi,n+1);
theta=theta(1:n);
[dx,dy] = pol2cart(theta,ones(1,n));

on_domain_vector = on_domain(x+dx*h,y+dy*h);

if all(on_domain_vector) || all(~on_domain_vector) 
    return
end
    
%change = diff([on_domain_vector(end),on_domain_vector]);
change = diff([on_domain_vector,on_domain_vector(1)]);

if nnz(change)>2
    warning('normal calculation not possible, h not small enough to dissolve the domain');
    return
end

index = round(mean(find(change)));

normal(1) = dx(index);
normal(2) = dy(index);

if on_domain_vector(index)
    normal = -normal;
end
%if on_domain(x+normal(1),y+normal(2))
%    normal = -normal;
%end

normal(abs(normal)<100*eps)=0;


end
