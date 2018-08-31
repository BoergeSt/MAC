function [ s ] = Find_Border_Distance(testfunc,x,y,xdir,ydir,iteration_number)

if testfunc(x+xdir,y+ydir)
    s = 1;
    return
end
u = 1;
l = 0;
for i = 1:iteration_number
    s = (u+l)/2;
    if testfunc(x+s*xdir,y+s*ydir)
        l = s;
    else
        u = s;
    end
end
s = (u+l)/2;

