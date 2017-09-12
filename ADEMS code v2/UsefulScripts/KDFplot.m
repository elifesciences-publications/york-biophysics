function [Dens, x]=KDFplot(data)
[Dens,x] = ksdensity(data,'npoints',100);
plot(x,Dens)
end