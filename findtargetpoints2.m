function [x2,y2]=findtargetpoints2(xc,yc,Mp,hws)

epsfilty=reshape(gpuArray.colon(-hws,hws),[1 2*hws+1 1]);
epsfiltx=reshape(gpuArray.colon(-hws,hws),[1 1 2*hws+1]);

x2=xc+Mp(:,1,3);
x2=x2+Mp(:,1,2).*epsfilty;
x2=x2+Mp(:,1,1).*epsfiltx;


y2=yc+Mp(:,2,3);
y2=y2+Mp(:,2,2).*epsfilty;
y2=y2+Mp(:,2,1).*epsfiltx;