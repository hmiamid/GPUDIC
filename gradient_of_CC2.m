function gC=gradient_of_CC2(d,dI0dx,dI0dy,epsfiltx,epsfilty)

dx=d.*dI0dx;
dy=d.*dI0dy;

gC(:,1)=squeeze(sum(dx,[2 3]));
gC(:,2)=squeeze(sum(dy,[2 3]));
gC(:,3)=squeeze(sum(epsfiltx.*sum(dx,2),3));
gC(:,4)=squeeze(sum(epsfilty.*sum(dx,3),2));
gC(:,5)=squeeze(sum(epsfiltx.*sum(dy,2),3));
gC(:,6)=squeeze(sum(epsfilty.*sum(dy,3),2));