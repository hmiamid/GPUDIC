function I0_info=calculate_I0_parameters(I0,hws,xc,yc,queuesize)
[dI0dx,dI0dy]=gradofI(I0);
I0_info.dI0dx=dI0dx;
I0_info.dI0dy=dI0dy;

filt=ones(2*hws+1,1,'double','gpuArray');N=length(filt)^2;

I0m=conv2(filt,filt,I0,'valid')/N;
I0std=single(sqrt(abs(conv2(filt,filt,double(I0).^2,'valid')-N*double(I0m).^2)));
I0_info.I0std=I0std;I0_info.I0m=I0m;

I0_info.H=calc_hessian(I0_info.dI0dx,I0_info.dI0dy,I0_info.I0std,hws,xc,yc,queuesize);
I0_info.I0=I0;