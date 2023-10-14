function [Mpout,CC,counterout]=IC_method(I0_info,I1,Mp,hws,xc,yc,conv_iter,tolerance,queuesize)
N=(2*hws+1)^2;
epsfilty=single(gpuArray.colon(-hws,hws));
epsfiltx=permute(epsfilty,[3 1 2]);
fprintf('5C\n')
% F=quintic_coeffs(I1);
[x,y]=meshgrid(gpuArray(single(1:size(I1,2))),gpuArray(single(1:size(I1,1))));
fprintf('5C done\n')

CC=NaN(size(xc),'single','gpuArray');

Mptemp=reshape(Mp,[size(Mp,1)*size(Mp,2) 2 3]);
Mpout=zeros(size(Mptemp),'single');
elements=numel(xc);
index=yc-hws+(xc-hws-1)*size(I0_info.I0m,1);

queuesize=min(queuesize,elements);
inqueue=(1:queuesize)';
outqueue=((queuesize+1):elements)';

% Load queue to GPU
Mp_b=gpuArray(Mptemp(inqueue,:,:));
H_b=gpuArray(I0_info.H(inqueue,:,:));
I0n=conv2block_arrayfun2(I0_info.I0,xc(inqueue),yc(inqueue),hws);
Is=I0_info.I0std(index(inqueue));
I0n=(I0n-I0_info.I0m(index(inqueue)))./Is;
dI0dx_b=conv2block_arrayfun2(I0_info.dI0dx,xc(inqueue),yc(inqueue),hws)./Is;
dI0dy_b=conv2block_arrayfun2(I0_info.dI0dy,xc(inqueue),yc(inqueue),hws)./Is;
dOld=Inf([queuesize 6],'single','gpuArray');
counter=zeros([queuesize 1]);
counterout=zeros(size(xc));

ndel=0;
progress=0;

while ~isempty(inqueue)
    [x2,y2]=findtargetpoints2(xc(inqueue),yc(inqueue),Mp_b,hws);

%     I1n=interp_quintic_arrayfun(F,x2,y2);
    I1n=interp2(x,y,I1,x2,y2,'cubic');
    
    I1im=sum(I1n,[2 3])/N;
    I1is=single(sqrt(abs(sum(double(I1n).^2,[2 3])-N*double(I1im).^2)));
    I1n=I1n-I1im;
    I1n=I1n./I1is;
    CC(inqueue)=sum(I1n.*I0n,[2,3]);
%     imagesc(CC)
%     caxis([-1 1])
%     drawnow
    gC=gradient_of_CC2(I1n-I0n,dI0dx_b,dI0dy_b,epsfiltx,epsfilty);
    d=gradient_descent(H_b,gC);
    Mp_b=inverseComp2(Mp_b,d);
    
    delta=(dOld(:,1)-d(:,1)).^2+(dOld(:,4)-d(:,4)).^2;
    dOld=d;
    counter=counter+1;
    finished=find(gather(delta<tolerance) | counter>=conv_iter);
    sumfinished=numel(finished);
    
    if sumfinished~=0
        Mpout(inqueue(finished),:,:)=gather(Mp_b(finished,:,:));
        counterout(inqueue(finished))=counter(finished);
        
        progress=progress+sumfinished;
        string=sprintf('Progress = %0.1f pc.\n',progress/elements*100);
        fprintf([repmat(8,[1 ndel]) string])
        ndel=numel(string);
        
        if length(outqueue)>=sumfinished
            newelements=outqueue(1:sumfinished);
            outqueue(1:sumfinished)=[];
            toreplace=finished;
            todelete=[];
        else
            len=length(outqueue);
            newelements=outqueue;
            outqueue=[];
            toreplace=finished(1:len);
            todelete=finished((len+1):end);
        end
        
        if ~isempty(toreplace)
            inqueue(toreplace)=newelements;
            Mp_b(toreplace,:,:)=gpuArray(Mptemp(newelements,:,:));
            H_b(toreplace,:,:)=gpuArray(I0_info.H(newelements,:,:));
            I0n(toreplace,:,:)=conv2block_arrayfun2(I0_info.I0,xc(newelements),yc(newelements),hws);
            Is=I0_info.I0std(index(newelements));
            I0n(toreplace,:,:)=(I0n(toreplace,:,:)-I0_info.I0m(index(newelements)))./Is;
            dI0dx_b(toreplace,:,:)=conv2block_arrayfun2(I0_info.dI0dx,xc(newelements),yc(newelements),hws)./Is;
            dI0dy_b(toreplace,:,:)=conv2block_arrayfun2(I0_info.dI0dy,xc(newelements),yc(newelements),hws)./Is;
            dOld(toreplace)=Inf;
            counter(toreplace)=0;
        end
        
        if ~isempty(todelete)
            inqueue(todelete)=[];
            Mp_b(todelete,:,:)=[];
            H_b(todelete,:,:)=[];
            I0n(todelete,:,:)=[];
            dI0dx_b(todelete,:,:)=[];
            dI0dy_b(todelete,:,:)=[];
            dOld(todelete,:)=[];
            counter(todelete)=[];
        end
    end
end
Mpout=reshape(Mpout,size(Mp));