function Hinv=calc_hessian(dI0dx,dI0dy,I0std,hws,xc,yc,queuesize)

elements=numel(xc);
queuesize=min(queuesize,elements);
sectlims=unique([0:queuesize:elements elements]);

ndel=0;

for sect=2:numel(sectlims)
    r=(sectlims(sect-1)+1:sectlims(sect))';
    dI0dx_b=conv2block_arrayfun2(dI0dx,xc(r),yc(r),hws);
    dI0dy_b=conv2block_arrayfun2(dI0dy,xc(r),yc(r),hws);

    % x2=dI0dx_b.^2;
    % y2=dI0dy_b.^2;
    % xy=dI0dx_b.*dI0dy_b;

    epsfiltx=reshape(gpuArray.colon(-hws,hws),[1 1 2*hws+1]);
    epsfilty=reshape(gpuArray.colon(-hws,hws),[1 2*hws+1 1]);
    H=zeros([length(r) 6 6],'single','gpuArray');
    for i=1:6
        for j=1:6
            if any(i==[1 3 4])
                d=dI0dx_b;
            else
                d=dI0dy_b;
            end
            if any(j==[1 3 4])
                d=d.*dI0dx_b;
            else
                d=d.*dI0dy_b;
            end
            if any(i==[3 5])
                d=d.*epsfiltx;
            elseif any(i==[4 6])
                d=d.*epsfilty;
            end
            if any(j==[3 5])
                d=d.*epsfiltx;
            elseif any(j==[4 6])
                d=d.*epsfilty;
            end
            H(:,i,j)=sum(d,[2 3]);
        end
    end
    c=yc(r)-hws+(xc(r)-hws-1)*size(I0std,1);
    H=H./I0std(c).^2;
    H=permute(pagefun(@inv,permute(H,[2 3 1])),[3 1 2]);
    Hinv(r,:,:)=gather(H);
    
    progress=sect/numel(sectlims);
    string=sprintf('Calculate hessians = %0.1f pc.\n',progress*100);
    fprintf([repmat(8,[1 ndel]) string])
    ndel=numel(string);
    
end