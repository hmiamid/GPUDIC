function Ib=conv2block_arrayfun2(I,xc,yc,hws)
    I=gpuArray(I);
    Ib=arrayfun(@bloc,xc,yc,reshape(gpuArray.colon(-hws,hws),[1 2*hws+1]),reshape(gpuArray.colon(-hws,hws),[1 1 2*hws+1]));

    function ib=bloc(x,y,i,j)
        ib=I(y+i,x+j);
    end
end