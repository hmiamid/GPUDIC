function [u,v,M]=int_disp(I0,I1,hws,hss_col,hss_row,max_mem_k)
filter=ones(2*hws+1,1,'single','gpuArray');
N=gpuArray(single((2*hws+1)^2));

I0m=conv2(filter,filter,I0,'valid');
I0s=conv2(filter,filter,I0.^2,'valid');
I0s=sqrt(abs(N*I0s-I0m.^2));
I1m=conv2(filter,filter,I1,'valid');
I1s=conv2(filter,filter,I1.^2,'valid');
I1s=sqrt(abs(N*I1s-I1m.^2));

ndel=0;
k=0;
kmax=(2*hss_row+1)*(2*hss_col+1);
max_mem_k=min(max_mem_k,kmax);
mem_k=0;
previous_k=0;

M=zeros(size(I0m),'single','gpuArray')-1;
uv=ones(size(I0m),'single','gpuArray');
CCrow=zeros([size(I0m) max_mem_k],'single','gpuArray');

I1=padarray(I1,[hss_row,hss_col],NaN);
I1m=padarray(I1m,[hss_row,hss_col],NaN);
I1s=padarray(I1s,[hss_row,hss_col],NaN);

for row=-hss_row:hss_row
    
    I1_sub_row=I1(hss_row+row+1:end-hss_row+row,:);
    I1m_sub_row=I1m(hss_row+row+1:end-hss_row+row,:);
    I1s_sub_row=I1s(hss_row+row+1:end-hss_row+row,:);
    
    for col=-hss_col:hss_col

        I1_sub=I1_sub_row(:,hss_col+col+1:end-hss_col+col);
        I1m_sub=I1m_sub_row(:,hss_col+col+1:end-hss_col+col);
        I1s_sub=I1s_sub_row(:,hss_col+col+1:end-hss_col+col);
        
        I0I1=I0.*I1_sub;
        
        CCsub=conv2(filter,filter,I0I1,'valid');
        CCsub=N*CCsub-I0m.*I1m_sub;
        CCsub=CCsub./I0s./I1s_sub;
        
        CCrow(:,:,mem_k+1)=CCsub;
        mem_k=mem_k+1;
        k=k+1;
        
        string=sprintf('Integer search progress = %0.1f pc.\n',k/kmax*100);
        fprintf([repmat(8,[1 ndel]) string])
        ndel=numel(string);
            
        
        if mem_k==max_mem_k || k==kmax
            [max_CCrow,ind_max]=max(CCrow(:,:,1:mem_k),[],3);
            M=max(M,max_CCrow);
            m=M==max_CCrow;
            uv(m)=ind_max(m)+previous_k;
            mem_k=0;
            previous_k=k;
            
%             imagesc(M)
%             drawnow
        end
    end
end

uv=gather(uv);
v=ceil(uv/(2*hss_col+1))-hss_row-1;
u=1+mod(uv-1,2*hss_col+1)-hss_col-1;
M=gather(int8(M*128));