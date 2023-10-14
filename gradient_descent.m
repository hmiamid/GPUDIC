function d=gradient_descent(Hinv,gC)
% d=zeros(size(gC),'single','gpuArray');
% for i=1:6
%     for j=1:6
%         d(:,i)=d(:,i)+Hinv(:,i,j).*gC(:,j);
%     end
% end
d=permute(pagefun(@mtimes,permute(Hinv,[2 3 1]),permute(gC,[2 3 1])),[3 1 2]);