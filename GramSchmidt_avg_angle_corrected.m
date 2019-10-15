function e = GramSchmidt_avg_angle(v,dim)
%minimum coefficient
[n,m] = size(v);
% p = randperm(n);
% v = v(p,:);
avg = mean(v);
v = [avg;v];
for i = 1:dim
    u(i,:) = v(i,:); 
    for j = 1:size(u,1)-1
        u(i,:) = u(i,:) - proj( u(j,:) , v(i,:) );
    end
    e(i,:) = u(i,:) / norm (u(i,:));
    for j=1:n-i+1 
        pro(j,:) = proj_basis( e, v(i+j,:) );
    end
    cosine = abs(sum(pro.*v(i+1:n+1,:),2)./(vecnorm(pro,2,2).*vecnorm(v(i+1:n+1,:),2,2)));
    cosine(isnan(cosine)) = 0;
    [M,k] = min(cosine);
    temp = v(k+i,:);
    v(k+i,:) = v(i+1,:);
    v(i+1,:) = temp;
    pro = [];
end
end