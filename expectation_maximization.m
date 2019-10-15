function out = expectation_maximization(data,k)
dim = size(data,2);         %Dimension
n = size(data,1);           %Number of observations
threshold = 0.0001;
%Initialize - kmeans
cluster_index = kmeans(data,k);
for i=1:k
    mu(i,:) = mean(data(find(cluster_index==i),:));
    sigma(:,:,i) = cov(data(find(cluster_index==i),:));
    p(i) = length(find(cluster_index==i))/n;
    den(i,:) = p(i) * mvnpdf(data,mu(i,:),sigma(:,:,i));
end
m = 1;
L(m) = sum(log(sum(den)))/n;     %Likelihood
z = cluster_index;            %Latent variable
while 1
    m = m+1;
    gamma = (den ./ sum(den))';
    nj = sum(gamma)';
    p = nj / n;
    mu = ((gamma'*data)'./ nj')';
    sigma = zeros(dim,dim,k);
    for j = 1:k
        for i = 1:n
            sigma(:,:,j) = sigma(:,:,j) + gamma(i,j)*(data(i,:)' - mu(j,:)')*(data(i,:)' - mu(j,:)')';
        end
        sigma(:,:,j) = sigma(:,:,j)/nj(j);
        den(j,:) = p(j) * mvnpdf(data,mu(j,:),sigma(:,:,j));
    end
    L(m) = sum(log(sum(den)))/n;     %Update likelihood
    if abs(L(m) - L(m-1)) <= threshold
        break;
    end
end
out = struct;
out.a = mu;
out.b = sigma;
out.c = p;
out.d = gamma;
end