function result=mixture_Normal_fitting_covariance_matrix_complex2(data,k,nit,nmc,nthin)
dim = 2*size(data,2);
n = size(data,1);           %Number of observations
%Rearrange the data
data = [real(data.'); imag(data.')]'; 
%prior of any mu_k is MVN(mu0,sigma2)
mu0 = zeros(dim,1);
% A = rand(dim,dim)*100+100;
% sigma0 = 0.5*(A+A');
sigma0 = eye(dim);
%prior for any sigma_k is InvWish_k(Psi,df0)
Psi = eye(dim);
df0 = dim + 1;
%prior for (p1,p2,...,pk) is Dirichlet (alpha_1,alpha_2,...,alpha_k)
alpha0=repmat(3,k,1);
p_prior=zeros(k,1);
%informative starting point - kmeans
cluster_index=kmeans(data,k);
%now based on what the cluster index is, compute initial values of mu_k, sigma_k, p_k
for i=1:k
    mu(i,:) = mean(data(find(cluster_index==i),:));
    sigma(:,:,i) = cov(data(find(cluster_index==i),:));
    p(i) = length(find(cluster_index==i))/n;
end
%define auxilliary variable z, one for each observation
z=cluster_index;
%rearrange for identifiability according to increasing value of mu
% c=tiedrank(vecnorm(mu,2,2));
% mu=mu(c,:);
% sigma=sigma(:,:,c);
% p=p(c);
% z_temp=z;
% for i=1:k
%     z_temp(find(z==i))=c(i);
% end
% z=z_temp;
%during the mcmc, you need to store the simulated samples of mu, sigma2, p and z, so define storage variables
A=[(1:n)',z];
prob=zeros(n,k);
mu_store=zeros(nmc/nthin,dim,k);
sigma_store=zeros(nmc/nthin,dim,k);
p_store=zeros(nmc/nthin,k);
z_store=zeros(n,k);         %every time simulated a sequence to specify from which distribution
%start the mcmc loop
for iter=1:(nit+nmc)
    for i=1:k
        pos = find(z==i);
        l_pos = length(pos);
        %posterior distribution of sigma is inverse wishart
        s1 = 0;
        s1 = (data(pos,:)-mu(i,:)).';
        s1 = reshape(s1,dim,1,size(s1,2));
        s1T = permute(s1,[2 1 3]);
        S = zeros(dim);
        for l = 1:size(s1,3)
            S = S + s1(:,:,l)*s1T(:,:,l);
        end
        sigma_post_scale = S + Psi;
        df = l_pos + df0;
        %then simulate sample of sigma2
        sigma(:,:,i) = iwishrnd(sigma_post_scale,df); %*(df-dim-1);
        %posterior distribution of mu is MVN
        mu_post_sigma = inv ( l_pos*inv(sigma(:,:,i)) + inv(sigma0) );
        mu_post_mean = ( ( l_pos .* mean(data(pos,:),1) / sigma(:,:,i) ) + ( mu0.' / sigma0 ) ) * mu_post_sigma;
        %then simulate sample of mu
        mu(i,:) = mvnrnd(mu_post_mean,mu_post_sigma);
        p_prior(i) = alpha0(i) + l_pos;
    end
    %Simultaneously sample the p-vector
    p=dirichletRnd(p_prior);
    %Simulate the samples of z-vector
    prob(:,1)=p(1)*mvnpdf(data,mu(1,:),sigma(:,:,1));
    for j=2:k 
        prob(:,j)=prob(:,j-1)+p(j)*mvnpdf(data,mu(j,:),sigma(:,:,j));
    end
    %z=indicator_simulation3(prob./repmat(prob(:,k),1,size(prob,2)));
    z=indicator_simulation4(prob./repmat(prob(:,k),1,size(prob,2)));
    %rearrange for identifiability according to increasing value of mu
    c=tiedrank(vecnorm(mu,2,2));
    mu(c,:);
%     sigma=sigma(:,:,c);
%     p=p(c);
%     z_temp=z;
%     for i=1:k
%         z_temp(find(z==i))=c(i);
%     end
%     z=z_temp;
    if(iter>nit && mod((iter-nit),nthin)==0)
        A(:,2)=z;
        for i = 1:k
            mu_store((iter-nit)/nthin,:,i)=mu(i,:);
            sigma_store((iter-nit)/nthin,:,i)=diag(sigma(:,:,i));
        end
        p_store((iter-nit)/nthin,:)=p;
        for i=1:size(A,1)
            z_store(A(i,1),A(i,2))=z_store(A(i,1),A(i,2))+1;
        end
%         for i=1:k
%             density_store=density_store+p(i)*normpdf(grid_values,mu(i),sigma2(i));
%         end
    end
end
result=struct;
result.a=mu_store;
result.b=sigma_store;
result.c=p_store;
result.d=z_store;
% result.e=grid_values;
% result.f=density_store/(nmc/nthin);
end