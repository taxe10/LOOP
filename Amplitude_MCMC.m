clear
close all
clc
%Load data
load('Mouse5B_Fresh.mat')   %Change to Mouse5B_Block.mat for block data
type_data = 'MOUSE';
type_tissue = 'FRESH';      %Change to BLOCK for block data
%Prepare the data
Fs = 1 / (trange(2)-trange(1)); % Sampling frequency       
L = length(trange);             % Length of signal
f = Fs*(0:((L/2)-1))/L;
[c,n1]=min(abs(f-0.1E12));      % Calculate the window range 0.1-4THz
[c,n2]=min(abs(f-4E12));
%Frequency domain
raw_fft = fft(ScanData,L,3);
response = abs(raw_fft(:,:,1:L/2)/L); %Normalized
data_case='5B';
distr ='NORMAL';
k=2;                            %Number of regions
data_red = 'LOOP';
%Dimension
dim = 2;
%Define save path
filename0=char(strcat(num2str(k),{'_'},lower(distr),{'_components_for_'},lower(type_data),{'_'},data_case,{'_'},lower(type_tissue)));
[filename,path] = uiputfile(strcat(filename0,'.fig'),'SAVE MESHES AS:');

rng default     % For reproducibility
%Prepare data for mixture
mask_original = flipud(matrix);
mask = reshape(mask_original,[],1);
pos_nonzero = find (mask > 0);
pos0 = find(mask <= 0);
true_mat_small = mask(pos_nonzero);
data = reshape(response,[size(ScanData,1)*size(ScanData,2) L/2]);
data = data(:,n1:n2);
response_time_summarized = response;
data(pos0,:) = [];      %Consider pixels inside the region of interest only
%Dimension reduction
e = GramSchmidt_avg_angle_corrected(data,dim)';
sum_data = e\data';
data = sum_data';
n_x = length(xrange);
n_y = length(yrange);
%Fit a mixture model with k components and import the output
nit=15000;          %Initial observations to discard
nmc=15000;          %Number of observations to consider
nthin=5;            %Interval to thin
nsample=nmc/nthin;  %Final number of samples retained
scale = 15;         %Scale definition
out=mixture_Normal_fitting_covariance_matrix_v2(scale*data,k,nit,nmc,nthin);
mc_samples_store_mean=out.a/scale;
mc_samples_store_var=out.b/(scale^2);
mc_samples_store_p=out.c;
indicator = out.d;
indicator_proportion=indicator/nsample;
[max_num, max_idx]=max(indicator_proportion'); %max per row
indicator_modal_class=ind2sub(size(indicator_proportion'),max_idx); %position instead of number
indicator_proportion_intermediate=repmat(NaN,n_x*n_y,k);
indicator_modal_class_intermediate=repmat(NaN,n_x*n_y,k);
indicator_proportion_intermediate(pos_nonzero,:)=indicator_proportion;
indicator_modal_class_intermediate(pos_nonzero,:)=repmat(indicator_modal_class',1,k);
indicator_restructured=zeros(n_y,n_x,k);
for h=1:k
    indicator_restructured(:,:,h)=reshape(indicator_proportion_intermediate(:,h),n_y,n_x);
end
indicator_modal_class_restructured = reshape(indicator_modal_class_intermediate(:,1),n_y,n_x); 
%ROC curves
h=figure(1);
for i=1:k
    subplot(1,k,i)
    prob_small=indicator_restructured(:,:,i);
    prob_small(pos0)=[];
    [X,Y,T,AUC]=perfcurve(double(true_mat_small==i),prob_small,1);
    plot(X,Y,'linewidth',2)
    hold on
    plot(0:1,0:1,'k','linewidth',2)
    grid on
    m=find(matrix_key_number==i);
    xlabel(['False Detection Rate ',matrix_key_tissue(m,:)])
    ylabel(['True Detection Rate ',matrix_key_tissue(m,:)])
    title(sprintf('ROC curve (AUC=%0.4f)',AUC))
end
set(gcf, 'Units', 'inches','Position', [0 0 5*k 5]);
filename1=strcat(path,strcat({'ROC_curve_'},filename));
savefig(h,char(filename1));
saveas(h,char(strcat(path,{'ROC_curve'},filename0,{'.eps'})),'epsc');
save(char(strcat(path,{'VARIABLES '},filename0)))
close(gcf)

x = xrange;
y = yrange;
h=figure(2);
for i=1:k
    subplot(k,2,2*i-1)
    true_mat_temp=int8(mask_original==i);
    pcolor(x,y,true_mat_temp) %plot for pathology
    colormap jet
    shading flat
    m=find(matrix_key_number==i);
    title(['Pathology Detection of ',matrix_key_tissue(m,:)])
    subplot(k,2,i*2)
    indicator_modal_class_restructured_new=int8(indicator_modal_class_restructured==i);
    pcolor(x,y,indicator_modal_class_restructured_new) %plot for model
    colormap jet
    shading flat
    title(['Model-based Detection of ',matrix_key_tissue(m,:)])
end
filename2=strcat(path,strcat({'Path_vs_Model_per_Region_'},filename));
savefig(h,char(filename2));
saveas(h,char(strcat(path,{'Path_vs_Model_per_Region_'},filename0,{'.eps'})),'epsc');
save(char(strcat(path,{'VARIABLES '},filename0)))
close(gcf)

h=figure(3);
for i=1:k
    subplot(1,k,i)
    prob_mat=indicator_restructured(:,:,i);
    prob_mat(pos0)=repmat(-0.05,size(pos0));
    pcolor(x,y,prob_mat)
    colorbar
    caxis([0,1])
    colormap jet
    shading flat
    m=find(matrix_key_number==i);
    title(['Probability of ',matrix_key_tissue(m,:)])
end
set(gcf, 'Units', 'inches','Position', [1 1 5*k 4]);
filename3=strcat(path,strcat({'Probability_of_detection_'},filename));
savefig(h,char(filename3))
saveas(h,char(strcat(path,{'Probability_of_detection_'},filename0,{'.eps'})),'epsc');
save(char(strcat(path,{'VARIABLES '},filename0)))
close(gcf)

h=figure(4);
subplot(1,2,1)
pcolor(x,y,indicator_modal_class_restructured)
caxis ([1-0.01,k+0.01])
colormap jet
shading flat
title('Model class')
subplot(1,2,2)
mask_to_plot=mask_original;
mask_to_plot(mask_original<=0)=NaN;
pcolor(x,y,abs(mask_to_plot))
caxis ([1-0.01,k+0.01])
colormap jet
shading flat
title('Pathology class')
set(gcf, 'Units', 'inches','Position', [1 1 10 4]);
filename4=strcat(path,strcat({'Model_VS_Path_'},filename));
savefig(h,char(filename4))
saveas(h,char(strcat(path,{'Model_VS_Path_'},filename0,{'.eps'})),'epsc');
saveas(h,char(strcat(path,{'Model_VS_Path_'},filename0,{'.png'})),'png')
save(char(strcat(path,{'VARIABLES '},filename0)))
close(gcf)

save(char(strcat(path,{'VARIABLES '},filename0)))
close(gcf)
