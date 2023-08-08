%%  prediction based on SVR and LOOCV
clc;clear all;close all;

X=indi_grad1(index_beha,:);
Y=behaviors(:,1);

%% K-fold cross validation
tic
data=X; % n*m matrix, n is the total sub number, m is the feature dimension
label=Y; %prediction label
SubNum=size(data,1);
K=5;  % set K = n for leave-one-out-cross-validation / or K=k for k-fold cross validation
indices=crossvalind('Kfold',SubNum,K);
predict_label_all=[];
testlabel_all=[];

for v=1:K
    test=[];
    trainall=[];
    test=(indices==v);
    trainall=~test;
    traindata=[];
    testdata=[];
    trainlabel=[];
    testlabel=[];
    predict_label=[];
    
    traindata= [traindata data(trainall,:)];
    testdata = [testdata data(test,:)];
    
    trainlabel = label(trainall,:);
    testlabel = label(test,:);
    
    
    %% SVR
    % model = fitcsvm(traindata,trainlabel);% -t 2 RBF   -t 0 Linear
    % [predict_label,mse,prob_estimates] = svmpredict(testlabel,testdata, model);

    model_svr = fitrsvm(traindata,trainlabel,'KernelFunction','linear');
    predict_label = predict(model_svr,testdata);

    predict_label_all=[predict_label_all; predict_label];
    testlabel_all=[testlabel_all; testlabel];
    
end

[r,p]=corr(predict_label_all,testlabel_all);

toc

real_r=r;
Yrand=zeros(26,1000);
for i=1:1000
index=randperm(26);
Yrand(:,i)=testlabel_all(index,1);
end
parfor i = 1:1000
  [r,p]= corr(predict_label_all,Yrand(:,i));
   surr_r(:,i) = r;
end
p_real=sum(gt(surr_r,real_r),2)/1000;

save predict_cognitive_g1_5_fold_result2 testlabel_all predict_label_all p_real real_r p model_svr


grad=zeros(1,length(mask_ind)); 
grad(mask_ind)=model_svr.Beta; % indi_gradient
[dim1,dim2,dim3]=size(mask_vol);
grad_nii=reshape(grad,dim1,dim2,dim3);
mask_hdr.fname='...\predict_cognitive_g1_5fold_comp_beta.nii';
mask_hdr.dt=[16,0]; 
spm_write_vol(mask_hdr,grad_nii);
