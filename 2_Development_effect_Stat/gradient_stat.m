%%examine the age-associated changes

% calculate the global gradient topography
indi_grad1(:,:)=realigned(:,1,:);
indi_grad2(:,:)=realigned(:,2,:);
g1_thres(:,1)=max(indi_grad1);
g1_thres(:,2)=min(indi_grad1);
range_indi_grad(:,1)=range(indi_grad1,2);
range_indi_grad(:,2)=range(indi_grad2,2);
std_indi_grad(:,1)=std(indi_grad1,0,2);
std_indi_grad(:,2)=std(indi_grad2,0,2);
save('indi_grad1','indi_grad1');
save('indi_grad2','indi_grad2');
save('range_indi_grad','range_indi_grad');
save('std_indi_grad','std_indi_grad');

gradient_order_resort=zeros(39,16); % 39: number of subjects
for i=1:39
    A=realigned_gradient.xfms{i};
    gradient_order=abs(A);
    [max_gradient_order,index]= max(gradient_order);
    gradient_order_resort(i,:)=index; 
end

gradient_indi_explan=zeros(39,77);
for i=1:39
load(['H:\Gradient_31-42\Results\third_analysis\indi_gradient\g',subname(i).name]);
explan=gradient.res.lambdas./sum(gradient.res.lambdas);
gradient_indi_explan(i,:)=explan;
end

gradient_indi_explan_resort=zeros(39,2); 
for i=1:39
gradient_indi_explan_resort(i,1)=gradient_indi_explan(i,gradient_order_resort(i,1));
gradient_indi_explan_resort(i,2)=gradient_indi_explan(i,gradient_order_resort(i,2));
end

save('gradient_order_resort','gradient_order_resort');
save('gradient_indi_explan','gradient_indi_explan');
save('gradient_indi_explan_resort','gradient_indi_explan_resort');

%sliding windows  
realigned_gradient_sort_age=zeros(6096,16,39);
age=table_model_resort.age;
for i=1:39
    index=order_sort_age(i);
    realigned_gradient_sort_age(:,:,i)=realigned_gradient.realigned(:,:,index);
    realigned_age(i,:)=age(index,:);
    realigned_explan(i,:)=gradient_indi_explan_resort(index,:);
end
save('realigned_age','realigned_age')
save('realigned_explan','realigned_explan')
save('realigned_gradient_sort_age','realigned_gradient_sort_age');

win_length=10;
win_step=2;
for i= 1:15
  gradient_window = realigned_gradient_sort_age(:,1:2,1+(i-1)*win_step:win_length+(i-1)*win_step);
  save(['gradient_window',mat2str(i)],'gradient_window');
  gradient1_window(:,i)=mean(gradient_window(:,1,:),3);
  gradient2_window(:,i)=mean(gradient_window(:,2,:),3);
  age_window = realigned_age(1+(i-1)*win_step:win_length+(i-1)*win_step,:);
  mean_age_window(i,:)=mean(age_window);
  explan_window = realigned_explan(1+(i-1)*win_step:win_length+(i-1)*win_step,:);
  mean_explan_window(i,:)=mean(explan_window);
  
end
save('gradient1_window','gradient1_window');
save('gradient2_window','gradient2_window');
save('mean_age_window','mean_age_window');
save('mean_explan_window','mean_explan_window');

% calculate maturation index
load('gradient_win15.mat');
g1=gradient.emb(:,1);
g2=gradient.emb(:,2);
corrcoef(indi_grad2',g2);
for i=1:39
[r,p]=corrcoef(indi_grad2(i,:)',g2);
matura_g2(i)=r(1,2);
end

%general linear model
X(:,1) = table_model_resort.age;
X(:,2) = table_model_resort.sex;
X(:,3) = table_model_resort.meanFD;
for i = 1:6
    Y=gradient_meas(:,i);
    [bb,dev,stats] = glmfit(X,Y);
    results{i}.stats = stats;
    results{i}.stats = stats;
    age_related_global_grad(i,1) = stats.t(2,1);
    age_related_global_grad(i,2) = stats.p(2,1);
end
save age_related_global_grad_t age_related_global_grad
save age_related_global_grad_results results

for i = 1:6096
    Y=indi_grad1(:,i);
    [bb,dev,stats] = glmfit(X,Y); 
    results{i}.stats = stats;
    results{i}.stats = stats;
    age_related_regional_grad(i,1) = stats.t(2,1);
    age_related_regional_grad(i,2) = stats.p(2,1);
end
save age_related_regional_grad2_t age_related_regional_grad
save age_related_regional_grad2_results results








