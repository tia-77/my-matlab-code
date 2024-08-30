clc;
clear;
close;
epoch=100; 
addpath('Functions')
run_no=1:epoch;
x = 5;
y = 8;
now2s1=zeros(epoch,900);
now2s2=zeros(epoch,900);
now2s3=zeros(epoch,900);
for ik=1:epoch
    HL = 1;
    P_k_store=DualControl(HL);
    step_no=size(P_k_store,1)-1;
    last_pos=P_k_store(1:end-1,1:2);
    now_pos=P_k_store(2:end,1:2);
    distance=sum(sqrt(sum((now_pos-last_pos).^2,2)));
    step_no1(ik)=step_no;
    distance1(ik)=distance;
    now2s1(ik,1:length(P_k_store))=sqrt(sum((P_k_store(:,1:2)-[x,y]).^2,2));
    HL = 2;
    P_k_store=DualControl(HL);
    step_no=size(P_k_store,1)-1;
    last_pos=P_k_store(1:end-1,1:2);
    now_pos=P_k_store(2:end,1:2);
    distance=sum(sqrt(sum((now_pos-last_pos).^2,2)));
    step_no2(ik)=step_no;
    distance2(ik)=distance;
    now2s2(ik,1:length(P_k_store))=sqrt(sum((P_k_store(:,1:2)-[x,y]).^2,2));
    HL = 3;
    P_k_store=DualControl(HL);
    step_no=size(P_k_store,1)-1;
    last_pos=P_k_store(1:end-1,1:2);
    now_pos=P_k_store(2:end,1:2);
    distance=sum(sqrt(sum((now_pos-last_pos).^2,2)));
    step_no3(ik)=step_no;
    distance3(ik)=distance;
    now2s3(ik,1:length(P_k_store))=sqrt(sum((P_k_store(:,1:2)-[x,y]).^2,2));
end
figure
plot(run_no,step_no1,'r',run_no,step_no2,'g',run_no,step_no3,'b','LineWidth',1);
legend('HL=1','HL=2','HL=3')
xlabel('No. of run');
ylabel('No. of step');
title('Number of Step Comparison For Different Horizon Length');

figure
plot(run_no,distance1,'r',run_no,distance2,'g',run_no,distance3,'b','LineWidth',1);
legend('HL=1','HL=2','HL=3')
xlabel('No. of run');
ylabel('distance');
title('Navigation Distance Comparison For Different Horizon Length');

step_flag1=logical(now2s1);
step_flag2=logical(now2s2);
step_flag3=logical(now2s3);
for i=1:epoch
    step_flag1(i,sum(step_flag1(i,:))+1)=1;
    step_flag2(i,sum(step_flag2(i,:))+1)=1;
    step_flag3(i,sum(step_flag3(i,:))+1)=1;
end

step_sum1=sum(step_flag1);step_sum1(step_sum1==0)=[];
step_sum2=sum(step_flag2);step_sum2(step_sum2==0)=[];
step_sum3=sum(step_flag3);step_sum3(step_sum3==0)=[];
now2s_average1=sum(now2s1(:,1:size(step_sum1,2)))./step_sum1;
now2s_average2=sum(now2s2(:,1:size(step_sum2,2)))./step_sum2;
now2s_average3=sum(now2s3(:,1:size(step_sum3,2)))./step_sum3;

figure
plot(1:size(step_sum1,2),now2s_average1,'r',1:size(step_sum2,2),now2s_average2,'g',1:size(step_sum3,2),now2s_average3,'b','LineWidth',1);
legend('HL=1','HL=2','HL=3')
xlabel('steps');
ylabel('distance to source');
title('Average Step Distance to Source Comparison For Different Horizon Length');
