clear all
clc


fname0='summary_data.xlsx';

%%%%%%%%%%%%%%%%%%%%%%%% load data


control_data = xlsread(fname0,'Sheet1','D3:G12');
cms_data = xlsread(fname0,'Sheet1','D19:G28');

figure
subplot(1,3,1)
q1=plot(1:4,nanmean(control_data),'b')
hold on
errorbar(1:4,nanmean(control_data),nanstd(control_data)./sqrt(10-sum(isnan(control_data),1)),'b')
title('Sham and cMS experiments')
q2=plot(1:4,nanmean(cms_data),'r')
errorbar(1:4,nanmean(cms_data),nanstd(cms_data)./sqrt(10-sum(isnan(cms_data),1)),'r')
lgd=legend([q1 q2],'sham','cMS')
lgd.FontSize=8
ylim([0 0.8])
xlim([0 5])

clear B Res
c1=0;
clear X Y

subplot(1,3,2)

xs=1:0.1:4;
for i=1:10
    y=cms_data(i,:);
    p = polyfit(1:1:4,y,3);
    Coffs_CMS(i,:)=p;
end

B=nanmean(Coffs_CMS);
txt2 = strcat('Average cMS fit: ',num2str(B(1)),'x^3+',num2str(B(2)),' x^2',num2str(B(3)),'x+',num2str(B(4)));
zs2=B(4)+B(3)*xs+B(2)*xs.*xs+B(1)*xs.*xs.*xs;
q1=plot(xs,zs2,'r','Linewidth',2)
hold on
C=nchoosek(1:10,5);

for w=1:size(C,1)
    F=randperm(10);
    clear Coffs_Control
    for ii=1:5
        i=squeeze(C(w,ii));
        y=control_data(i,:);
        p = polyfit(1:1:4,y,3);
        Coffs_Control(ii,:)=p;
    end
     
    A=nanmean(Coffs_Control);
    
    zs1=A(4)+A(3)*xs+A(2)*xs.*xs+A(1)*xs.*xs.*xs;
    allzs(w,:)=zs1;
    if w==1
        q2=plot(xs,zs1,'b','Linewidth',0.1)
    else
        plot(xs,zs1,'b','Linewidth',0.1)
    end

end

plot(xs,zs2,'r','Linewidth',2)
hold on
lgd=legend([q1 q2],'cms fit','sham fits')
lgd.FontSize=8
t2=text(0.3, 0.1,txt2,'FontSize',8)
t2.Color='r'
t2.FontSize=8
ylim([0 0.8])
xlim([0 5])

subplot(1,3,3)
p1=plot(xs,zs2,'r','Linewidth',2)
hold on
p2=plot(xs,max(allzs),'b')
plot(xs,min(allzs),'b')

for i=1:length(xs)
    [C,ix]=sort(allzs(:,i),'descend');
    upper(i)=C(length(ix)-12);
    [C,ix]=sort(allzs(:,i),'ascend');
    lower(i)=C(12);
    [C,ix]=sort(allzs(:,i),'ascend');
    med(i)=C(126);
end
p3=plot(xs,upper,'b--')
plot(xs,lower,'b--')
p4=plot(xs,med,'b','Linewidth',2)
ylim([0 0.8])
xlim([0 5])

lgd=legend([p1 p2 p3 p4],'cms fit','sham max and min fits','5 and 95 percentile','sham median')
lgd.FontSize=8

