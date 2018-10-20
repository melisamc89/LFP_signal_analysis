%% data load RAT14570
clear all
load('EC1_3','-mat')
%load('EC1_3_gamma','-mat')
load('EC1_3_vx','-mat')
load('EC1_3_amplitude','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        matrix_amplitude_low(:,counter)=amp_low{index}{trial};
        matrix_amplitude_high(:,counter)=amp_high{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
    end
end

n_fs=400;
init=[1,2,7,42,47,48,54,89,94]*n_fs;

for i=1:length(init)-1
    for j=1:floor(size(matrix_low,2)/50)
        x=matrix_low(init(i):init(i+1)-1,(j-1)*50+1:j*50);
        y=matrix_low(init(i)+1:init(i+1),(j-1)*50+1:j*50);
        z=matrix_high(init(i):init(i+1)-1,(j-1)*50+1:j*50);
        x=reshape(x,[size(x,1)*size(x,2) 1]);
        y=reshape(y,[size(y,1)*size(y,2) 1]);
        index=find(x<0 & y>0);
        %index=init(i)+index;
        if isempty(index)~=1
            [phase amplitude]=MakeMIHistogram(x(index(1):index(end)),z(index(1):index(end)),20);
            mi(i,j)=ModularityIndex(phase,amplitude);
            mi_phase(i,j,:)=phase;
            mi_amplitude(i,j,:)=amplitude;
        end
    end
end

number=[7,6,2,3];
for i=1:length(number)
    x=[];
    y=[];
    for j=1:size(matrix_low,2)
        index=find(matrix_low(init(number(i)):init(number(i)+1)-1,j)<0 & matrix_low(init(number(i))+1:init(number(i)+1),j)>0);
        index=init(number(i))+index;
        if isempty(index)~=1
            x=[x,matrix_low(index(1):index(end),j)'];
            y=[y,matrix_high(index(1):index(end),j)'];
        end
    end
    [phase amplitude]=MakeMIHistogram(x,y,20);
    mi_all(i)=ModularityIndex(phase,amplitude);
    mi_phase_all(i,:)=phase;
    mi_amplitude_all(i,:)=amplitude;
end

%% MI in speeds significance
mi_1(1,:)=[mi(7,:)];%,mi(7,:)];2
mi_1(2,:)=[mi(6,:)];%,mi(8,:)];4
mi_1(3,:)=[mi(2,:)];%,mi(9,:)];3
mi_1(4,:)=[mi(3,:)];%,mi(9,:)];3

subplot(2,6,[1:3])
for i=1:4
    histogram(mi_1(i,~isnan(mi_1(i,:))),'Binwidth',0.0005,'Normalization','probability')
    hold on
    mu(i)=mean(mi_1(i,~isnan(mi_1(i,:))));
    sigma(i)=var(mi_1(i,~isnan(mi_1(i,:))));
    x=0:0.0005:0.05;
    y=exp(-(x-mu(i)).*(x-mu(i))/(2*sigma(i)))/sqrt(2*pi*sigma(i));
    y=y/sum(y);
    plot(x,y,'k','Linewidth',2)
end
ylabel('Probability')
xlabel('MI')

subplot(2,6,7)
plot(mi_phase_all(1,:),reshape(mi_amplitude(7,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(1,:),'k','Linewidth',2)
subplot(2,6,8)
plot(mi_phase_all(1,:),reshape(mi_amplitude(6,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(1,:),'k','Linewidth',2)
subplot(2,6,9)
plot(mi_phase_all(1,:),reshape(mi_amplitude(2,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(1,:),'k','Linewidth',2)
subplot(2,6,10)
plot(mi_phase_all(1,:),reshape(mi_amplitude(3,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(1,:),'k','Linewidth',2)


subplot(2,6,[4,5,6])
errorbar([-10,-0.5,0.5,10],mu,sqrt(sigma))
hold on
%plot([-36,-6,0.5,0.5,6,36],mu)
plot([-10,-0.5,0.5,10],mi_all)

subplot(2,6,[11,12])
for i=1:4
    for j=1:4
        [h p]=ttest2(mi_1(i,:),mi_1(j,:));
        significance(i,j)=p;
    end
end
imagesc(log(significance)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

%% learning speed

subplot(3,3,1)
histogram(mi_1(1,~isnan(mi_1(1,1:20))),'Binwidth',0.0002,'Normalization','probability')
hold on
mu_1(1)=mean(mi_1(1,~isnan(mi_1(1,1:20))));sigma_1(1)=var(mi_1(1,~isnan(mi_1(1,1:20))));
x=0:0.0002:0.05;
y=exp(-(x-mu_1(1)).*(x-mu_1(1))/(2*sigma_1(1)))/sqrt(2*pi*sigma_1(1));y=y/sum(y);
plot(x,y,'b','Linewidth',2)
histogram(mi_1(1,end-20:end),'Binwidth',0.0002,'Normalization','probability')
hold on
mu_2(1)=mean(mi_1(1,end-20:end));sigma_2(1)=var(mi_1(1,end-20:end));
y=exp(-(x-mu_2(1)).*(x-mu_2(1))/(2*sigma_2(1)))/sqrt(2*pi*sigma_2(1));y=y/sum(y);
plot(x,y,'r','Linewidth',2)
xlim([0 0.005])

subplot(3,3,4)
histogram(mi_1(4,~isnan(mi_1(4,1:20))),'Binwidth',0.0002,'Normalization','probability')
hold on
mu_1(2)=mean(mi_1(4,~isnan(mi_1(4,1:20))));sigma_1(2)=var(mi_1(4,~isnan(mi_1(4,1:20))));
y=exp(-(x-mu_1(2)).*(x-mu_1(2))/(2*sigma_1(2)))/sqrt(2*pi*sigma_1(2));y=y/sum(y);
plot(x,y,'b','Linewidth',2)
histogram(mi_1(4,end-20:end),'Binwidth',0.0002,'Normalization','probability')
hold on
mu_2(2)=mean(mi_1(4,end-20:end));sigma_2(2)=var(mi_1(4,end-20:end));
y=exp(-(x-mu_2(2)).*(x-mu_2(2))/(2*sigma_2(2)))/sqrt(2*pi*sigma_2(2));y=y/sum(y);
plot(x,y,'r','Linewidth',2)
xlim([0 0.005])

subplot(3,3,[2,3])
errorbar([-10,10],mu_1,sqrt(sigma_1))
hold on
errorbar([-10,10],mu_2,sqrt(sigma_2))

for i=1:2
    [h p]=ttest2(mi_1(i,1:20),mi_1(i,end-20:end));
    significance(i)=p;
    for j=1:2
        [h p]=ttest2(mi_1(i,1:20),mi_1(j,1:20));
        significance1(i,j)=p;
        [h p]=ttest2(mi_1(i,end-20:end),mi_1(j,end-20:end));
        significance2(i,j)=p;
    end
end
subplot(3,3,[5 6])
imagesc(log(significance)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

subplot(3,3,8)
imagesc(log(significance1)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

subplot(3,3,9)
imagesc(log(significance2)/log(10))
colorbar
colormap(gray)
caxis([-3 0])


%% left and right
mi_1(1,:)=[mi(2,:)];%,mi(7,:)];2
mi_1(2,:)=[mi(4,:)];%,mi(8,:)];4
mi_1(3,:)=[mi(3,:)];%,mi(9,:)];3

mi_2(1,:)=mi(7,:);
mi_2(2,:)=mi(8,:);
mi_2(3,:)=mi(9,:);

for i=1:3
    subplot(3,5,[(i-1)*5+1:(i-1)*5+3])
    histogram(mi_1(i,~isnan(mi_1(i,:))),'Binwidth',0.0005,'Normalization','probability')
    hold on
    histogram(mi_2(i,~isnan(mi_2(i,:))),'Binwidth',0.0005,'Normalization','probability')
    mu_1(i)=mean(mi_1(i,~isnan(mi_1(i,:))));
    sigma_1(i)=var(mi_1(i,~isnan(mi_1(i,:))));
    mu_2(i)=mean(mi_2(i,~isnan(mi_2(i,:))));
    sigma_2(i)=var(mi_2(i,~isnan(mi_2(i,:))));
    x=0:0.0005:0.05;
    y=exp(-(x-mu_1(i)).*(x-mu_1(i))/(2*sigma_1(i)))/sqrt(2*pi*sigma_1(i));
    y=y/sum(y);
    plot(x,y,'k','Linewidth',2)
    y=exp(-(x-mu_2(i)).*(x-mu_2(i))/(2*sigma_2(i)))/sqrt(2*pi*sigma_2(i));
    y=y/sum(y);
    plot(x,y,'k','Linewidth',2)
    ylabel('Probability')
    xlabel('MI')
end


subplot(2,5,[4,5])
errorbar([0,6,36],mu_1,sigma_1)
hold on
errorbar([0,6,36],mu_2,sigma_2)

subplot(2,5,[9,10])
for i=1:3
        [h p]=ttest2(mi_1(i,:),mi_2(i,:));
        significance(i)=p;
end
imagesc(log(significance)/log(10))
colorbar
colormap(gray)


%% rat exponential RAT14566
clear all
load('EC1_1','-mat')
%load('EC1_1_gamma','-mat')
load('EC1_1_vx','-mat')
load('EC1_1_amplitude','-mat')

counter=0;
for index=1:length(signal_low)
    for trial=1:length(signal_low{index})
        counter=counter+1;
        matrix_low(:,counter)=signal_low{index}{trial};
        matrix_high(:,counter)=signal_high{index}{trial};
        matrix_amplitude_low(:,counter)=amp_low{index}{trial};
        matrix_amplitude_high(:,counter)=amp_high{index}{trial};
        matrix_pos(:,counter)=signal_pos{index}{trial};
    end
end

%n_fs=400;
%init=[1,2,7,11,15,19,23,27,31,32,38,42,46,50,54,58,62]*n_fs;

n_fs=400;
init=[1,2,7,27,31,32,38,58,62]*n_fs;

for i=1:length(init)-1
    for j=1:floor((size(matrix_low,2)-200)/10)
        x=matrix_low(init(i):init(i+1)-1,(j-1)*10+1:j*10);
        y=matrix_low(init(i)+1:init(i+1),(j-1)*10+1:j*10);
        z=matrix_high(init(i):init(i+1)-1,(j-1)*10+1:j*10);
        x=reshape(x,[size(x,1)*size(x,2) 1]);
        y=reshape(y,[size(y,1)*size(y,2) 1]);
        index=find(x<0 & y>0);
        %index=init(i)+index;
        if isempty(index)~=1
            [phase amplitude]=MakeMIHistogram(x(index(1):index(end)),z(index(1):index(end)),20);
            mi(i,j)=ModularityIndex(phase,amplitude);
            mi_phase(i,j,:)=phase;
            mi_amplitude(i,j,:)=amplitude;
        end
    end
end

for i=1:8
    subplot(2,4,i)
    plot(reshape(mi_amplitude(i,1,:),[20 1]))
    ylim([0.04 0.065])
    hold on
    plot(reshape(mi_amplitude(i,7,:),[20 1]))
end

number=[3,2,6,7];
for i=1:length(number)
    x=[];
    y=[];
    for j=1:size(matrix_low,2)
        index=find(matrix_low(init(number(i)):init(number(i)+1)-1,j)<0 & matrix_low(init(number(i))+1:init(number(i)+1),j)>0);
        index=init(number(i))+index;
        if isempty(index)~=1
            x=[x,matrix_low(index(1):index(end),j)'];
            y=[y,matrix_high(index(1):index(end),j)'];
        end
    end
    [phase amplitude]=MakeMIHistogram(x,y,20);
    mi_all(i)=ModularityIndex(phase,amplitude);
    mi_phase_all(i,:)=phase;
    mi_amplitude_all(i,:)=amplitude;
end
%% plot first part speed
    mi_1(1,:)=mi(3,:); 
    mi_1(2,:)=mi(2,:); 
    mi_1(3,:)=mi(6,:); 
    mi_1(4,:)=mi(7,:); 
    

    subplot(2,6,[1:3])
    for i=1:4
        histogram(mi_1(i,~isnan(mi_1(i,:))),'Binwidth',0.0002,'Normalization','probability')
        hold on
        mu(i)=mean(mi_1(i,~isnan(mi_1(i,:))));
        sigma(i)=var(mi_1(i,~isnan(mi_1(i,:))));
        x=0:0.0002:0.005;
        y=exp(-(x-mu(i)).*(x-mu(i))/(2*sigma(i)))/sqrt(2*pi*sigma(i));
        y=y/sum(y);
        plot(x,y,'k','Linewidth',2)
    end
    ylabel('Probability')
    xlabel('MI')

subplot(2,6,[4:6])
errorbar([-10,-1,1,10],mu,sqrt(sigma))
hold on
plot([-10,-1,1,10],mi_all)

subplot(2,6,7)
plot(mi_phase_all(1,:),reshape(mi_amplitude(3,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(1,:),'k','Linewidth',2)
subplot(2,6,8)
plot(mi_phase_all(1,:),reshape(mi_amplitude(2,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(2,:),'k','Linewidth',2)
subplot(2,6,9)
plot(mi_phase_all(1,:),reshape(mi_amplitude(6,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(3,:),'k','Linewidth',2)
subplot(2,6,10)
plot(mi_phase_all(1,:),reshape(mi_amplitude(7,:,:),[size(mi_amplitude,2) 20])')
hold on
plot(mi_phase_all(1,:),mi_amplitude_all(4,:),'k','Linewidth',2)

subplot(2,6,12)
for i=1:4
    for j=1:4
        [h p]=ttest2(mi_1(i,:),mi_1(j,:));
        significance(i,j)=p;
    end
end
imagesc(log(significance)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

%% learning

color=colormap(jet(7))
for i=1:10
    plot([-10,-1,1,10],[mi(3,i),mi(2,i),mi(6,i),mi(7,i)],'color',color(i,:))
    hold on
end



%% for learning

%subplot(5,5,[1:3])
%histogram(mi_1(1,~isnan(mi_1(1,1:75))),'Binwidth',0.0005,'Normalization','probability')
%hold on
mu_1(1)=mean(mi_1(1,~isnan(mi_1(1,1:15))));sigma_1(1)=var(mi_1(1,~isnan(mi_1(2,1:15))));
mu_2(1)=mean(mi_1(1,end-10:end));sigma_2(1)=var(mi_1(1,end-10:end));

mu_1(2)=mean(mi_1(2,~isnan(mi_1(1,1:15))));sigma_1(2)=var(mi_1(2,~isnan(mi_1(2,1:15))));
mu_2(2)=mean(mi_1(2,end-10:end));sigma_2(2)=var(mi_1(2,end-10:end));

mu_1(3)=mean(mi_1(3,~isnan(mi_1(3,1:15))));sigma_1(3)=var(mi_1(3,~isnan(mi_1(3,1:15))));
mu_2(3)=mean(mi_1(3,end-10:end));sigma_2(3)=var(mi_1(3,end-10:end));

mu_1(4)=mean(mi_1(4,~isnan(mi_1(4,1:15))));sigma_1(4)=var(mi_1(4,~isnan(mi_1(4,1:15))));
mu_2(4)=mean(mi_1(4,end-10:end));sigma_2(4)=var(mi_1(4,end-10:end));

mu_1(5)=mean(mi_1(5,~isnan(mi_1(5,1:15))));sigma_1(5)=var(mi_1(5,~isnan(mi_1(5,1:15))));
mu_2(5)=mean(mi_1(5,end-10:end));sigma_2(5)=var(mi_1(5,end-10:end));

mu_1(6)=mean(mi_1(6,~isnan(mi_1(1,1:15))));sigma_1(6)=var(mi_1(6,~isnan(mi_1(6,1:15))));
mu_2(6)=mean(mi_1(6,end-10:end));sigma_2(6)=var(mi_1(6,end-10:end));

mu_1(7)=mean(mi_1(7,~isnan(mi_1(7,1:15))));sigma_1(7)=var(mi_1(7,~isnan(mi_1(7,1:15))));
mu_2(7)=mean(mi_1(7,end-10:end));sigma_2(7)=var(mi_1(7,end-10:end));

mu_1(8)=mean(mi_1(8,~isnan(mi_1(8,1:15))));sigma_1(8)=var(mi_1(8,~isnan(mi_1(8,1:15))));
mu_2(8)=mean(mi_1(8,end-10:end));sigma_2(8)=var(mi_1(8,end-10:end));

mu_1(9)=mean(mi_1(9,~isnan(mi_1(9,1:15))));sigma_1(9)=var(mi_1(9,~isnan(mi_1(9,1:15))));
mu_2(9)=mean(mi_1(9,end-10:end));sigma_2(9)=var(mi_1(9,end-10:end));

mu_1(10)=mean(mi_1(10,~isnan(mi_1(10,1:15))));sigma_1(10)=var(mi_1(10,~isnan(mi_1(10,1:15))));
mu_2(10)=mean(mi_1(10,end-10:end));sigma_2(10)=var(mi_1(10,end-10:end));

mu_1(11)=mean(mi_1(11,~isnan(mi_1(11,1:15))));sigma_1(11)=var(mi_1(11,~isnan(mi_1(11,1:15))));
mu_2(11)=mean(mi_1(11,end-10:end));sigma_2(11)=var(mi_1(11,end-10:end));

mu_1(12)=mean(mi_1(12,~isnan(mi_1(12,1:15))));sigma_1(12)=var(mi_1(12,~isnan(mi_1(12,1:15))));
mu_2(12)=mean(mi_1(12,end-10:end));sigma_2(12)=var(mi_1(12,end-10:end));

subplot(3,5,[4,5])
errorbar([-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5],mu_1,sqrt(sigma_1))
hold on
errorbar([-5,-4,-3,-2,-1,-0.5,0.5,1,2,3,4,5],mu_2,sqrt(sigma_2))

for i=1:12
    [h p]=ttest2(mi_1(i,1:15),mi_1(i,end-10:end));
    significance(i)=p;
    for j=1:12
        [h p]=ttest2(mi_1(i,1:15),mi_1(j,1:15));
        significance1(i,j)=p;
        [h p]=ttest2(mi_1(i,end-10:end),mi_1(j,end-10:end));
        significance2(i,j)=p;
    end
end
subplot(3,5,9)
imagesc(log(significance)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

subplot(3,5,14)
imagesc(log(significance1)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

subplot(3,5,15)
imagesc(log(significance2)/log(10))
colorbar
colormap(gray)
caxis([-3 0])

subplot(3,5,13)
imagesc(significance2-significance1)
colorbar
colormap(gray)


%% left and right
mi_1(1,:)=[mi(3,:)];%,mi(7,:)];2
mi_1(2,:)=[mi(4,:)];%,mi(8,:)];4
mi_1(3,:)=[mi(5,:)];%,mi(9,:)];3
mi_1(4,:)=[mi(6,:)];%,mi(8,:)];4
mi_1(5,:)=[mi(7,:)];%,mi(9,:)];3

mi_2(1,:)=mi(15,:);
mi_2(2,:)=mi(14,:);
mi_2(3,:)=mi(13,:);
mi_2(4,:)=[mi(12,:)];%,mi(8,:)];4
mi_2(5,:)=[mi(11,:)];%,mi(9,:)];3

for i=1:5
    subplot(5,5,[(i-1)*5+1:(i-1)*5+3])
    histogram(mi_1(i,~isnan(mi_1(i,:))),'Binwidth',0.0005,'Normalization','probability')
    hold on
    histogram(mi_2(i,~isnan(mi_2(i,:))),'Binwidth',0.0005,'Normalization','probability')
    mu_1(i)=mean(mi_1(i,~isnan(mi_1(i,:))));
    sigma_1(i)=var(mi_1(i,~isnan(mi_1(i,:))));
    mu_2(i)=mean(mi_2(i,~isnan(mi_2(i,:))));
    sigma_2(i)=var(mi_2(i,~isnan(mi_2(i,:))));
    x=0:0.0005:0.05;
    y=exp(-(x-mu_1(i)).*(x-mu_1(i))/(2*sigma_1(i)))/sqrt(2*pi*sigma_1(i));
    y=y/sum(y);
    plot(x,y,'k','Linewidth',2)
    y=exp(-(x-mu_2(i)).*(x-mu_2(i))/(2*sigma_2(i)))/sqrt(2*pi*sigma_2(i));
    y=y/sum(y);
    plot(x,y,'k','Linewidth',2)
    ylabel('Probability')
    xlabel('MI')
end


subplot(5,5,[4,5])
errorbar([1,2,3,4,5],mu_1,sigma_1)
hold on
errorbar([1,2,3,4,5],mu_2,sigma_2)

subplot(5,5,[9,10])
for i=1:5
        [h p]=ttest2(mi_1(i,:),mi_2(i,:));
        significance(i)=p;
end
imagesc(log(significance)/log(10))
colorbar
colormap(gray)
