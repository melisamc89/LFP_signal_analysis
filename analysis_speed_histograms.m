clear all
load('EC1_3','-mat')
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
init=[1,2,7,12,42,47,48,54,84,89,94]*n_fs;

%rat=1
%n_fs=400;
%init=[1,2,7,11,15,19,23,27,31,32,38,42,46,50,54,58,62]*n_fs;


for i=1:length(init)-1
    signal1{i}=matrix_amplitude_low(init(i):init(i+1)-1,:);
    signal2{i}=matrix_amplitude_high(init(i):init(i+1)-1,:);
end


signal1_1{1}=[signal1{1};signal1{6}];
signal1_1{2}=[signal1{2};signal1{7}];
signal1_1{3}=[signal1{4};signal1{8}];
signal1_1{4}=[signal1{3};signal1{9}];
signal1_1{5}=[signal1{5};signal1{10}];


color=['b','m','r','g','y'];
for i=1:length(signal1_1)
    subplot(2,6,[1:3])
    data1{i}=mean(signal2_1{i});
    a=histogram(mean(signal2_1{i}),'Binwidth',0.05,'Normalization','probability')
    mu=mean(mean(signal2_1{i}));sigma=var(mean(signal2_1{i}));
    x=0:0.05:5;
    y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
    y=y/sum(y);
    hold on
    plot(x,y,'Linewidth',2,'color',color(i))
    mu_1(i)=mu;
    var_1(i)=sigma;
    title('Theta Amplitude')
    xlabel('Amplitude')
    ylabel('Probability')
    
    subplot(2,6,[7:9])
    data2{i}=mean(signal2_1{i});
    a=histogram(mean(signal1_1{i}),'Binwidth',0.05,'Normalization','probability')
    mu=mean(mean(signal1_1{i}));sigma=var(mean(signal1_1{i}));
    x=0:0.05:5;
    hold on
    y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
    y=y/sum(y);
    plot(x,y,'Linewidth',2,'color',color(i))
    mu_2(i)=mu;
    var_2(i)=sigma;
    title('Delta amplitude')
        xlabel('Amplitude')
    ylabel('Probability')
end
speed1={'0_s_t_i_l_l',strcat('\mu=',num2str(mu_1(1))),'0_a_l_a_r_m',strcat('\mu=',num2str(mu_1(2))),...
        '6',strcat('\mu=',num2str(mu_1(3))),'36',strcat('\mu=',num2str(mu_1(4))),...
        '0_a_f_t_e_r',strcat('\mu=',num2str(mu_1(5)))};
speed2={'0_s_t_i_l_l',strcat('\mu=',num2str(mu_2(1))),'0_a_l_a_r_m',strcat('\mu=',num2str(mu_2(2))),...
        '6',strcat('\mu=',num2str(mu_2(3))),'36',strcat('\mu=',num2str(mu_2(4))),...
        '0_a_f_t_e_r',strcat('\mu=',num2str(mu_2(5)))};
subplot(2,6,[1:3])
legend(speed1)
subplot(2,6,[7:9])
legend(speed2)

    speeds=[0.1,0.2,6,36,0.3];
    subplot(2,6,4)
    errorbar(speeds,mu_1,var_1,'Linewidth',2)
    xlabel('Speed (m/s)')
    ylabel('MeanAmp')
    
    subplot(2,6,10)
    errorbar(speeds,mu_2,var_2,'Linewidth',2)
    xlabel('Speed (m/s)')
    ylabel('MeanAmp')
    
    for i=1:5
       for j=1:5
           [h p]=ttest2(sort(data1{i}),sort(data1{j}));
           significance1(i,j)=p;
           [h p]=ttest2(sort(data2{i}),sort(data2{j}));
           significance2(i,j)=p;
       end
    end
    subplot(2,6,[5:6])
    imagesc(log(significance1)/log(10))
    colormap(gray)
    colorbar
    caxis([-3 0])
    title('Significance')
    %caxis([0 0.05])
   
    subplot(2,6,[11:12])
    imagesc(log(significance2)/log(10))
    colormap(gray)
    caxis([-3 0])
    colorbar
    title('Significance')
    
    %% learning
       
signal1_1{1}=[signal1{1};signal1{6}];
signal1_1{2}=[signal1{2};signal1{7}];
signal1_1{3}=[signal1{4};signal1{8}];
signal1_1{4}=[signal1{3};signal1{9}];
signal1_1{5}=[signal1{5};signal1{10}];

signal2_1{1}=[signal2{1};signal2{6}];
signal2_1{2}=[signal2{2};signal2{7}];
signal2_1{3}=[signal2{4};signal2{8}];
signal2_1{4}=[signal2{3};signal2{9}];
signal2_1{5}=[signal2{5};signal2{10}];
    
    
states={'Still','Alarm','6 m/s','36 m/s'};
for i=1:length(signal1_1)-1
    subplot(6,2,i*2-1)
    data1{i}=mean(signal2_1{i});
    a=histogram(data1{i}(1:100),'Binwidth',0.05,'Normalization','probability')
    mu=mean(data1{i}(1:100));sigma=var(data1{i}(1:100));
    x=0:0.05:5;
    y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
    y=y/sum(y);
    hold on
    plot(x,y,'Linewidth',2,'color','b')
    mu_1(i,1)=mu;
    var_1(i,1)=sigma;
    
    a=histogram(data1{i}(end-100:end),'Binwidth',0.05,'Normalization','probability')
    mu=mean(data1{i}(end-100:end));sigma=var(data1{i}(end-100:end));
    x=0:0.05:5;
    y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
    y=y/sum(y);
    hold on
    plot(x,y,'Linewidth',2,'color','r')
    mu_1(i,2)=mu;
    var_1(i,2)=sigma;
    xlim([0 3])
    %title('Theta Amplitude')
    %xlabel('Amplitude')
    %ylabel('Probability')
    title(states{i})
    legend('F',strcat('\mu=',num2str(mu_1(i,1)),',\sigma=',num2str(var_1(i,1))),'L',...
        strcat('\mu=',num2str(mu_1(i,2)),',\sigma=',num2str(var_1(i,2))))

    
    subplot(6,2,(i*2))
    data2{i}=mean(signal1_1{i});
    a=histogram(data2{i}(1:100),'Binwidth',0.05,'Normalization','probability')
    mu=mean(data2{i}(1:100));sigma=var(data2{i}(1:100));
    hold on
    y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
    y=y/sum(y);
    plot(x,y,'Linewidth',2,'color','b')
    mu_2(i,1)=mu;
    var_2(i,1)=sigma;
    
    data2{i}=mean(signal1_1{i});
    a=histogram(data2{i}(end-100:end),'Binwidth',0.05,'Normalization','probability')
    mu=mean(data2{i}(end-100:end));sigma=var(data2{i}(end-100:end));
    hold on
    y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
    y=y/sum(y);
    plot(x,y,'Linewidth',2,'color','r')
    mu_2(i,2)=mu;
    var_2(i,2)=sigma;
    xlim([0 3])
    title(states{i})
    legend('F',strcat('\mu=',num2str(mu_2(i,1)),',\sigma=',num2str(var_2(i,1))),'L',...
        strcat('\mu=',num2str(mu_2(i,2)),',\sigma=',num2str(var_2(i,2))))
    %title('Delta amplitude')
     %   xlabel('Amplitude')
    %ylabel('Probability')
end

subplot(6,2,9)
speeds=[0.1,0.2,6,36];
errorbar(speeds,mu_1(:,1),var_1(:,1),'Linewidth',2)
hold on
errorbar(speeds,mu_1(:,2),var_1(:,2),'Linewidth',2)
xlabel('Speed (m/s)')
ylabel('MeanAmp')
legend('F','L')

subplot(6,2,10)
speeds=[0.1,0.2,6,36];
errorbar(speeds,mu_2(:,1),var_2(:,1),'Linewidth',2)
hold on
errorbar(speeds,mu_2(:,2),var_2(:,2),'Linewidth',2)
xlabel('Speed (m/s)')
ylabel('MeanAmp')
legend('F','L')

for i=2:4
    [h p]=ttest2(data1{i}(1:100),data1{i}(end-100:end));    
    matrix1(i-1)=p;
    [h p]=ttest2(data2{i}(1:100),data2{i}(end-100:end));    
    matrix2(i-1)=p;    
   for j=2:4
        [h p]=ttest2(data1{i}(1:100),data1{j}(1:100));
        significance1_f(i-1,j-1)=p;
        [h p]=ttest2(data1{i}(end-100:end),data1{j}(end-100:end));
        significance1_l(i-1,j-1)=p;
        
        [h p]=ttest2(data2{i}(1:100),data2{j}(1:100));
        significance2_f(i-1,j-1)=p;
        [h p]=ttest2(data2{i}(end-100:end),data2{j}(end-100:end));
        significance2_l(i-1,j-1)=p;
   end
end

subplot(6,6,31)
imagesc(log(matrix1)/log(10))
colormap(gray)
colorbar
title('Speed')

subplot(6,6,32)
imagesc(log(significance1_f)/log(10))
colormap(gray)
colorbar
title('First')

subplot(6,6,33)
imagesc(log(significance1_l)/log(10))
colormap(gray)
colorbar
title('Last')


subplot(6,6,34)
imagesc(log(matrix2)/log(10))
colormap(gray)
colorbar
caxis([-3 0])
title('Speed')

subplot(6,6,35)
imagesc(log(significance2_f)/log(10))
colormap(gray)
colorbar
caxis([-3 0])
title('First')

subplot(6,6,36)
imagesc(log(significance2_l)/log(10))
colormap(gray)
colorbar
caxis([-3 0])
title('Last')

%% leaerning ks

   
states={'Still','Alarm','6 m/s','36 m/s'};
for i=1:length(signal1_1)-1
    subplot(6,2,i*2-1)
    data1{i}=mean(signal2_1{i});
    scatter(1:100,sort(data1{i}(1:100)),'Filled')
    hold on
    scatter(1:100,sort(data1{i}(end-99:end)),'Filled')
    title(states{i})
    legend('F','L')
    
    subplot(6,2,(i*2))
     scatter(1:100,sort(data2{i}(1:100)),'Filled')
    hold on
    scatter(1:100,sort(data2{i}(end-99:end)),'Filled')
    title(states{i})
    legend('F','L')
end

for i=1:4
    [h p]=kstest2(data1{i}(1:100),data1{i}(end-100:end));    
    matrix1(i)=p;
    [h p]=kstest2(data2{i}(1:100),data2{i}(end-100:end));    
    matrix2(i)=p;    
   for j=1:4
        [h p]=kstest2(data1{i}(1:100),data1{j}(1:100));
        significance1_f(i,j)=p;
        [h p]=kstest2(data1{i}(end-100:end),data1{j}(end-100:end));
        significance1_l(i,j)=p;
        
        [h p]=kstest2(data2{i}(1:100),data2{j}(1:100));
        significance2_f(i,j)=p;
        [h p]=kstest2(data2{i}(end-100:end),data2{j}(end-100:end));
        significance2_l(i,j)=p;
   end
end

subplot(6,6,[25,31])
imagesc(log(reshape(matrix1,[2 2]))/log(10))
colormap(gray)
colorbar
title('Speed')

subplot(6,6,[26,32])
imagesc(log(significance1_f)/log(10))
colormap(gray)
colorbar
title('First')

subplot(6,6,[27,33])
imagesc(log(significance1_l)/log(10))
colormap(gray)
colorbar
title('Last')


subplot(6,6,[28,34])
imagesc(log(reshape(matrix2,[2 2]))/log(10))
colormap(gray)
colorbar
title('Speed')

subplot(6,6,[29,35])
imagesc(log(significance2_f)/log(10))
colormap(gray)
colorbar
title('First')

subplot(6,6,[30,36])
imagesc(log(significance2_l)/log(10))
colormap(gray)
colorbar
title('Last')
