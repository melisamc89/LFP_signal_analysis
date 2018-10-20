clear all
ratas=[16262,16347,16456,16653];
count=ones(5,1);
for rat=1:length(ratas)
    %name=strcat(int2str(ratas(rat)),'_R');
    name=strcat(int2str(ratas(rat)),'_L');
    load(name,'-mat')

    for i=1:length(ecL_phase)
        for k=1:length(ecL_phase{i})
                if isempty(ecL_phase{i}{k})~=1
                    x=ecL_phase{i}{k}(1:end-1);
                    y=ecL_phase{i}{k}(2:end);    
                    index=find(x<0 & y>0);
                    if isempty(index)~=1
                        amp1{i}{count(i)}=mean(ecR_amplitude_high{i}{k}(1:index(end)));
                        amp2{i}{count(i)}=mean(ecR_amplitude_low{i}{k}(1:index(end))); 
                        count(i)=count(i)+1;
                    end
                end
        end
    end
end
    figure(rat)
    speed=[0,7,14,21,28];
    for i=2:5
        for k=1:100
            vector1(k)=amp1{i}{k};
            vector2(k)=amp1{i}{end-k};
            vector3(k)=amp2{i}{k};
            vector4(k)=amp2{i}{end-k};
        end
        subplot(6,2,(i-1)*2-1)
        a=histogram(vector1,'Binwidth',0.05,'Normalization','probability','Facecolor','b')
        ylabel('Probability')
        mu=mean(vector1);sigma=var(vector1');
        x=0:0.01:5;
        hold on
        y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
        y=y*max(a.Values)/max(y);
        plot(x,y,'Linewidth',2,'color','b')
        mu_1(i,1)=mu;sigma_1(i,1)=sigma;
        histogram(vector2,'Binwidth',0.05,'Normalization','probability','Facecolor','r')  
        mu=mean(vector2);sigma=var(vector2');
        y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
        y=y*max(a.Values)/max(y);
        plot(x,y,'Linewidth',2,'color','r')
        xlim([0 3])
        ylim([0 0.3])
        mu_1(i,2)=mu;sigma_1(i,2)=sigma;
        legend(strcat(int2str(speed(i)),' m/s F'),...
            strcat('\mu=',num2str(mu_1(i,1)),',\sigma=',num2str(sigma_1(i,1))),...
            strcat(int2str(speed(i)),' m/s L'),...
            strcat('\mu=',num2str(mu_1(i,2)),',\sigma=',num2str(sigma_1(i,2))))
        
        subplot(6,2,(i-1)*2)
        a=histogram(vector3,'Binwidth',0.05,'Normalization','probability','Facecolor','b')
        ylabel('Probability')
        mu=mean(vector3);sigma=var(vector3');
        x=0:0.01:5;
        hold on
        y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
        y=y*max(a.Values)/max(y);
        plot(x,y,'Linewidth',2,'color','b')
        mu_2(i,1)=mu;sigma_2(i,1)=sigma;
        histogram(vector4,'Binwidth',0.05,'Normalization','probability','Facecolor','r')  
        mu=mean(vector4);sigma=var(vector4');
        x=0:0.01:5;
        hold on
        y=exp(-(x-mu).*(x-mu)/(2*sigma))/sqrt(2*pi*sigma);
        y=y*max(a.Values)/max(y);
        plot(x,y,'Linewidth',2,'color','r')
        xlim([0 3])
        ylim([0 0.3])
        mu_2(i,2)=mu;sigma_2(i,2)=sigma;
        legend(strcat(int2str(speed(i)),' m/s F'),...
            strcat('\mu=',num2str(mu_2(i,1)),',\sigma=',num2str(sigma_2(i,1))),...
            strcat(int2str(speed(i)),' m/s L'),...
            strcat('\mu=',num2str(mu_2(i,2)),',\sigma=',num2str(sigma_2(i,2))))
        
        [h p1]=ttest2(vector1,vector2);
        [h p2]=ttest2(vector3,vector4);
        matrix1(i-1)=p1;
        matrix2(i-1)=p2;
        data1(i,:)=vector1;
        data2(i,:)=vector2;
        data3(i,:)=vector3;
        data4(i,:)=vector4;
    end    
    subplot(6,2,7)
    xlabel('Theta Amplitude')
    subplot(6,2,8)
    xlabel('Delta Amplitude')
    subplot(6,2,1)
    title('Theta')
    subplot(6,2,2)
    title('Delta')
    
    subplot(6,2,9)
    errorbar(speed(2:5),mu_1(2:5,1),sigma_1(2:5,1),'Linewidth',2)
    hold on
    errorbar(speed(2:5),mu_1(2:5,2),sigma_1(2:5,2),'Linewidth',2)
    xlabel('Speed (m/s)')
    ylabel('MeanAmp')
    legend('F','L')

    subplot(6,2,10)
    errorbar(speed(2:5),mu_2(2:5,1),sigma_2(2:5,1),'Linewidth',2)
    hold on
    errorbar(speed(2:5),mu_2(2:5,2),sigma_2(2:5,2),'Linewidth',2)
    xlabel('Speed (m/s)')
    ylabel('MeanAmp')
    legend('F','L')

    for i=2:5
        for j=2:5
        [h p]=ttest2(data1(i,:),data1(j,:));
        significance1_f(i-1,j-1)=p;
        [h p]=ttest2(data2(i,:),data2(j,:));
        significance1_l(i-1,j-1)=p;
        
        [h p]=ttest2(data3(i,:),data3(j,:));
        significance2_f(i-1,j-1)=p;
        [h p]=ttest2(data4(i,:),data4(j,:));
        significance2_l(i-1,j-1)=p;
        end
    end
    
    subplot(6,6,31)
    imagesc(log(reshape(matrix1,[2 2]))/log(10))
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
    imagesc(log(reshape(matrix2,[2 2]))/log(10))
    colormap(gray)
    colorbar
    title('Speed')

    subplot(6,6,35)
    imagesc(log(significance2_f)/log(10))
    colormap(gray)
    colorbar
    title('First')

    subplot(6,6,36)
    imagesc(log(significance2_l)/log(10))
    colormap(gray)
    colorbar
    title('Last')
   