
load('16262_L','-mat')

for k=1:5
signal1{k}=[];
signal2{k}=[];
signal3{k}=[];
signal4{k}=[];
end

for i=1:length(ecL_phase)
    for k=1:length(ecL_phase{i})
            if isempty(ecL_phase{i}{k})~=1
                x=ecL_phase{i}{k}(1:end-1);
                y=ecL_phase{i}{k}(2:end);    
                index=find(x<0 & y>0);
                if isempty(index)~=1
                    signal1{i}=[signal1{i} ecL_phase{i}{k}(1:index(end))'];
                    signal2{i}=[signal2{i} ecL_amplitude{i}{k}(1:index(end))']; 
                    signal3{i}=[signal3{i} ecR_amplitude_low{i}{k}(1:index(end))']; 
                    signal4{i}=[signal4{i} ecR_amplitude_high{i}{k}(1:index(end))']; 
                end
            end
    end
end

for i=1:length(signal1)
    [phase amplitude]=MakeMIHistogram(signal1{i},signal2{i},20);
    mi(i)=ModularityIndex(phase,amplitude);
    amp_low(i)=mean(signal3{i});
    amp_high(i)=mean(signal4{i});
end


for i=1:length(signal1)
    subplot(3,2,i)
    histogram(signal3{i},'Binwidth',0.1,'Normalization','probability')
    hold on
    histogram(signal4{i},'Binwidth',0.1,'Normalization','probability')
    xlim([0 5])
end

%% Analysis with pooling

clear signal1 signal2 signal3 signal4
pools=9;

for l=1:3
    for k=1:5
    signal1{k}{l}=[];
    signal2{k}{l}=[];
    signal3{k}{l}=[];
    signal4{k}{l}=[];
    end
    for m=1:9
        i=(l-1)*9+m;
        for k=1:length(ecL_phase{i})
            for j=1:length(ecL_phase{i}{k})
                if isempty(ecL_phase{i}{k}{j})~=1
                    x=ecL_phase{i}{k}{j}(1:end-1);
                    y=ecL_phase{i}{k}{j}(2:end);    
                    index=find(x<0 & y>0);
                    if isempty(index)~=1
                        signal1{k}{l}=[signal1{k}{l} ecL_phase{i}{k}{j}(1:index(end))'];
                        signal2{k}{l}=[signal2{k}{l} ecL_amplitude{i}{k}{j}(1:index(end))']; 
                        signal3{k}{l}=[signal3{k}{l} ecR_amplitude_low{i}{k}{j}(1:index(end))']; 
                        signal4{k}{l}=[signal4{k}{l} ecR_amplitude_high{i}{k}{j}(1:index(end))']; 
                    end
                end
            end
        end
    end
end

for i=1:5
    for j=1:3
        mean_amplitude(i,j)=mean(signal4{i}{j});
    end
end

for i=1:5
    for j=1:3
        [phase amplitude]=MakeMIHistogram(signal1{i}{j},signal2{i}{j},20);
        mi(i,j)=ModularityIndex(phase,amplitude);
    end
end
    