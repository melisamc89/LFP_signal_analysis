%% Before starting, load everything

clear all
close all
files=Data_Listing();

n_fs=400;
sf=50;
n_osc=0.5;

load('1CA1','-mat');
load('1EC1','-mat');
load('1EC2','-mat');
rat=1;
area={CA1,EC1,EC2};
lfp={Days(CA1),Days(EC1),Days(EC2)};
name=['1CA1';'1EC1';'1EC2'];


load('2CA1','-mat');
load('2EC1','-mat');
load('2EC2','-mat');
rat=2;
area={CA1_2,EC1_2,EC2_2};
lfp={Days(CA1_2),Days(EC1_2),Days(EC2_2)};
name=['2CA1';'2EC1';'2EC2'];

load('3CA1','-mat');
load('3EC1','-mat');
load('3EC2','-mat');
rat=3;
area={CA1_3,EC1_3,EC2_3};
lfp={Days(CA1_3),Days(EC1_3),Days(EC2_3)};
name=['3CA1';'3EC1';'3EC2'];

protocol={'A','B','C'};
[n m]=size(area);
a=2;
count=1;
inf_limit=[0,5,15,20,25];
sup_limit=[5,15,20,25,50];
for prot=1:2
for i=1:60
    st={0,0,0,0,0};
    et={0,0,0,0,0};
    session_number=RecordingSession(files,rat,area{a},lfp{a}(i,:),protocol{prot});
       if length(session_number)
           for se=1:length(session_number)
           order=Wished_Register_Order(area{a},lfp{a}(i,:));
           register=Rat_Register(files,area{a},rat,order);
           [eeg,fs]=ReadEEG(register,session_number(se));
           [qual4,qual,qual2]=egf_theta_quality_2016(eeg);
           if qual>0.5
                 n_fs=floor(fs/floor(fs/n_fs));
                 eeg=resample(eeg,1,floor(fs/n_fs)); 
                 pos=Loading_Pos(register,session_number(se));
                 vel=pos.data.v;
                 time=floor(pos.data.t*n_fs)+1;
                 counter=zeros(5,1);
                 episodes=ones(5,1);
                 for k=1:length(vel)
                     for limit=1:5
                         if vel(k)>inf_limit(limit) && vel(k)<sup_limit(limit)
                                if counter(limit)==0
                                    last=find(counter);
                                    if counter(last)>n_osc*sf;
                                         st{last}(episodes(last))=start;
                                         et{last}(episodes(last))=time(k);
                                         episodes(last)=episodes(last)+1;
                                    end
                                    counter=zeros(5,1);
                                    counter(limit)=1;
                                    start=time(k);
                                else
                                    counter(limit)=counter(limit)+1;
                                end
                         end
                     end
                 end
                 for limit=1:5
                     Y=[];
                     if (st{limit})
                     for ep=1:length(st{limit})
                         x1=st{limit}(ep); x2=et{limit}(ep);
                         Y=[Y,eeg(x1:x1+(x2-x1))'];
                     end
                     end
                     n = 2^10;
                     sig=fft(Y,n);
                     L=length(sig);
                     f = n_fs*(0:(L-1))/L;
                     matrix(count,limit,:)=abs(sig/sqrt(L));   
                 end
                 count=count+1;
              end
          end
       end
end
end

c = colormap(jet(5));
[n1 n2 n3]=size(matrix);
L=n3;
for i=1:5
    P2 = mean(matrix(:,i,:));
    vP2=sqrt(var(matrix(:,i,:)));
    mP2=vP2/sqrt(n1);
    P1 = P2(1,1,1:L/2+1);
    vP1=vP2(1,1,1:L/2+1);
    mP1=mP2(1,1,1:L/2+1);
    P1(1,1,2:end-1)= 2*P1(1,1,2:end-1);
    f = n_fs*(0:(L/2))/L;
    %errorbar(f,P1,vP1,'Color',c(i,:))
    %hold on
    %errorbar(f,P1,mP1,'Color',[128 128 128]/225)
    plot(f,reshape(P1,[length(P1) 1]),'Color',c(i,:),'Linewidth',2);
    xlabel('Frequency [Hz]','FontSize',15);
    ylabel('Mean PSD','FontSize',15);
    axis([0 150 0 max(P1(10:end))])
    hold on
    box on
end
axes('Position',[.5 0.6 .3 .3])
    box on
for i=1:5
    P2 = mean(matrix(:,i,:));
    vP2=sqrt(var(matrix(:,i,:)));
    mP2=vP2/sqrt(n1);
    P1 = P2(1,1,1:L/2+1);
    vP1=vP2(1,1,1:L/2+1);
    mP1=mP2(1,1,1:L/2+1);
    P1(1,1,2:end-1)= 2*P1(1,1,2:end-1);
    f = n_fs*(0:(L/2))/L;
    %errorbar(f,P1,vP1,'Color',[192 192 192]/225)
    %errorbar(f,P1,mP1,'Color',[128 128 128]/225)
    plot(f,reshape(P1,[length(P1) 1]),'Color',c(i,:),'Linewidth',2);
    xlabel('Frequency [Hz]','FontSize',10);
    ylabel('Mean PSD','FontSize',10);
        hold on
            axis([0.5 15 0 max(P1(10:end))])

end

 axes('Position',[.5 0.2 .3 .3])
    box on
    
for i=1:5
   P2 = mean(matrix(:,i,:));
    vP2=sqrt(var(matrix(:,i,:)));
    mP2=vP2/sqrt(n1);
    P1 = P2(1,1,1:L/2+1);
    vP1=vP2(1,1,1:L/2+1);
    mP1=mP2(1,1,1:L/2+1);
    P1(1,1,2:end-1)= 2*P1(1,1,2:end-1);
    f = n_fs*(0:(L/2))/L;
    %errorbar(f,P1,vP1,'Color',[192 192 192]/225)
    %errorbar(f,P1,mP1,'Color',[128 128 128]/225)
    plot(f,reshape(P1,[length(P1) 1]),'Color',c(i,:),'Linewidth',2);
    xlabel('Frequency [Hz]','FontSize',10);
    ylabel('Mean PSD','FontSize',10);
        hold on
            axis([40 100 0.05 0.3])

end

