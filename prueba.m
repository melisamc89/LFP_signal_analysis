for i=2:length(lfp{2})
PhaseAmplitudeCoupling2(files,rat,EC1,lfp{2}(i,:),prot,n_fs,1,20,1,100,2)
end


f1 = lowf1:fstep/2:highf1;
f2 = lowf2:fstep/2:highf2;
imagesc(f1,f2,abs(mi_4))
set(gca,'YDir','normal'), colormap(jet)
%axis xy
axis tight
xlabel('Frequency[s]');
ylabel('Frequency [Hz]');
colorbar
box on


clear Y1_1 Y2_1 Y3_1 Y4_1 Y5_1 Y6_1 Y7_1 Y8_1
clear Y1_2 Y2_2 Y3_2 Y4_2 Y5_2 Y6_2 Y7_2 Y8_2

trial=1;
Y1_1=phase1(:,srr(trial):err(trial));
Y2_1=phase1(:,srl(trial):erl(trial));
Y3_1=phase1(:,s1(trial):sr(trial));
Y4_1=phase1(:,s2(trial):sl(trial));

Y5_1=phase1(:,sr(trial):er1(trial));
Y6_1=phase1(:,el1(trial):el(trial));
Y7_1=phase1(:,er1(trial):er(trial));
Y8_1=phase1(:,sl(trial):el1(trial));

trial=1;
Y1_2=amplitude2(:,srr(trial):err(trial));
Y2_2=amplitude2(:,srl(trial):erl(trial));
Y3_2=amplitude2(:,s1(trial):sr(trial));
Y4_2=amplitude2(:,s2(trial):sl(trial));

Y5_2=amplitude2(:,sr(trial):er1(trial));
Y6_2=amplitude2(:,el1(trial):el(trial));
Y7_2=amplitude2(:,er1(trial):er(trial));
Y8_2=amplitude2(:,sl(trial):el1(trial));

for trial=2:length(sl)-1
    
        Y1_1=cat(2,Y1_1,phase1(:,srr(trial):err(trial)));
        Y2_1=cat(2,Y2_1,phase1(:,srl(trial):erl(trial)));
        Y3_1=cat(2,Y3_1,phase1(:,s1(trial):sr(trial)));
        Y4_1=cat(2,Y4_1,phase1(:,s2(trial):sl(trial)));

        Y5_1=cat(2,Y5_1,phase1(:,sr(trial):er1(trial)));
        Y6_1=cat(2,Y6_1,phase1(:,el1(trial):el(trial)));
        Y7_1=cat(2,Y7_1,phase1(:,er1(trial):er(trial)));
        Y8_1=cat(2,Y8_1,phase1(:,sl(trial):el1(trial)));

        Y1_2=cat(2,Y1_2,amplitude2(:,srr(trial):err(trial)));
        Y2_2=cat(2,Y2_2,amplitude2(:,srl(trial):erl(trial)));
        Y3_2=cat(2,Y3_2,amplitude2(:,s1(trial):sr(trial)));
        Y4_2=cat(2,Y4_2,amplitude2(:,s2(trial):sl(trial)));

        Y5_2=cat(2,Y5_2,amplitude2(:,sr(trial):er1(trial)));
        Y6_2=cat(2,Y6_2,amplitude2(:,el1(trial):el(trial)));
        Y7_2=cat(2,Y7_2,amplitude2(:,er1(trial):er(trial)));
        Y8_2=cat(2,Y8_2,amplitude2(:,sl(trial):el1(trial)));
end

signal_phase_non_mov=cat(2,Y1_1,Y2_1);
signal_phase_wait=cat(2,Y3_1,Y4_1);
signal_phase_fast=cat(2,Y5_1,Y6_1);
signal_phase_slow=cat(2,Y7_1,Y8_1);

signal_amp_non_mov=cat(2,Y1_2,Y2_2);
signal_amp_wait=cat(2,Y3_2,Y4_2);
signal_amp_fast=cat(2,Y5_2,Y6_2);
signal_amp_slow=cat(2,Y7_2,Y8_2);


 mi_1=zeros(count1-1,count2-1);
for i=1:count1-1
     for j=1:count2-1
         [f_phases mean_amplitude]=DistributionforMudularityIndex(signal_phase_non_mov(i,:),signal_amp_non_mov(j,:),40);
         mi_1(i,j)=ModularityIndex(f_phases,mean_amplitude);
     end
end

 mi_2=zeros(count1-1,count2-1);
for i=1:count1-1
     for j=1:count2-1
         [f_phases mean_amplitude]=DistributionforMudularityIndex(signal_phase_wait(i,:),signal_amp_wait(j,:),40);
         mi_2(i,j)=ModularityIndex(f_phases,mean_amplitude);
     end
end

mi_3=zeros(count1-1,count2-1);
for i=1:count1-1
     for j=1:count2-1
         [f_phases mean_amplitude]=DistributionforMudularityIndex(signal_phase_fast(i,:),signal_amp_fast(j,:),40);
         mi_3(i,j)=ModularityIndex(f_phases,mean_amplitude);
     end
end

 mi_4=zeros(count1-1,count2-1);
for i=1:count1-1
     for j=1:count2-1
         [f_phases mean_amplitude]=DistributionforMudularityIndex(signal_phase_slow(i,:),signal_amp_slow(j,:),40);
         mi_4(i,j)=ModularityIndex(f_phases,mean_amplitude);
     end
end
