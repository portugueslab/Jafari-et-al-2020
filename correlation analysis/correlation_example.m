clear all
clc

% Each mouse was imaged at 4 timepoints. The same three sites were imaged
% at each timpoint and each site was imaged three times in each timepoint.
% Each timepoint has its own excel file and within each file there is a
% sheet for each site and repetition: the site is the second number and the
% imaging repetition the fourth. 

fname0='sample1_sysEAE.xlsx';
fname1='sample1_cyto+3.xlsx';
fname2='sample1_cyto+10.xlsx';
fname3='sample1_cyto+17.xlsx';
mtitle='sample1_cMS';

%%%%%%%%%%%%%%%%%%%%%%%% fname0 

[status,sheets]=xlsfinfo(fname0);
numsite=0;
sessions=0;
for j=1:length(sheets)
    Q=char(sheets{j});
    numsite=max(str2num(Q(2)),numsite);
    sessions=max(str2num(Q(4)),sessions);
end
numsite
sessions

% We now go through the different imaging sites and see
% what is the maximum number of cells recorded in each site across all
% three repetitions. This number is m1(i) where i labels the sites.

for i=1:numsite
    m1(i)=0;
    for TP=1:4
        for j=1:sessions
            tSITE=(i-1)*sessions+j;
            M=xlsread(fname0,char(sheets{tSITE}));
            if size(M,1)>0
                m1(i)=max(m1(i),max(M(:,1)));
            end
        end
    end
end
m1;

% All the data for this mouse will be saved in a structure 
% SITE({1,2,3}).TP{1,2,3,4}.session{1,2,3}
% where SITE referes to location
% TP to timepoint
% session to repetition within location and timepoint
% with fields 
% exist = 1 if the experiment exists (0 otherwise)
% corr is the correlation matrix
% SITE_id is a list of the rois that were identified in this imaging
%       session

for i=1:numsite
    toremove=[];
    for j=1:sessions
        tSITE=(i-1)*sessions+j;
        M=xlsread(fname0,char(sheets{tSITE}));   
        % M is the correlation matrix from the excel sheet, 
        % where the first row and column labels the id of the SITE
        clear tcorr
        SITE(i).TP(1).session(j).exist=0;
        SITE_id=[];
        tcorr=zeros([m1(i),m1(i)]);
        if size(M,1)>0
            for m=2:size(M,1)
                for n=2:size(M,1)
                    tcorr(M(1,m),M(1,n))=M(m,n);
                    if m==n
                        tcorr(M(1,m),M(1,n))=1;
                    end
                end
            end
            
            for m=1:size(tcorr,1)
                if sum(squeeze(tcorr(m,:)))==0
                    toremove=[toremove m];
                else
                    SITE_id=[SITE_id m];
                end
            end
            SITE(i).TP(1).session(j).corr=tcorr;
            SITE(i).TP(1).session(j).exist=1;
            SITE(i).TP(1).session(j).SITEs=SITE_id;
        end
    end
    toremove=unique(toremove);
end

%%%%%%%%%%%%%%%%%%%%%%%% fname1 

[status,sheets]=xlsfinfo(fname1);
numsite=0;
sessions=0;
for j=1:length(sheets)
    Q=char(sheets{j});
    numsite=max(str2num(Q(2)),numsite);
    sessions=max(str2num(Q(4)),sessions);
end
numsite
sessions

for i=1:numsite
    toremove=[];
    for j=1:sessions
        tSITE=(i-1)*sessions+j;
        M=xlsread(fname1,char(sheets{tSITE}));
        clear tcorr
        SITE(i).TP(2).session(j).exist=0;
        SITE_id=[];
        tcorr=zeros([m1(i),m1(i)]);
        if size(M,1)>0
            for m=2:size(M,1)
                for n=2:size(M,1)
                    tcorr(M(1,m),M(1,n))=M(m,n);
                    if m==n
                        tcorr(M(1,m),M(1,n))=1;
                    end
                end
            end
            
            for m=1:size(tcorr,1)
                if sum(squeeze(tcorr(m,:)))==0
                    toremove=[toremove m];
                else
                    SITE_id=[SITE_id m];
                end
            end
            SITE(i).TP(2).session(j).corr=tcorr;
            SITE(i).TP(2).session(j).exist=1;
            SITE(i).TP(2).session(j).SITEs=SITE_id;
        end
    end
    toremove=unique(toremove);
end


%%%%%%%%%%%%%%%%%%%%%%%% fname2 

[status,sheets]=xlsfinfo(fname2);
numsite=0;
sessions=0;
for j=1:length(sheets)
    Q=char(sheets{j});
    numsite=max(str2num(Q(2)),numsite);
    sessions=max(str2num(Q(4)),sessions);
end
numsite
sessions

for i=1:numsite
    toremove=[];
    for j=1:sessions
        tSITE=(i-1)*sessions+j;
        M=xlsread(fname2,char(sheets{tSITE}));
        clear tcorr
        SITE(i).TP(3).session(j).exist=0;
        SITE_id=[];
        tcorr=zeros([m1(i),m1(i)]);
        if size(M,1)>0
            for m=2:size(M,1)
                for n=2:size(M,1)
                    tcorr(M(1,m),M(1,n))=M(m,n);
                    if m==n
                        tcorr(M(1,m),M(1,n))=1;
                    end
                end
            end
            
            for m=1:size(tcorr,1)
                if sum(squeeze(tcorr(m,:)))==0
                    toremove=[toremove m];
                else
                    SITE_id=[SITE_id m];
                end
            end
            SITE(i).TP(3).session(j).corr=tcorr;
            SITE(i).TP(3).session(j).exist=1;
            SITE(i).TP(3).session(j).SITEs=SITE_id;
        end
    end
    toremove=unique(toremove);
end


%%%%%%%%%%%%%%%%%%%%%%%% fname3 

[status,sheets]=xlsfinfo(fname3);
numsite=0;
sessions=0;
for j=1:length(sheets)
    Q=char(sheets{j});
    numsite=max(str2num(Q(2)),numsite);
    sessions=max(str2num(Q(4)),sessions);
end
numsite
sessions

for i=1:numsite
    toremove=[];
    for j=1:sessions
        tSITE=(i-1)*sessions+j;
        M=xlsread(fname3,char(sheets{tSITE}));
        clear tcorr
        SITE(i).TP(4).session(j).exist=0;
        SITE_id=[];
        tcorr=zeros([m1(i),m1(i)]);
        if size(M,1)>0
            for m=2:size(M,1)
                for n=2:size(M,1)
                    tcorr(M(1,m),M(1,n))=M(m,n);
                    if m==n
                        tcorr(M(1,m),M(1,n))=1;
                    end
                end
            end
            
            for m=1:size(tcorr,1)
                if sum(squeeze(tcorr(m,:)))==0
                    toremove=[toremove m];
                else
                    SITE_id=[SITE_id m];
                end
            end
            SITE(i).TP(4).session(j).corr=tcorr;
            SITE(i).TP(4).session(j).exist=1;
            SITE(i).TP(4).session(j).SITEs=SITE_id;
        end
    end
    toremove=unique(toremove);
end

figure('Renderer','painters','Position',[100 100 1400 600],'NumberTitle', 'off', 'Name', mtitle)

for nsite=1:numsite     % labels the site
    for i=1:sessions    % labels the repetition 
        for j=1:4       % labels the timepoint
            if SITE(nsite).TP(j).session(i).exist==1
                Q=SITE(nsite).TP(j).session(i).corr;
                subplot(6,14,5*(nsite-1)+j+14*(i-1))
                imagesc(Q)
                size(Q)
                for a=1:size(Q,1)
                    for b=1:size(Q,2)
                        if a>=b
                            Q(a,b)=1000;
                        end
                    end
                end
                A2=reshape(Q,[1 size(Q,1)*size(Q,2)]);
                A2(A2==1000)=[];
                SITE(nsite).TP(j).session(i).corr_vector=A2;
                drawnow
            end
        end
    end
    
end


for nsite=1:numsite
    for k=1:4
        values(k).correlations=[];
    end
    for j=1:1
        for k=1:4
            for i=1:sessions
                for i2=1:sessions
                    if k==1
                        if i2>i
                            if SITE(nsite).TP(j).session(i).exist*SITE(nsite).TP(k).session(i2).exist==1
%                                 C1=setdiff(SITE(nsite).TP(j).session(i).SITEs,SITE(nsite).TP(k).session(i2).SITEs);
%                                 idx=[];
%                                 for q=1:length(C1)
%                                     idx=[idx find(SITE(nsite).TP(j).session(i).SITEs==C1(q))];
%                                 end
%                                 C2=setdiff(SITE(nsite).TP(k).session(i2).SITEs,SITE(nsite).TP(j).session(i).SITEs);
%                                 idx2=[];
%                                 for q=1:length(C2)
%                                     idx2=[idx2 find(SITE(nsite).TP(k).session(i2).SITEs==C2(q))];
%                                 end
                                Q1=SITE(nsite).TP(j).session(i).corr;
                                Q2=SITE(nsite).TP(k).session(i2).corr;
%                                 Q1(idx,:)=[];
%                                 Q1(:,idx)=[];
%                                 Q2(idx2,:)=[];
%                                 Q2(:,idx2)=[];
                                if sum(size(Q1)==size(Q2))==2
                                    A=SITE(nsite).TP(j).session(i).corr_vector;
                                    B=SITE(nsite).TP(k).session(i2).corr_vector;
                                    C=corrcoef(A,B);
                                    values(k).correlations=[values(k).correlations C(1,2)];
                                end
                            end
                        end
                    else
                        if SITE(nsite).TP(j).session(i).exist*SITE(nsite).TP(k).session(i2).exist==1
                            Q1=SITE(nsite).TP(j).session(i).corr;
                            Q2=SITE(nsite).TP(k).session(i2).corr;
                            if sum(size(Q1)==size(Q2))==2
                                A=SITE(nsite).TP(j).session(i).corr_vector;
                                B=SITE(nsite).TP(k).session(i2).corr_vector;
                                C=corrcoef(A,B);
                                values(k).correlations=[values(k).correlations C(1,2)];
                            end
                        end
                    end
                end
            end
        end
    end
    if nsite==1
        subplot(6,14,[43 44 45 46 57 58 59 60 71 72 73 74])
    end
    if nsite==2
        subplot(6,14,[48 49 50 51 62 63 64 65 76 77 78 79])
    end
    if nsite==3
        subplot(6,14,[53 54 55 56 67 68 69 70 81 82 83 84])
    end
    
    for ii=1:4
        errorbar(ii,nanmean(values(ii).correlations),nanstd(values(ii).correlations)/sqrt(length(values(ii).correlations)),'k')
        datas(nsite,ii)=nanmean(values(ii).correlations);
        hold on
        xlim([0 5])
        ylim([0 1])
    end
    drawnow
end

subplot(6,14,[43 44 45 46 57 58 59 60 71 72 73 74])
ylabel('Average correlation value Site 1')
xlabel('Timepoint')

subplot(6,14,[48 49 50 51 62 63 64 65 76 77 78 79])
ylabel('Average correlation value Site 2')
xlabel('Timepoint')

subplot(6,14,[53 54 55 56 67 68 69 70 81 82 83 84])
ylabel('Average correlation value Site 3')
xlabel('Timepoint')


