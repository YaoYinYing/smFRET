%% To analyze FRET data in a rigorous way. ALEX-like data generation and anlysis.
%  a 10r-ng-10r illumination pattern should be strictly followed.

clear;
close all;
fclose('all');

%% Read data, a movie at a time
pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
   	pth='C:\User\tir data\yyyy\New Folder';
end
cd(pth);
fList=dir;
nf=size(fList,1);

%three baseline intensities (use s_tr.m to get these baseline intensities)
donorBlank=5.0;
acceptorBlank=5.0;
donor2Blank=0.0;
acceptor2Blank=0.0;
%theta
reflection=0.075;
%alpha
leakage=0.13;
%delta
directEx=0.0;
%gamma
detection=1.0;
%beta
normalEx=1.0;

TOTAL=[];
intensities1_10=[];
% for Chris
intensityTOTAL11_20=[];
EDGES={(-0.5:0.025:1.5),(-0.5:0.025:1.5)};
N=zeros(length(EDGES{1,1}),length(EDGES{1,2}));
for n = 3:nf
    s=fList(n).name;
    if fList(n).isdir || ~strcmp(s(end-5:end), 'traces')
        continue;
    end
    fid=fopen(s,'r');

    %first line of binary file specifies length of trace
    len=fread(fid,1,'int32');
    disp('The len of the time traces is: ')
    disp(len);

    %number of traces
    Ntraces=fread(fid,1,'int16');
    disp('The number of traces is: ')
    disp(Ntraces/2);

    %raw is a linear array, looks like so
    raw=fread(fid,Ntraces*len,'int16');
    disp('Done reading data.');
    fclose(fid);

    %convert into traces
    index=(1:Ntraces*len);
    Data=zeros(Ntraces,len);
    Data(index)=raw(index);
    
    donor=zeros(Ntraces/2,len);
    acceptor=zeros(Ntraces/2,len);
    EFRET=zeros(Ntraces/2,1);
    STOIC=zeros(Ntraces/2,1);
    
    startFrame=12;
    for i=1:(Ntraces/2)
        donor(i,:)=Data(i*2-1,:);
        acceptor(i,:)=Data(i*2,:);
        %correct for baseline signals
        tempDD=mean(donor(i,(startFrame:startFrame+7)),2)-donorBlank;
        tempAD=mean(acceptor(i,(startFrame:startFrame+7)),2)-acceptorBlank;
        tempDA=mean(donor(i,(2:9)),2)-donor2Blank;
        tempAA=mean(acceptor(i,(2:9)),2)-acceptor2Blank;
        intensities1_10=[intensities1_10;tempDA tempAA];
        %correct for dichroic mirror reflection
        tempDD=tempDD-tempAD*reflection;
        tempAD=tempAD+tempAD*reflection;
        %correct for donor leakage and acceptor direct excitation
        tempAD=tempAD-leakage*tempDD-directEx*tempAA;
        %correct for gamma
        tempDD=detection*tempDD;
        %correct for beta
        tempAA=tempAA/normalEx;
        %calculate corrected EFRET and STOIC
        EFRET(i)=tempAD/(tempDD+tempAD);
        STOIC(i)=(tempDD+tempAD)/(tempDD+tempAD+tempAA);
        %for Chris
        intensityTOTAL11_20=[intensityTOTAL11_20;tempDD+tempAD];
    end
    
    tempTOTAL=[EFRET STOIC];
    tempN = hist3(tempTOTAL,'Edges',EDGES);
    TOTAL=[TOTAL;tempTOTAL];
    N = N+tempN;
end
% for Chris
hist(intensityTOTAL11_20,-100:10:max(intensityTOTAL11_20));
xlabel('Intensity of the 11-20 frames');
ylabel('Counts');

%% Use 2-D intensity plot to determine the reflection factor
hdl0 = figure;
plot(intensities1_10(:,2),intensities1_10(:,1),'o');
hold on;
p = fit(intensities1_10(:,2),intensities1_10(:,1),'poly1');
plot(p);
disp(p.p1);
xlabel('Intensity in Cy5 channel');
ylabel('Intensity in Cy3 channel');
hold off;

%% Use the E-S plot to determine the alpha and delta factors
dumN=max(max(N))-N;
hdl1 = figure;
pc = pcolor(EDGES{1,1},EDGES{1,2},dumN');
set(pc, 'EdgeColor', 'none');
axis square;
zoom on;
colormap([zeros(16,3);hot(8)]);
xlabel('E_{FRET}');
ylabel('Stoichiometry');

warning('if you wish to specify the alpha and delta factors, stop here (ctrl+c)');

input('continue?  ');
%% Reject the donor- and acceptor-only species
scutoff1=input('stoichiometry low cutoff: ','s');
cutoff1=str2num(scutoff1);
if isempty(cutoff1)
    cutoff1=0.3;
end
scutoff2=input('stoichiometry high cutoff: ','s');
cutoff2=str2num(scutoff2);
if isempty(cutoff2)
    cutoff2=0.7;
end

index=logical((TOTAL(:,2)>cutoff1).*(TOTAL(:,2)<cutoff2));
selectedTOTAL=TOTAL(index,:);
selectedTOTAL(:,2)=1./selectedTOTAL(:,2);

%% When multiple FRET subpopulations are present, use linear fitting to
%  determine the gamma and beta factors
hdl2 = figure;
plot(selectedTOTAL(:,1),selectedTOTAL(:,2),'o');
hold on;
xlabel('E_{FRET}');
ylabel('1/Stoichiometry');

pts=[];
[ptsX, ~] = ginput();
for n = 1:length(ptsX)
    index=logical((selectedTOTAL(:,1)>ptsX(n)-0.025).*(selectedTOTAL(:,1)<ptsX(n)+0.025));
    dumTOTAL=selectedTOTAL(index,:);
    pts=[pts;mean(dumTOTAL,1)];
end

plot(pts(:,1),pts(:,2),'rx');
if size(pts,1)>=2
    p = fit(pts(:,1),pts(:,2),'poly1');
    plot(p);
    disp('the gamma and beta factors are: ');
    disp((p.p2-1)/(p.p1+p.p2-1));
    disp(p.p1+p.p2-1);
end
hold off;

warning('if you wish to specify the gamma and beta factors, stop here (ctrl+c)');

input('continue?  ');
%% Plot FRET histograms and save data 
binSize = EDGES{1,1}(2)-EDGES{1,1}(1);
hist(selectedTOTAL(:,1), EDGES{1,1}(1:end-1)+binSize/2.0);
temp=axis;
temp(1)=-0.5;
temp(2)=1.5;
axis(temp);

[counts,centers]=hist(selectedTOTAL(:,1), EDGES{1,1}(1:end-1)+binSize/2.0);
centers=centers';
counts=counts';
xlabel('E_{FRET}');
ylabel('Counts');

save('FRET histogram.dat','counts','-ascii');

% f = fit(centers, counts, 'gauss1');
% hold on;
% plot(f);
% hold off;


