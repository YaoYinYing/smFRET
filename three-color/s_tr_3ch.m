%% Use this code to go through all traces and select the good ones.
%  This is the three-color version.
%  A 10r-ng(ir)-10r illumination pattern during data acquisition should be strictly followed.
%  Use 's' to save individual RAW traces as .dat files.
%  Use 'b' to go to the previous molecule, and 'g' to go to a specific molecule.
%  Use 'k' to click on the Cy3/Cy5 baselines, and 'l' to on the Cy3/Cy5 intensities.

function s_tr_3ch()
close all;
fclose('all');

%% Read data
pth=input('directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
    disp('error');
end
cd(pth);
save_file=pth;

fname=input('index # of filename [default=1]  ');
if isempty(fname)
    fname=1;
end
fname=num2str(fname);
disp(['hel' fname '.traces']);

timeunit=input('time unit [default=0.1 sec]  ');
if isempty(timeunit)
    timeunit=0.1;
end

%select the folder to which the files need to be saved
newfolder = [fname ' selected traces'];
mkdir(newfolder);

fid=fopen(['hel' fname '.traces'],'r');

%first line of binary file specifies length of trace
len=fread(fid,1,'int32');
disp('The len of the time traces is: ')
disp(len);

%number of traces
Ntraces=fread(fid,1,'int16');
disp('The number of traces is: ')
disp(Ntraces/3);

%raw is a linear array, looks like it
raw=fread(fid,Ntraces*len,'int16');
disp('Done reading data.');
fclose(fid);

%read .pks files
rawPeakPos = dlmread(['hel' fname '.pks']);
rawPeakPos = rawPeakPos(:,2:3);

%select which frames are used in IDL for spot picking
selection=input('spots selected using frame 1-10-[1], frame 11-20-[2], or frame 1-20-[12]  ');
if selection == 1
    startFrame = 12;
else
    startFrame = 2;
end

%set alpha correction
LEAKAGE=0.0;

%convert raw into traces
index=(1:Ntraces*len);
Data=zeros(Ntraces,len);
Data(index)=raw(index);

rdonor=zeros(Ntraces/3,len);
racceptor=zeros(Ntraces/3,len);
rpresenter=zeros(Ntraces/3,len);
peakPos=zeros(Ntraces/3,6);
total=zeros(Ntraces/3,1);
for i=1:(Ntraces/3)
    rdonor(i,:)=Data(i*3-2,:);
    racceptor(i,:)=Data(i*3-1,:);
    rpresenter(i,:)=Data(i*3,:);
    peakPos(i,1:2)=rawPeakPos(i*3-2,:);
    peakPos(i,3:4)=rawPeakPos(i*3-1,:);
    peakPos(i,5:6)=rawPeakPos(i*3,:);
    tempD=sum(rdonor(i,(startFrame:startFrame+7)),2);
    tempA=sum(racceptor(i,(startFrame:startFrame+7)),2);
    total(i)=(tempD+tempA)/8.; 
end

%% Rough rejection of singly labelled species
disp('use intensity to remove singly labelled species');
figure;
hist(total,-100:10:2000);
grid on;
zoom on;
xlabel('Intensity of the selected frames');
ylabel('Counts');

fcutoff1=input('intensity low cutoff: ','s');
cutoff1=str2num(fcutoff1);
if isempty(cutoff1)
    cutoff1=300;
end
fcutoff2=input('intensity high cutoff: ','s');
cutoff2=str2num(fcutoff2);
if isempty(cutoff2)
    cutoff2=1000;
end

index=logical((total < cutoff1)+(total > cutoff2));
rdonor(index,:)=[];
racceptor(index,:)=[];
rpresenter(index,:)=[];
peakPos(index,:)=[];
peakPos=round(peakPos);

N_mol=size(rdonor,1);
disp(['there are ' num2str(N_mol) ' traces']);

%% View remaining baseline intensities after photobleaching
% for i=1:N_mol
%     tD(i)=sum(rdonor(i,(end-18:end-11)),2)/8.;
%     tA(i)=sum(racceptor(i,(end-18:end-11)),2)/8.;
% end
% 
% figure;
% plot(tD,tA,'x');
% axis square;
% temp=axis;
% temp(1)=0;
% temp(3)=0;
% temp(2)=max(temp(2),temp(4));
% temp(4)=temp(2);
% axis(temp);
% grid on;
% zoom on;
% [donor_blank,acceptor_blank]=ginput(1);
% disp(donor_blank);
% disp(acceptor_blank);
donor_blank=0;
acceptor_blank=0;
presenter_blank=0;

donor=rdonor-donor_blank;
acceptor=racceptor-acceptor_blank;
presenter=rpresenter-presenter_blank;
CorB=[];
CorL=[];

%% View traces
hdl=figure;
i=0;
time=(0:(len-1))*timeunit;
while (N_mol-i) > 0
    i = i+1 ;
    
    %trace window
    figure(hdl);
    ax1=subplot(3,1,1);
    plot(time,donor(i,:),'g', time,acceptor(i,:)-LEAKAGE*donor(i,:),'r', time,donor(i,:)+acceptor(i,:)-LEAKAGE*donor(i,:)+400,'c');
    title(['  Molecule ' num2str(i) ' of ' num2str(N_mol)]);
    axis tight;
    temp=axis;
    %temp(3)=-temp(4)*0.2;
    temp(4)=temp(4)*1.1;
    if temp(4) < 0.5*cutoff2
        temp(4)=0.5*cutoff2;
    elseif temp(4) > 2*cutoff2
        temp(4)=2*cutoff2;
    end
    axis(temp);
    grid on;
    zoom on;
    
    ax2=subplot(3,1,2);
    fretE=(acceptor(i,:)-LEAKAGE*donor(i,:))./(donor(i,:)+acceptor(i,:)-LEAKAGE*donor(i,:));
    %this is to avoid undefined fretE
    for m=1:len
        if donor(i,:)+acceptor(i,:)-LEAKAGE*donor(i,:)<=0
            fretE(m)=-0.2;
        end
    end
    fretE(fretE>1.2)=1.2;
    fretE(fretE<-0.2)=-0.2;
    plot(time,fretE,'b');
    axis tight;
    temp=axis;
    temp(3)=-0.2;
    temp(4)=1.2;
    axis(temp);
    linkaxes([ax1,ax2],'x');
    grid on;
    zoom on;
    
    ax3=subplot(3,1,3);
    plot(time,presenter(i,:),'k');
    axis tight;
    temp=axis;
    %temp(3)=-temp(4)*0.2;
    temp(4)=temp(4)*1.1;
    if temp(4) < 300
        temp(4)=300;
    end
    axis(temp);
    linkaxes([ax1,ax3],'x');
    grid on;
    zoom on;
    
    answer=input('press b-back,g-go,k-baseline,l-leakage,s-save  ','s');
    
    if answer=='b'
        i=i-2;
        if i<0
            i=0;
        end
    end

    if answer=='g'
        mol= input('which molecule do you choose:  ');
        i= mol-1;
    end
    
    if answer=='k'
        [~,Y]=ginput(2);
        correction=zeros(2,1);
        correction(1)=Y(1);
        correction(2)=Y(2);
        CorB=[CorB,correction];
        i=i-1;
    end
    
    if answer=='l'
        [~,Y]=ginput(2);
        correction=zeros(2,1);
        correction(1)=Y(1);
        correction(2)=Y(2);
        CorL=[CorL,correction];
        i=i-1;
    end
    
    %to save individual traces
    if answer=='s'
        output=[time' rdonor(i,:)' racceptor(i,:)' rpresenter(i,:)'];
        save([save_file '\' newfolder  '\hel' fname '_tr' num2str(i) ...
            '_' num2str(peakPos(i,1)) '_' num2str(peakPos(i,2))...
            '_' num2str(peakPos(i,3)) '_' num2str(peakPos(i,4))...
            '_' num2str(peakPos(i,5)) '_' num2str(peakPos(i,6))...
            '.dat'],'output','-ascii');
    end
end

disp('the baseline intensities are: ');
disp(mean(CorB,2));
blank=mean(CorB,2);
CorL(1,:)=CorL(1,:)-blank(1);
CorL(2,:)=CorL(2,:)-blank(2);
LEAKAGE=mean(CorL,2);
disp('the leakage is: ');
disp(LEAKAGE(2)/LEAKAGE(1));

close all;
fclose('all');


