%% Use this code to further analyze the selected traces by s_tr.m. A median filter can be applied to EFRET.
%  a 10r-ng-10r illumination pattern should be strictly followed.
%  Use 'f' to analyze traces. A full-set of correction will be applied.
%  Use 'b' to go to the previous molecule, and 'g' to go to a specific molecule.

function s_tr2_E()
close all;
fclose('all');

%read data
pth=input('Directory [default=C:\\User\\tir data\\yyyy\\New Folder]  ');
if isempty(pth)
   	disp('error');
end
cd(pth);

fname=input('index # of filename [default=1]  ');
if isempty(fname)
    fname=1;
end
fname=num2str(fname);

fid = fopen(['hel' fname '.pma'],'r');
if ( fid == -1 )
	disp('file open failed');
end
fileinfo = dir(['hel' fname '.pma']);
film_x = fread( fid, 1, 'uint16' );
film_y = fread( fid, 1, 'uint16' );
headersize = 4;
len = uint32( ( fileinfo.bytes - headersize ) * 1.0 / film_x / film_y );

frame = zeros(film_x,film_y,len,'uint8');
for t=1:len
    frame(:,:,t) = fread(fid,[film_x,film_y],'uint8');
end
fclose(fid);

pth=[pth '\' fname ' selected traces'];
cd(pth);
A=dir;
[nf,~]=size(A);

donor=[];
acceptor=[];
trName={};
trNum=0;

for i=1:nf
    if A(i).isdir == 0
        s=A(i).name;
        if strcmp(s(end-3:end), '.dat')
            disp(s);
            Data=dlmread(s);
            donor=[donor,Data(:,2)];
            acceptor=[acceptor,Data(:,3)];
            trName{end+1}=s;
            trNum=trNum+1;
        end
    end
end

donor=donor';
acceptor=acceptor';

%len=size(Data,1);
timeunit=Data(2,1);
time=Data(:,1);

%% Apply correction factors
%three baselines
donorBlank=0.0; %Ewelina, please change this
acceptorBlank=0.0; %Ewelina, please change this
acceptor2Blank=0.0;
%theta
reflection=0.0;
%alpha
leakage=0.0; %Ewelina, please change this
%delta
directEx=0.0;
%gamma
detection=1.0;
%beta
normalEx=1.0;

%correct for baseline signals
donor(:,11:end-10)=donor(:,11:end-10)-donorBlank;
acceptor(:,11:end-10)=acceptor(:,11:end-10)-acceptorBlank;
tempAA=mean(acceptor(:,2:9),2)-acceptor2Blank;
%correct for dichroic mirror reflection
donor=donor-acceptor*reflection;
acceptor=acceptor+acceptor*reflection;
%correct for donor leakage and acceptor direct excitation
for n = 1:trNum
    acceptor(n,11:end-10)=acceptor(n,11:end-10)-leakage*donor(n,11:end-10)-directEx*tempAA(n);
end
%correct for gamma
donor(:,11:end-10)=detection*donor(:,11:end-10);
%correct for beta
tempAA=tempAA/normalEx;

%% Analyze traces
newfolder = 'HaMMy traces';
mkdir(newfolder);
cd([pth '\' newfolder]);

hdl=figure;
i=0;
countsT=[];
while (trNum-i) > 0
    i = i+1;
    
    %trace window
    figure(hdl);
    ax1=subplot(2,10,[1 9]);
    plot(time,donor(i,:),'g', time,acceptor(i,:),'r', time,donor(i,:)+acceptor(i,:)+400,'c');
    title(['  Molecule ' num2str(i) ' of ' num2str(trNum)]);
    axis tight;
    temp=axis;
    %temp(3)=-temp(4)*0.2;
    temp(4)=temp(4)*1.1;
    if temp(4) < 500
        temp(4)=500;
    end
    axis(temp);
    grid on;
    zoom on;
    
    ax2=subplot(2,10,[11 19]);
    %optional median filter to EFRET
%     fretE = medfilt1(acceptor(i,:),3)./(medfilt1(donor(i,:),3)+medfilt1(acceptor(i,:),3));
    %this is to avoid undefined fretE
    fretE = acceptor(i,:)./(donor(i,:)+acceptor(i,:));
    for m=1:len
        if acceptor(i,m)+donor(i,m)<=0
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

    answer=input('press b-back,g-go,,m-movie,f-FRET histogram  ','s');
    disp(answer);
    
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
    
    if answer=='m'
        pos=regexp(trName{i},'_(\d+)_(\d+)_(\d+)_(\d+).dat','tokens');
        pos=pos{1};
        center_x1=uint16(str2double(pos{1}));
        center_y1=uint16(str2double(pos{2}));
        center_x2=uint16(str2double(pos{3}));
        center_y2=uint16(str2double(pos{4}));
        
        I=display_movie(frame,len,center_x1,center_y1,...
            center_x2,center_y2);
        implay(I);
        input('enter-to continue ','s');
        i=i-1;
    end
    
    if answer=='f'
        %build a FRET histogram using the segment in between two clicks
        [X,~]=ginput(2);
        X=round(X/timeunit);
        output=[time(X(1):X(2)) donor(i,(X(1):X(2)))' acceptor(i,(X(1):X(2)))'];
        save(['HaMMy_' trName{i}],'output','-ascii');
        
        xbins=(-0.4875:0.025:1.4875);
        [counts,centers] = hist(fretE(X(1):X(2)),xbins);
        counts=counts/(X(2)-X(1)+1);
        subplot(2,10,20);
        plot(counts,centers,'bo');
        title('FRET Histogram');
        temp=axis;
        temp(3) = -0.2;
        temp(4) = 1.2;
        axis(temp);
        countsT=[countsT,counts'];
        
        %show the cross-corelation between the donor and acceptor signals
        Ave_donor=mean(donor(i,X(1):X(2)));
        Ave_acceptor=mean(acceptor(i,X(1):X(2)));
        delta_donor=donor(i,X(1):X(2))-Ave_donor;
        delta_acceptor=acceptor(i,X(1):X(2))-Ave_acceptor;
        Ave_donor_cross_acceptor=zeros(1,X(2)-X(1)+1);
        for tau=0:(X(2)-X(1))
            normalize_number=0;
            sum_tau_1=0;
            for t=1:(X(2)-X(1)+1-tau)
                sum_tau_1=sum_tau_1+(delta_donor(t)*delta_acceptor(t+tau));
                normalize_number=normalize_number+1;
            end
            Ave_donor_cross_acceptor(tau+1)=sum_tau_1/normalize_number/(Ave_donor+Ave_acceptor);
        end
        subplot(2,10,10)
        t=0:(X(2)-X(1));
        plot(t,Ave_donor_cross_acceptor,'-b');
        title('Cross Correlation');
        temp=axis;
        temp(1)=-1;
        temp(2)=50;
        temp(4)=2;
        axis(temp);
        grid on;

        input('enter-to continue ','s');
        clf(hdl,'reset');
    end
    
end

cd(pth);
save('FRETResult_tr.dat','countsT','-ascii');

close all;
fclose('all');
end







%% Movie display
function I=display_movie(frame,len,center_x1,center_y1,...
            center_x2,center_y2)

img_size = 17;
I = zeros(img_size,2*img_size,len);
temp1 = zeros(img_size,img_size,len,'uint8');
temp2 = zeros(img_size,img_size,len,'uint8');

for n=1:len
    temp1(:,:,n) = frame(center_x1-(img_size-1)/2:center_x1+(img_size-1)/2,...
        center_y1-(img_size-1)/2:center_y1+(img_size-1)/2,n);
    temp2(:,:,n) = frame(center_x2-(img_size-1)/2:center_x2+(img_size-1)/2,...
        center_y2-(img_size-1)/2:center_y2+(img_size-1)/2,n);
end
temp1 = mat2gray(temp1);
temp2 = mat2gray(temp2);
for n=1:len
    I(:,:,n) = [temp1(:,:,n) temp2(:,:,n)];
end

return;
end
