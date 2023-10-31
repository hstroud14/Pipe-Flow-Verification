clear; clc; close all
tag='ver_Erode_SmallU';
%Geometry
Ri_i = 0.005;                     %m
Ro_i = 0.01;                     %m
L = 2;                        %m

%Fluid properties
ubar_i = 6;                   %m/s, inital average velocity
Q = pi()*Ri_i.^2*ubar_i;        %m3/s, remains constant
mu = 1.793e-3;                  %Pa*s, dynamic viscosity of water
Pinlet_i = 8*ubar_i*mu*L/(Ri_i.^2);  %Pa, initial inlet pressure
rho = 997;                      %kg/m3, density of water
Re = rho*ubar_i*2*Ri_i./mu;     %Reynolds number


%Material properties
k = 1e-6;                                           %m/Pa*s, erosion constant
E = 1e6;                       %Pa, youngs modulus
nu = 0.3;                       %Poisson ratio
cp = 0;
ctau = 0;
alph=1e-7;

%Simulation params
tnum = 11;
l = linspace(0,L,322);            %length increment
t = linspace(0,100,tnum);           %s, duration of flow
tinc = abs(t(2)-t(1));
delL = abs(l(2)-l(1));
div=11;
plot10 = tnum./div;
ramp1=linspace(0,1,ceil(tnum./2));
ramp2=linspace(1,1,floor(tnum./2));
ramp3=ones(1,length(ramp1)+length(ramp2));
ramp3(2:length(ramp1)+1)=ramp1;
ramp3(length(ramp1)+2:end-1)=ramp2(1:end-2);
ramp3(1)=0;

%initialize constants
P = zeros(length(l)+1,length(t));
Tau = zeros(length(l),length(t));
SigRi = zeros(length(l),length(t));
SigTi = zeros(length(l),length(t));
SigRo = zeros(length(l),length(t));
SigTo = zeros(length(l),length(t));
SigA = zeros(length(l),length(t));
uri = zeros(length(l),length(t));
uro = zeros(length(l),length(t));
Ri = zeros(length(l),length(t));
Ro = zeros(length(l),length(t));
Riref = zeros(length(l),length(t));
Roref = zeros(length(l),length(t));
ubar = zeros(length(t),1);
delP = zeros(length(l),length(t));
delPL = zeros(length(l),length(t));
Q=Q.*ramp3;

%at time = 0:
Ri(:,1) = Ri_i;
Ro(:,1) = Ro_i;
Riref(:,1) = Ri_i;
Roref(:,1) = Ro_i;

%changing in time
sumtinc=0;
for j=1:length(t)
    %changing in space
    for i = 0:length(l)-1
      z = length(l)-i;

      %solve pressure
      delP(z,j) = 8*mu*Q(j)*delL./(pi().*Ri(z,j).^4);
      delPL(z,j) = 8*mu*Q(j)./(pi().*Ri(z,j).^4);
      P(z,j) = P(z+1,j) + delP(z,j);

      %solve shear stress due to fluid flow
      Tau(z,j) = 4*mu*Q(j)./(pi().*Ri(z,j).^3);

      %solve stresses along length 
%       SigRi(z,j) = P(z,j).*Ri(z,j).^2./(Ro(z,j).^2-Ri(z,j).^2).*(1-Ro(z,j).^2./Ri(z,j).^2);
%       SigRo(z,j) = P(z,j).*Ri(z,j).^2./(Ro(z,j).^2-Ri(z,j).^2).*(1-Ro(z,j).^2./Ro(z,j).^2);
%       SigTi(z,j) = P(z,j).*Ri(z,j).^2./(Ro(z,j).^2-Ri(z,j).^2).*(1+Ro(z,j).^2./Ri(z,j).^2);
%       SigTo(z,j) = P(z,j).*Ri(z,j).^2./(Ro(z,j).^2-Ri(z,j).^2).*(1+Ro(z,j).^2./Ro(z,j).^2);
      SigRi(z,j) = -P(z,j);
      SigRo(z,j) = 0;
      SigTi(z,j) = P(z,j).*(Ri(z,j).^2+Ro(z,j).^2)./(Ro(z,j).^2-Ri(z,j).^2);
      SigTo(z,j) = 2.*P(z,j).*Ri(z,j).^2./(Ro(z,j).^2-Ri(z,j).^2);
      SigA(z,j)  = P(z,j).*Ri(z,j).^2./(Ro(z,j).^2-Ri(z,j).^2);
      
      %solve displacements due to structure
      uri(z,j) = Ri(z,j)./E.*(SigTi(z,j)-nu.*(SigRi(z,j)+SigA(z,j)));
      uro(z,j) = Ro(z,j)./E.*(SigTo(z,j)-nu.*(SigRo(z,j)+SigA(z,j)));
      
      Riref(z,j+1) = Riref(z,j)+k*Tau(z,j)*tinc;
      Roref(z,j+1) = Ro_i;

      Ri(z,j+1) = Riref(z,j+1)+uri(z,j);
      Ro(z,j+1) = Roref(z,j+1)+uro(z,j);
      
      if Ri(z,j+1)>=Ro(z,j+1)
          Ri(z,j+1)=NaN;
          Ro(z,j+1)=NaN;
      end

    end
    sumtinc=sumtinc+tinc;
end

a=linspace(1,0,length(1:plot10:length(t)));
mycolors=jet(div);

Ri_calc(:,:)=Ri(:,1:end-1);

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Radius')
xlabel('Z (mm)')
axis([0,length(l),5e-3,6e-3])

ylabel('R (m)')
hold on
% for k=plot10*6:plot10*7
%     plot(Ri(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:length(t)
    plot(Ri_calc(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
%plot(Ri(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Ri.pdf'))
figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Ref Radius')
xlabel('Z (mm)')
%axis([0,1000,5e-3,5.25e-3])

ylabel('R (m)')
hold on
% for k=plot10*6:plot10*7
%     plot(Riref(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:length(t)
    plot(Riref(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(Riref(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_RiRef.pdf'))
figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Pressure')
xlabel('Z (mm)')
axis([0,length(l),0,1e4])

ylabel('P (Pa)')
hold on
% for k=plot10*6:plot10*7
%     plot(P(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:length(t)
    plot(P(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(P(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Press.pdf'))
figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Shear Stress')
xlabel('Z (mm)')
%axis([0,1000,0,4e-4])
ylabel('Tau (Pa)')
hold on
% for k=plot10*6:plot10*7
%     plot(Tau(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:length(t)
    plot(Tau(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(Tau(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Tau.pdf'))
figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Elastic Displacement')
xlabel('Z (mm)')
axis([0,length(l),0,0.5e-3])

ylabel('Ur (m)')
hold on
% for k=plot10*6:plot10*7
%     plot(uri(:,k),'Color',mycolors(6,:))
% end

count=1;
for k=1:plot10:length(t)
    
    plot(uri(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(uri(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Disp.pdf'))

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Total Displacement')
xlabel('Z (mm)')
axis([0,length(l),0,1e-3])

ylabel('Ur (m)')
hold on
% for k=plot10*6:plot10*7
%     plot(uri(:,k),'Color',mycolors(6,:))
% end

count=1;
for k=1:plot10:length(t)
    
    plot(Ri(:,k)-Ri_i,'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(Ri(:,end)-Ri_i,'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_TotDisp.pdf'))

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Thickness')
xlabel('Z (mm)')
axis([0,length(l),4e-3,5e-3])
thick = Ro(:,1:end-1) - Ri_calc;
hold on
% for k=plot10*6:plot10*7
%     plot(thick(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:length(t)
    plot(thick(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
%plot(thick(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Thick.pdf'))
% fh = findall(0,'Type','Figure');
% set( findall(fh, '-property', 'fontsize'), 'fontsize', 16)
figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Q')
xlabel('Time (s)')
hold on
plot(t,Q,'k','Linewidth',2)
count=1;
for k=1:plot10:length(t)
    plot(t(k),Q(k), 'ok', 'MarkerFaceColor',mycolors(count,:))
    count=count+1;
end
plot(t(end),Q(end), 'ok', 'MarkerFaceColor',mycolors(end,:))
exportgraphics(gca,strcat(tag,'_Q.pdf'))
%%
%clear; clc; %close all
tag='abq_Erode_SmallU';

% 
% fid=fopen('Job-1-INP.txt');
% formatSpec='%d %f %f %f %f %f';
% N=320; %number of elements reported in this file
% 
% C = [];
% linenum=0;
% while ~feof(fid)
%     linenum=linenum+1;
%     if linenum>24 & linenum<347
%         s = textscan(fid,formatSpec,N,'Delimiter','\t');
%         C = [C; s{1,1}, s{1,2}, s{1,3}, s{1,4}, s{1,5}, s{1,6}];
%     end
% end
file=readmatrix('JobErode.csv');
file2=readmatrix('JobErode-Ro.csv');

%N=322;

time=file(:,4);
R_i=file(:,13);
Z=file(:,15);
Ur_i=file(:,16);
S_i=file(:,17);

R_o=file2(:,13);
Ur_o=file2(:,16);
S_o=file2(:,17);


thick=R_o-R_i;

i=1;
N=1;
while Z(i+1)>=Z(i)
    i=i+1;
    N=N+1;
end

numinc=length(Z)/N;
fid = fopen('check.txt') ;
data = textscan(fid,'%f %f',N,'HeaderLines',2) ;
temp = cell2mat(data);
P_abq(:,2)=temp(:,1);
Pvec_abq(:,:,2)=temp(:,:);

for i=2:numinc-1
    clear temp data
    textscan(fid,'%s  %*d',1);
    data=textscan(fid,'%f %f',N-1);
    temp=cell2mat(data);
    P_abq(:,i+1)=temp(:,1);
    Pvec_abq(:,:,i+1)=temp(:,:);
end
fclose(fid);

div=11;
tnum=length(time)/N;
plot10 = ceil(tnum./div);

Ri_abq=zeros(N,floor(length(R_i)/N));

for i=1:length(R_i)/N
    Ri_abq(N:-1:1,i)=R_i(N*(i-1)+1:N*i);
    thick_abq(N:-1:1,i)=thick(N*(i-1)+1:N*i);
end

% mycolors=jet(div);

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Radius')
xlabel('Z (mm)')
axis([0,N,5e-3,6e-3])

ylabel('R (m)')
hold on
% for k=plot10*6:plot10*7
%     plot(Ri(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:tnum
    plot(Ri_abq(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
%plot(N:-1:1,R_i(end-N+1:end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Ri.pdf'))

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Pressure')
xlabel('Z (mm)')
axis([0,N,0,1e4])

ylabel('P (Pa)')
hold on
% for k=plot10*6:plot10*7
%     plot(Ri(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:tnum
    plot(P_abq(end:-1:1,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
%plot(N:-1:1,R_i(end-N+1:end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_P.pdf'))

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Elastic Displacement')
xlabel('Z (mm)')
axis([0,300,0,1e-3])

ylabel('Ur (m)')
hold on
% for k=plot10*6:plot10*7
%     plot(uri(:,k),'Color',mycolors(6,:))
% end

count=1;
for k=1:plot10:tnum
    plot(N:-1:1,Ur_i(N*(k-1)+1:N*k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(N:-1:1,Ur_i(end-N+1:end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Disp.pdf'))

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Thickness')
xlabel('Z (mm)')
axis([0,300,4e-3,5e-3])
hold on
% for k=plot10*6:plot10*7
%     plot(thick(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:tnum
    plot(thick_abq(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(N:-1:1,thick(end-N+1:end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Thick.pdf'))

%% Calculating Errors

Ri_abq=zeros(N,floor(length(R_i)/N));

for i=1:length(R_i)/N
    Ri_abq(N:-1:1,i)=R_i(N*(i-1)+1:N*i);
    thick_abq(N:-1:1,i)=thick(N*(i-1)+1:N*i);
end


if length(l)>N
    sample=N;
    inc=floor(length(l)./N);
    index=1;
    for i=1:N-1
        Ri_ver(i,:)=Ri(index,:);
        thick_ver(i,:)=Ro(index,:)-Ri(index,:);
        index=index+inc;
    end
    Ri_ver(N,:)=Ri(end,:);
    thick_ver(N,:)=Ro(end,:)-Ri(end,:);
elseif length(l)==N
    Ri_ver=Ri_calc;
    thick_ver=Ro(:,1:end-1)-Ri_calc;
    sample=N;
else
    sample=length(l);
    inc=floor(N./length(l));
    index=1;
    for i=1:N-1
        Ri_abq(i,:)=Ri_abq(index,:);
        thick_abq(i,:)=thick_abq(index,:);
        index=index+inc;
    end
    Ri_abq(end,:)=Ri_abq(end,:);
    thick_abq(end,:)=thick_abq(end,:);
end

error_Ri=abs(Ri_abq-Ri_ver)./Ri_ver(1,:);
error_thick=abs(thick_abq-thick_ver)./thick_ver(1,:);
len_plot=linspace(0,L,sample);

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Error R_i')
xlabel('Length')
axis([0,L,0,0.01])
count=1
hold on
for k=1:plot10:length(t)
    
    plot(len_plot,error_Ri(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(len_plot,error_Ri(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_ErrorRi.pdf'))

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Error Thick')
xlabel('Length')
axis([0,L,0,0.01])
count=1
hold on
for k=1:11
    
    plot(len_plot,error_thick(:,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(len_plot,error_thick(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_ErrorThick.pdf'))

%% Pressure Error

% get interpolated analytical pressure

%Stack matlab pressure with length
Pvec=zeros(N,2,numinc);
for i=1:N
    for j=1:numinc
      Pvec(N+1-i,:,j)=[P(i+1,j),l(N+1-i)];
    end
end

%interpolate
% for i=1:length(Pvec_abq)
%     for j=1:11
%         index1=find(Pvec_abq(i,2,j)>=Pvec(:,2,j),1, 'last');
%         index2=find(Pvec_abq(i,2,j)<=Pvec(:,2,j),1, 'first');
%         pslope=(Pvec(index2,1,j)-Pvec(index1,1,j))./(Pvec(index2,2,j)-Pvec(index1,2,j));
%         Pcomp(i,j)=pslope*(Pvec_abq(i,2,j))+Pvec(index1,1,j)-pslope*Pvec(index1,2,j);
%         
%     end
% end
for i=1:numinc
    Pcomp(:,i) = interp1(Pvec(:,2,2),Pvec(:,1,i),Pvec_abq(:,2,2));
end

%plot
figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Pressure')
xlabel('Z (mm)')
axis([0,N,0,1e4])

ylabel('P (Pa)')
hold on
% for k=plot10*6:plot10*7
%     plot(Ri(:,k),'Color',mycolors(6,:))
% end
count=1;
for k=1:plot10:tnum
    plot(Pcomp(end:-1:1,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
%plot(N:-1:1,R_i(end-N+1:end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_Pcomp.pdf'))

%calculate error between Pcomp (which is the analytical version) and Pabq
error_P=abs(Pcomp-P_abq)./Pinlet_i;

error_P(isnan(error_P))=0;

figure
set(gca,'FontSize',16,'LineWidth',1.5)
title('Error P')
xlabel('Length')
axis([0,length(error_P),0,0.01])
count=1
hold on
for k=1:plot10:length(t)
    
    plot(error_P(end:-1:2,k),'Color',mycolors(count,:),'Linewidth',2)
    count=count+1;
end
plot(error_Ri(:,end),'Color',mycolors(end,:),'Linewidth',2)
exportgraphics(gca,strcat(tag,'_ErrorP.pdf'))