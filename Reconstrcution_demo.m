clear all;clc;
%Non-contact full-field dynamic strain reconstruction of rotating blades under multi-mode vibration
%Chunyan Ao 2022.6.20
%Please cite if used
%% load demo data
load demo_data.mat;% 
% fs- frequency in simulation
% t- time in simulation
% t_BTT- sampling time of BTT sensor
% y_BTT_no_noise- sampling data by BTT sensor without noise
% y_BTT_noise- sampling data by BTT sensor with 15dB SNR white noise
% region- resonating region of revolution index
% X_strain_S- normal X strain of simulation in resonating region
% Z_strain_S- normal Z strain of simulation in resonating region
% Dis_Mode- blade displacement mode shape 
% Node_arrayd- displacement node serial number (with midside node)
% Strain_Mode- strain mode shape
% Node_array- strain node serial number (with out midside node)
% BTT_node- BTT node serial number in simulation
%% BTT mode shape and rotatinal speed
N1=find(Node_arrayd==BTT_node);
sm=2; % mode number
EO=[4 12]; % engine order
Node_numd=length(Node_arrayd);
Node_num=length(Node_array);
Node_y=[1*Node_numd]+[N1];
BTT_dismode=Dis_Mode(Node_y,1:sm);%BTT node displacement mode shape in Y direction
Time=t_BTT(1,:)';
RPM=zeros(length(Time),1);% rotational speed calculation 
for i=1:1:(length(Time)-1)
    ii=i+1;
    rpm=1/(Time(ii)-Time(i))*60;
    RPM(i)=rpm;
end
RPM(length(Time))=RPM(length(Time)-1);
Freq=[RPM/60*EO(1) RPM/60*EO(2)];
%% Design matrix 
angle_Sensor=[0 5 15 20 30 50 290 325 330 335 350]; % BTT sensor layout
angle=angle_Sensor/360*(2*pi);
for i=1:length(EO)
    E=EO(i); 
    for j=1:length(angle)
    An=angle(j);
    h(j,:)=[sin(E*An) cos(E*An)];    
    end
    H1(:,[2*i-1,2*i])=h;
end
S=ones(length(H1(:,1)),1);
H=[H1]; % fitting without DC offset 
% H=[H1 S];% fitting with DC offset 
%% Response,BTT sampling and rotational speed without noise
y_BTT=y_BTT_no_noise;
figure (1)
[Ax,H1,H2]=plotyy(t,Y,Time,RPM);
hold on
for i=1:11
plot(t_BTT(1,:),y_BTT(i,:),'-')
hold on
end
hold off
legend('Displacement','sensor1','sensor2','sensor3','...')
xlabel('Time (s)');
ylabel('Amplitude (um)');
%% BTT sampling with noise
y_BTT=y_BTT_noise;
figure (2)
for i=1:11
plot(t_BTT(i,:),y_BTT(i,:),'*')
hold on
end
hold off
legend('sensor1','sensor2','sensor3','...')
xlabel('Time (s)');
ylabel('Amplitude (um)');
%% LS model fitting
for i=1:length(y_BTT(1,:))
    P(:,i)=pinv(H)*y_BTT(:,i);
end
for i=1:length(EO)
eval( ['A_',num2str(i),'=P(2*i-1,:);']);
eval( ['B_',num2str(i),'=P(2*i,:);']);
AMP(i,:)=sqrt(P(2*i-1,:).^2+P(2*i,:).^2);
ANG(i,:)=atan( P(2*i,:)./P(2*i-1,:) );
end
%% Fitting results (multi-mode amplitude and phase)
figure (3)
subplot(1,2,1)
plot(Time,AMP);
legend('4EO1M','12EO2M' );
xlabel('Time (s)');
ylabel('Amplitude (um)');
subplot(1,2,2)
plot(Time,ANG);
legend('4EO1M','12EO2M');
xlabel('Time (s)');
ylabel('Phase (rad)');
%% Multi-mode displacement response decoupling
Time_res=Time(region);
det=1/fs;
T_res=Time_res(1):det:Time_res(end);
Freq_res=Freq(region,:);
AMP_res=AMP(:,region);
ANG_res=ANG(:,region);
t_orin=t;
m=0;tt=0;
for i=1:length(Time_res)
    tre=Time_res(i);
    t=T_res(1)+tt;
    k=m+1;
    while t<=tre
        Y_mode1(k)=AMP_res(1,i)*sin( 2*pi*Freq_res(1,1)*t+ANG_res(1,1) );
        Y_mode2(k)=AMP_res(2,i)*sin( 2*pi*Freq_res(1,2)*t+ANG_res(2,1) );
        YY(k)=Y_mode1(k)+Y_mode2(k);
        m=k;
        k=k+1;
        tt=t-T_res(1)+det;
        t=t+det;  
    end
end
t=t_orin;
res_pos=find(t>=Time_res(1)&t<=Time_res(end));
t_res=t(res_pos);% resonating time in simulation
Y_res=Y(res_pos);% resonating displacment in simulation
%% Recovered and simulated displacement comparison
figure (4)
plot(t_res,Y_res);
hold on
plot(T_res,YY);
xlabel('Time (s)');
ylabel('Amplitude (um)');
legend('Simulation','BTT Recovered');
xlim([10.2,10.32])
%%  Full-field displacement-to-strain transform matrix (normal X and Z strains as examples)
X_strainmode=Strain_Mode(1:Node_num,:);
Z_strainmode=Strain_Mode(2*Node_num+1:3*Node_num,:);
T_ds_X=[X_strainmode(:,1)/BTT_dismode(1) X_strainmode(:,2)/BTT_dismode(2)];
T_ds_Z=[Z_strainmode(:,1)/BTT_dismode(1) Z_strainmode(:,2)/BTT_dismode(2)];
%% Full-field strain reconstruction
X_Strain_R=T_ds_X*[Y_mode1;Y_mode2];
Z_Strain_R=T_ds_Z*[Y_mode1;Y_mode2];
%% Strain comparison of each node 
Node=find(Node_array==21599);% find the node number in Node_array (e.g.,6143,21599...)

figure (5) % normal X strain comparison
plot(t_res,X_strain_S(Node,:));
hold on
plot(T_res,X_Strain_R(Node,:),'--');
legend('Simulation','Reconstruction');
xlabel('Time (s)');
ylabel('Strain (\mu\epsilon)');
xlim([10.2,10.32])

figure (6) % normal Z strain comparison
plot(t_res,Z_strain_S(Node,:));
hold on
plot(T_res,Z_Strain_R(Node,:),'--');
xlabel('Time (s)');
ylabel('Strain (\mu\epsilon)');
legend('Simulation','Reconstruction');
xlim([10.2,10.32])
%% Mac values of full-field reconstructed and simulated strain at all the instants
MAC=mac(X_strain_S',X_Strain_R');
MAC_X=diag(MAC)

MAC=mac(Z_strain_S',Z_Strain_R');
MAC_Z=diag(MAC)


