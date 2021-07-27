% MAIN CODE
clc
clear
% DATA FROM ABAQUS READ FROM EXCEL SHEETS
CCORD=xlsread('3_Node_Iso_Para_Tri.xlsx',1,'C5:E132');
NCA=xlsread('3_Node_Iso_Para_Tri.xlsx',1,'L5:O212');
NNODES=xlsread('3_Node_Iso_Para_Tri.xlsx',1,'I5');                                                                                                                                                                                              
NELEMENTS=xlsread('3_Node_Iso_Para_Tri.xlsx',1,'I6');
DOFPN=2;

% PLATE IS ASSUMED AS STEEL
t=2;
E=2*10^5 % N/mm^2
nu=0.3
Traction= 10 % N/mm^2
% PLAIN STRESS PROBLEM

GSTIFF=zeros(DOFPN*NNODES,DOFPN*NNODES);

Bn=[];
RNn=[];


for EN=1:NELEMENTS
    [B,AREA,SFx,SFy] = Iso_Para_Tri_ele(EN,NCA,CCORD);
    B1=B;
    B1={B1};
    Bn=[Bn;B1];
    
    [C] = C_plain_stress( E,nu);
    [ K ] = K_ele(B,C,t,AREA);
    
    
    [ RN ] = Row_3Node( EN,NCA );
    Rn=RN;
    Rn={RN};
    RNn=[RNn;Rn];
    
    for i=1:6
        for  j=1:6
            GSTIFF(RN(i),RN(j))=GSTIFF(RN(i),RN(j))+K(i,j);
        end
    end
    
end


% DISPLACEMENT VECTOR 

% UNI AXIAL TENSION 
u_uni=zeros(DOFPN*NNODES,1);
for i=1:DOFPN*NNODES
    u_uni(i)=1;
end

% BI AXIAL TENSION

u_bi=zeros(DOFPN*NNODES,1);
for i=1:DOFPN*NNODES
    u_bi(i)=1;
end

% UNDER TENSION IN ONE DIRECTION AND COMPRESSION IN ANOTHER

u_tc=zeros(DOFPN*NNODES,1);
for i=1:DOFPN*NNODES
    u_tc(i)=1;
end

% TRACTION 

% DATA FROM ABAQUS READ FROM EXCEL SHEETS
ELEMENT_EDGE=xlsread('3_Node_Iso_Para_Tri.xlsx',2,'C6:G15');
NODE_EDGE=xlsread('3_Node_Iso_Para_Tri.xlsx',2,'J6:N16');


NR=NODE_EDGE(:,2);
NL=NODE_EDGE(:,3);
NT=NODE_EDGE(:,4);
NB=NODE_EDGE(:,5);

L=10;
c=(t*L)/2;
Traction=Traction*c;
% UNIAXIAL TRACTION
F_uni=zeros(DOFPN*NNODES,1);
for i=1:length(NR)
    a=2*NR(i)-1;
    F_uni(a)=Traction;
    b=2*NL(i)-1;
    F_uni(b)=-Traction;
end

% BIAXIAL TRACTION
F_bi=zeros(DOFPN*NNODES,1);
for i=1:length(NR)
    a=2*NR(i)-1;
    F_bi(a)=Traction;
    b=2*NL(i)-1;
    F_bi(b)=-Traction;
    c=2*NT(i);
    F_bi(c)=Traction;
    d=2*NB(i);
    F_bi(d)=-Traction;
end

% UNDER TENSION IN ONE DIRECTION AND COMPRESSION IN ANOTHER

F_tc=zeros(DOFPN*NNODES,1);

for i=1:length(NR)
    a=2*NR(i)-1;
    F_tc(a)=Traction;
    b=2*NL(i)-1;
    F_tc(b)=-Traction;
    c=2*NT(i);
    F_tc(c)=-Traction;
    d=2*NB(i);
    F_tc(d)=Traction;
end


% ESSENTIAL BOUNDARY CONDITIONS 

u_uni(2*2)=0;
u_uni(2*9)=0;
u_uni(2*7-1)=0;
u_uni(2*4-1)=0;
u_uni(2*1)=0;
u_uni(2*10)=0;
u_uni(2*6-1)=0;
u_uni(2*5-1)=0;

u_bi(2*2)=0;
u_bi(2*9)=0;
u_bi(2*7-1)=0;
u_bi(2*4-1)=0;
u_bi(2*1)=0;
u_bi(2*10)=0;
u_bi(2*6-1)=0;
u_bi(2*5-1)=0;


u_tc(2*2)=0;
u_tc(2*9)=0;
u_tc(2*7-1)=0;
u_tc(2*4-1)=0;
u_tc(2*1)=0;
u_tc(2*10)=0;
u_tc(2*6-1)=0;
u_tc(2*5-1)=0;


% PENALTY CONDITION for UNIAXIAL TENSION

m=max(max(GSTIFF))*10^50;

Gsu=GSTIFF;
Gsb=GSTIFF;
Gstc=GSTIFF;
for i=1:length(u_uni)
    if u_uni(i)~=1
        Gsu(i,i)=m;
        F_uni(i)=u_uni(i)*m;
        
    end
    
end

U_uni=inv(Gsu)*F_uni;

% PENALTY CONDITION for BIAXIAL TENSION

for i=1:length(u_bi)
    if u_bi(i)~=1
        Gsb(i,i)=m;
        F_bi(i)=u_bi(i)*m;
        
    end
    
end

U_bi=inv(Gsb)*F_bi;




% Reduced Matrix for UNDER TENSION IN ONE DIRECTION AND COMPRESSION IN ANOTHER


for i=1:length(u_tc)
    if u_tc(i)~=1
        Gstc(i,i)=m;
        F_tc(i)=u_tc(i)*m;
        
    end
    
end

U_tc=inv(Gstc)*F_tc;





% NODE & ELEMENT IN VICINITY OF HOLE

ELEMENT_HOLE=xlsread('3_Node_Iso_Para_Tri.xlsx',3,'E6:F21');  
Element_Hole=ELEMENT_HOLE(:,2);
NODE_HOLE=xlsread('3_Node_Iso_Para_Tri.xlsx',3,'J6:K21');


% STRESS CALCULATION FOR EACH ELEMENT IN VICINITY OF HOLE

[eps,sig,epsb,sigb,epstc,sigtc]=Stress_calc(Element_Hole,C,Bn,RNn,U_uni,U_bi,U_tc);
[sig_mean_x,sigb_mean_x,sigtc_mean_x,sig_mean_y,sigb_mean_y,sigtc_mean_y]=Stress_avg(NELEMENTS,C,Bn,RNn,U_uni,U_bi,U_tc);


disp('Elements in vicinity of Hole')
disp(Element_Hole.')


disp('Stresses in vicinity of hole - UNIAXIAL')
disp(sig);

disp('Stresses in vicinity of hole - BIAXIAL')
disp(sigb)

disp('Stresses in vicinity of hole - TENSION IN ONE DIRECTION AND COMPRESSION IN ANOTHER')
disp(sigtc)



% STRESS CONCENTRATION FACTOR

% UNIAXIAL TRACTION
disp('Stress Concentration Factor- Uniaxial')

stconc_uni=max(sig(1,:))/sig_mean_x;
disp(stconc_uni)

disp('Stress Concentration Factor- Biaxial')

stconc_bi=max(sigb(1,:))/sigb_mean_x;
disp(stconc_bi)

disp('Stress Concentration Factor- TENSION IN ONE DIRECTION AND COMPRESSION IN ANOTHER')

stconc_tc=max(sigtc(1,:))/sigtc_mean_x;
disp(stconc_tc)

