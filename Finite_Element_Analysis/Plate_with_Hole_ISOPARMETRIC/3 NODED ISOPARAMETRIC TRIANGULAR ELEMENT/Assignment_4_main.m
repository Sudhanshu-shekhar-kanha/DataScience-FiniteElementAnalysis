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


rn=[];

for EN=1:NELEMENTS
    [B,AREA,SFx,SFy] = Iso_Para_Tri_ele(EN,NCA,CCORD);
 
    
    [C] = C_plain_stress( E,nu);
    [ K ] = K_ele(B,C,t,AREA);
    
    
    [ RN ] = Row_3Node( EN,NCA );
    rn=[rn;RN];
    
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



% SOLUTION FOR EACH ELEMENT

UE_uni=[];
VE_uni=[];
UE_bi=[];
VE_bi=[];
UE_tc=[];
VE_tc=[];
for EN=1:NELEMENTS
    [d_uni,d_bi,d_tc]=D_tri(EN,U_uni,U_bi,U_tc,rn);
 
    Uf_uni=SFx*d_uni;
    UE_uni=[UE_uni;Uf_uni];
    Vf_uni=SFy*d_uni;
    VE_uni=[VE_uni;Vf_uni];
    
    Uf_bi=SFx*d_bi;
    UE_bi=[UE_bi;Uf_bi];
    Vf_bi=SFy*d_bi;
    VE_bi=[VE_bi;Vf_bi];
    
    Uf_tc=SFx*d_tc;
    UE_tc=[UE_tc;Uf_tc];
    Vf_tc=SFy*d_tc;
    VE_tc=[VE_tc;Vf_tc];
    
    
end

vpa(UE_uni(1:5));
vpa(VE_uni(1:5));

% NODE & ELEMENT IN VICINITY OF HOLE

ELEMENT_HOLE=xlsread('3_Node_Iso_Para_Tri.xlsx',3,'E6:F21');                                                                                                                                                                                              
NODE_HOLE=xlsread('3_Node_Iso_Para_Tri.xlsx',3,'J6:K21');

Node_disp_uni=[];
Node_disp_bi=[];
Node_disp_tc=[];

for i=1:length(NODE_HOLE)
    ax=2*NODE_HOLE(i)-1;
    by=2*NODE_HOLE(i);
    axx=U_uni(ax);
    byy=U_uni(by);
    ab=[axx;byy];
    Node_disp_uni=[Node_disp_uni,ab];
    
    abx=U_bi(ax);
    bby=U_bi(by);
    abb=[abx;bby];
    Node_disp_bi=[Node_disp_bi,abb];
    
    atx=U_tc(ax);
    bty=U_tc(by);
    abt=[atx;bty];
    Node_disp_tc=[Node_disp_tc,abt];
    
end

disp('Nodes in vicinity of Hole')
NODE_HOLE(:,2).'
disp('Nodal Displacements of nodes in vicinity of hole')

Node_disp_uni.'
Node_disp_bi.'
Node_disp_tc.'


Element_uni=[];
Element_bi=[];
Element_tc=[];

for i=1:length(ELEMENT_HOLE)
    ex=ELEMENT_HOLE(i);
    exx=UE_uni(ex);
    eyy=VE_uni(ex);
    aeb=[exx;eyy];
    Element_uni=[Element_uni,aeb];
    
    ebx=UE_bi(ex);
    eby=VE_bi(ex);
    ebb=[ebx;eby];
    Element_bi=[Element_bi,ebb];
    
    etx=UE_tc(ex);
    ety=VE_tc(ex);
    ebt=[etx;ety];
    Element_tc=[Element_tc,ebt];
    
end

disp('ELEMENT NEAR HOLE')
ELEMENT_HOLE(:,2).'
disp('Displacement distribution over the element[ Near Hole] ')

vpa(Element_uni.')
vpa(Element_bi.')
vpa(Element_tc.')
