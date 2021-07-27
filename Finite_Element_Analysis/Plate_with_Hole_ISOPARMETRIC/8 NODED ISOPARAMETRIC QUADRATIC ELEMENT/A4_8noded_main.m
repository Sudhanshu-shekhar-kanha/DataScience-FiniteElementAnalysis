% MAIN CODE
clc
clear
% DATA FROM ABAQUS READ FROM EXCEL SHEETS
CCORD=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',1,'C5:E364');
NCA=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',1,'L5:T108');
NNODES=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',1,'I5');                                                                                                                                                                                              
NELEMENTS=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',1,'I6');
DOFPN=2;

% PLATE IS ASSUMED AS STEEL
h=2;
E=2*10^5 % N/mm^2
nu=0.3
Traction= 10 % N/mm^2

% PLAIN STRESS PROBLEM

GSTIFF=zeros(DOFPN*NNODES,DOFPN*NNODES);
[C] = C_plain_stress( E,nu);

XXC8=[];
YYC8=[];
rn8=[];


for EN=1:NELEMENTS
    [K8,xxc8,yyc8] = K_8_Node(EN,NCA,CCORD,C,h);
    XXC8=[XXC8,xxc8];
    YYC8=[YYC8,yyc8];
    
    [RN8] = Row_8(EN,NCA);
    rn8=[rn8;RN8];
    
    for i=1:16
        for  j=1:16
            
            GSTIFF(RN8(i),RN8(j))=GSTIFF(RN8(i),RN8(j))+K8(i,j);
       
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

NODE_EDGE=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',2,'E6:Y9');
ELEMENT_HOLE=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',2,'E6:F21');                                                                                                                                                                                              
NODE_HOLE=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',3,'E5:AR5');

NR=NODE_EDGE(1,:);
NL=NODE_EDGE(2,:);
NT=NODE_EDGE(3,:);
NB=NODE_EDGE(4,:);

L=10;
c=(h*L)/2;
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

m=max(max(GSTIFF))*10^100;

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
    [Nxc8,Nyc8]=N8_element(EN,XXC8,YYC8);
    [d_uni,d_bi,d_tc]=D8_element(EN,U_uni,U_bi,U_tc,rn8);
 
    Uf_uni=Nxc8*d_uni;
    UE_uni=[UE_uni;Uf_uni];
    Vf_uni=Nyc8*d_uni;
    VE_uni=[VE_uni;Vf_uni];
    
    Uf_bi=Nxc8*d_bi;
    UE_bi=[UE_bi;Uf_bi];
    Vf_bi=Nyc8*d_bi;
    VE_bi=[VE_bi;Vf_bi];
    
    Uf_tc=Nxc8*d_tc;
    UE_tc=[UE_tc;Uf_tc];
    Vf_tc=Nyc8*d_tc;
    VE_tc=[VE_tc;Vf_tc];
    
    
end

vpa(UE_uni(1:5));
vpa(VE_uni(1:5));

% NODE & ELEMENT IN VICINITY OF HOLE

ELEMENT_HOLE=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',3,'E10:L10');                                                                                                                                                                                              
NODE_HOLE=xlsread('8Noded_IsoPara_Quadratic_element.xlsx',3,'E5:AR5');


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
NODE_HOLE(1,:).'
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
ELEMENT_HOLE(1,:).'
disp('Displacement distribution over the element[ Near Hole] ')

vpa(Element_uni.')
vpa(Element_bi.')
vpa(Element_tc.')














