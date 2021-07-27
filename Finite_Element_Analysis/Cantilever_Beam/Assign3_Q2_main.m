% MAIN CODE
clc
clear
% DATA FROM ABAQUS READ FROM EXCEL SHEETS
CCORD4=xlsread('4_Node_element.xlsx',1,'B5:D109');
NCA4=xlsread('4_Node_element.xlsx',1,'K5:O84');
NNODES4=xlsread('4_Node_element.xlsx',1,'H6');                                                                                                                                                                                              
NELEMENTS4=xlsread('4_Node_element.xlsx',1,'H7');

CCORD8=xlsread('8_Node_element.xlsx',1,'B6:D294');
NCA8=xlsread('8_Node_element.xlsx',1,'K7:S86');
NNODES8=xlsread('8_Node_element.xlsx',1,'G7');                                                                                                                                                                                              
NELEMENTS8=xlsread('8_Node_element.xlsx',1,'G8');

DOFPN=2;
NELEMENTS=NELEMENTS4;
% INPUT DATA : [N,mm Unit]

E=200000;
nu=0.3;
L=300;
A=50*50;
q=-0.01; % N/mm (UDL) downward
P=2.5; % N (Point load upward)
h=2;

% PLAIN STRESS PROBLEM

GSTIFF4=zeros(DOFPN*NNODES4,DOFPN*NNODES4);
GSTIFF8=zeros(DOFPN*NNODES8,DOFPN*NNODES8);
[C] = C_plain_stress( E,nu);
XXC4=[];
YYC4=[];
XXC8=[];
YYC8=[];
rn4=[];
rn8=[];
for EN=1:NELEMENTS
    [K4,xxc4,yyc4] = K_4node(EN,NCA4,CCORD4,C,h);
    XXC4=[XXC4,xxc4];
    YYC4=[YYC4,yyc4];
    
    [K8,xxc8,yyc8] = K_8node(EN,NCA8,CCORD8,C,h);
    XXC8=[XXC8,xxc8];
    YYC8=[YYC8,yyc8];
    
    [RN4,RN8] = ROW( EN,NCA4,NCA8 );
    rn4=[rn4;RN4];
    rn8=[rn8;RN8];
    
    for i=1:8
        for  j=1:8
            
            GSTIFF4(RN4(i),RN4(j))=GSTIFF4(RN4(i),RN4(j))+K4(i,j);
       
        end
    end
    
    
    
     for i=1:16
        for  j=1:16
            
            GSTIFF8(RN8(i),RN8(j))=GSTIFF8(RN8(i),RN8(j))+K8(i,j);
       
        end
    end
    
       
end



% Element having UDL 
E_UDL=[40 39 38 37 36 35 34 33 32 31 2 4 6 8 10 12 14 16 18 20];
N_E_UDL4=[6 40 41 42 43 44 45 46 47 48 2 11 12 13 14 15 16 17 18 19 3];


N_E_UDL8=[6 205 40 200 41 195 42 190 43 185 44 180 45 175 46 170 47 165 48 160 2 111 11 116 12 121 13 126 14 131 15 136 16 141 17 146 18 151 19 156 3];

l=L/20;
co=(h*l)/2;
qy=q*co;



% Forces
%UDL
rq4=zeros(DOFPN*NNODES4,1);
rq8=zeros(DOFPN*NNODES8,1);

for i=1:length(N_E_UDL4)
    
    r1=N_E_UDL4(i);
    rq4(2*r1)=qy;

end

for i=1:length(N_E_UDL8)
    
    e1=N_E_UDL8(i);
    rq8(2*e1)=qy;
  
    
end



% Point Load
F_N=[9];
F4=zeros(DOFPN*NNODES4,1);
F4(2*9)=P;

F8=zeros(DOFPN*NNODES8,1);
F8(2*9)=P;

R4=rq4+F4;
R8=rq8+F8;

% Displacement Vector

D4=zeros(DOFPN*NNODES4,1);
D8=zeros(DOFPN*NNODES8,1);

for i=1:DOFPN*NNODES4
    D4(i)=1;
end

for i=1:DOFPN*NNODES8
    D8(i)=1;
end


% BOUNDARY CONDITIONS 
% Fixed End
N_EBC4=[3 20 4 50 8];
N_EBC8=[3 157 20 154 4 249 50 247 8];

% Roller Support at midspan
N_EBC_M=[7];

for i=1:length(N_EBC4)
    r1=N_EBC4(i);
    
    D4(2*r1-1)=0;
    D4(2*r1)=0;
      
end

D4(2*7)=0;

for i=1:length(N_EBC8)
    
    e1=N_EBC8(i);

    D8(2*e1-1)=0;
    D8(2*e1)=0;

    
end

D8(2*7)=0;


% PENALTY CONDITION

m4=max(max(GSTIFF4))*10^50;
m8=max(max(GSTIFF8))*10^100;

G4=GSTIFF4;
G8=GSTIFF8;
for i=1:length(D4)
    if D4(i)~=1
        G4(i,i)=m4;
        R4(i)=D4(i)*m4;
        
    end
    
end

U4=inv(G4)*R4;

for i=1:length(D8)
    if D8(i)~=1
        G8(i,i)=m8;
        R8(i)=D8(i)*m8;
        
    end
    
end


U8=inv(G8)*R8;


% SOLUTION FOR EACH ELEMENT
Uf4=[];
Vf4=[];
Uf8=[];
Vf8=[];
for EN=1:NELEMENTS
    [Nxc4,Nxc8,Nyc4,Nyc8]=N_element(EN,L,XXC4,XXC8,YYC4,YYC8);
    [d4,d8]=d_element(EN,U4,U8,rn4,rn8);
    Uff4=Nxc4*d4;
    Uf4=[Uf4;Uff4];
    Vff4=Nyc4*d4;
    Vf4=[Vf4;Vff4];
    
    Uff8=Nxc8*d8;
    Uf8=[Uf8;Uff8];
    Vff8=Nyc8*d8;
    Vf8=[Vf8;Vff8];
    
    
end




% PLOTS
% U vs y [taking x= constant=300]

Ex_Last=[80 79 40 49];
Xm_last=300;

cz=[];
for i=1:length(Ex_Last)
    a=Ex_Last(i);
    b=Uf4(a);
    syms x;
    c=subs(b,{x},{Xm_last});
    cz=[cz,c];
    
end

yp=[1,2,3,4];
yl=50/4;
ylin=linspace(0,yl,50);
ccb=[];

for i=1:length(yp)
    cc=cz(i);
    syms y;
    cb=subs(cc,{y},{(i-1)*yl+ylin});
    ccb=[ccb;cb];
end



yf4=[ccb(1,:),ccb(2,:),ccb(3,:),ccb(4,:)];
xf4=[ylin,yl+ylin,(2*yl)+ylin,(3*yl)+ylin];

length(yf4);
length(xf4);


% 8 NODE
Ex_Last=[80 79 40 49];
Xm_last=300;
cz8=[];
for i=1:length(Ex_Last)
    a8=Ex_Last(i);
    b8=Uf8(a8);
    syms x;
    c8=subs(b8,{x},{Xm_last});
    cz8=[cz8,c8];
    
end

yp=[1,2,3,4];
yl=50/4;
ylin=linspace(0,yl,50);
ccb8=[];

for i=1:length(yp)
    cc8=cz8(i);
    syms y;
    cb8=subs(cc8,{y},{(i-1)*yl+ylin});
    ccb8=[ccb8;cb8];
end



yf8=[ccb8(1,:),ccb8(2,:),ccb8(3,:),ccb8(4,:)];
xf8=[ylin,yl+ylin,(2*yl)+ylin,(3*yl)+ylin];

length(yf8);
length(xf8);






% V vs x [taking y= constant=0]

Ey_bottom=[80 78 76 74 72 70 68 66 64 62 41 43 45 47 49 51 53 55 57 59];
y_bottom=0;
ccz=[];
for i=1:length(Ey_bottom)
    aa=Ey_bottom(i);
    bb=Vf4(aa);
    syms y;
    cccc=subs(bb,{y},{y_bottom});
    ccz=[ccz,cccc];
    
end   
    


xl=300/20;
ylinn=linspace(0,xl,50);
cxb=[];

for i=1:20
    c1c=ccz(i);
    syms x;
    cccb=subs(c1c,{x},{((i-1)*xl)+ylinn});
    cxb=[cxb;cccb];
end
    
yyf4=[];
xxf4=[];
for i =1:20
    asd=cxb(i,:);
    yyf4=[yyf4,asd];
    sdf=(i-1)*xl+ylinn;
    xxf4=[xxf4,sdf];
    
end
xxf4=300-xxf4;

yyf4_2=yyf4(1:(length(yyf4)/2));

yyf4_4=yyf4((length(yyf4)/2)+1:length(yyf4));
yyf4=[-yyf4_2,yyf4_4];

% 8 NODE

Ey_bottom=[59 57 55 53 51 49 47 45 43 41 62 64 66 68 70 72 74 76 78 80];
%Ey_bottom=[80 78 76 74 72 70 68 66 64 62 41 43 45 47 49 51 53 55 57 59];
y_bottom=0;
ccz8=[];
for i=1:length(Ey_bottom)
    aa8=Ey_bottom(i);
    bb8=Vf8(aa8);
    syms y;
    cccc8=subs(bb8,{y},{y_bottom});
    ccz8=[ccz8,cccc8];
    
end   
    


xl=300/20;
ylinn=linspace(0,xl,50);
cxb8=[];

for i=1:20
    c1c8=ccz8(i);
    syms x;
    cccb8=subs(c1c8,{x},{((i-1)*xl)+ylinn});
    cxb8=[cxb8;cccb8];
end

yyf8=[];
xxf8=[];
for i =1:20
    asd8=cxb8(i,:);
    yyf8=[yyf8,asd8];
    sdf8=(i-1)*xl+ylinn;
    xxf8=[xxf8,sdf8];
    
end
figure(1)
plot(xf4,yf4)
title('U vs y , 4 NODE , with x =300 as constant')
xlabel('y (mm) - Variation along height')
ylabel('U (mm)- Displacement in X direction')

figure(2)
plot(xf8,yf8)
title('U vs y , 8 NODE with x =300 as constant')
xlabel('y (mm) - Variation along height')
ylabel('U (mm)- Displacement in X direction')

figure(3)
plot(xxf4,yyf4)
title('V vs x , 4 NODE , with y =0 as constant')
xlabel('x (mm) - Variation along Length')
ylabel('V (mm)- Displacement in Y direction')


figure(4)
plot(xxf8,yyf8)
title('V vs x , 8 NODE , with y =0 as constant')
xlabel('x (mm) - Variation along Length')
ylabel('V (mm)- Displacement in Y direction')









