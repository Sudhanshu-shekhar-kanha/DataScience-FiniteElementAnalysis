function [B,AREA,SFx,SFy] = Iso_Para_Tri_ele(EN,NCA,CCORD)
    B=zeros(3,6);
    x=[];
    y=[];
    NN=[];
    SFx=[];
    SFy=[];
    for j=2:4
        N=NCA(NCA(:,1)==EN,j);
        NN=[NN,N];
        xc=[CCORD(CCORD(:,1)==N,2)];
        x=[x,xc];
        yc=[CCORD(CCORD(:,1)==N,3)];
        y=[y,yc];
    end
   
    
    NN;
    alp1=x(2)*y(3)-y(2)*x(3);
    bet1=y(2)-y(3);
    gam1=x(3)-x(2);
    alp2=x(3)*y(1)-y(3)*x(1);
    bet2=y(3)-y(1);
    gam2=x(1)-x(3);
    alp3=x(1)*y(2)-y(1)*x(2);
    bet3=y(1)-y(2);
    gam3=x(2)-x(1);
    B(1,1)=bet1;
    B(1,3)=bet2;
    B(1,5)=bet3;
    B(2,2)=gam1;
    B(2,4)=gam2;
    B(2,6)=gam3;
    B(3,1)=gam1;
    B(3,2)=bet1;
    B(3,3)=gam2;
    B(3,4)=bet2;
    B(3,5)=gam3;
    B(3,6)=bet3;
    AREA=0.5*(alp1+alp2+alp3);
    
    B=(1/(2*AREA))*B;
    
    syms x y
    
    % SHAPE FUNCTIONS
    N1=(1/2*AREA)*(x*bet1+y*gam1+alp1);
    N2=(1/2*AREA)*(x*bet2+y*gam2+alp2);
    N3=(1/2*AREA)*(x*bet3+y*gam3+alp3);
    
    Nx=[N1 0 N2 0 N3 0];
    Ny=[0 N1 0 N2 0 N3];
    
    SFx=[SFx,Nx];
    SFy=[SFy,Ny];
    
    
    
end