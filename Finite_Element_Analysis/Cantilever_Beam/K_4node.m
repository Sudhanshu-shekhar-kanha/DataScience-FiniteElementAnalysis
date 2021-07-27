function [K4,xxc4,yyc4] = K_4node(EN,NCA4,CCORD4,C,h)
    B=zeros(3,8);
    x=[];
    y=[];
    NN=[];
    syms t s
    
    
    for j=2:5
        N=NCA4(NCA4(:,1)==EN,j);
        NN=[NN,N];
        xc=[CCORD4(CCORD4(:,1)==N,2)];
        x=[x,xc];
        
        yc=[CCORD4(CCORD4(:,1)==N,3)];
        y=[y,yc];
    end
    xxc4=mean(x);
    yyc4=mean(y);
    xs=0.25*[-1+t 1-t 1+t -1-t]*x.';
    xt=0.25*[-1+s -1-s 1+s 1-s]*x.';
    ys=0.25*[-1+t 1-t 1+t -1-t]*y.';
    yt=0.25*[-1+s -1-s 1+s 1-s]*y.';
    
    detJ=xs*yt-xt*ys;
    A=(1/detJ)*[yt -ys 0 0;
                0 0 -xt xs;
                -xt xs yt -ys];
    G=0.25*[-1+t 0 1-t 0 1+t 0 -1-t 0;
            -1+s 0 -1-s 0 1+s 0 1-s 0;
            0 -1+t 0 1-t 0 1+t 0 -1-t;
            0 -1+s 0 -1-s 0 1+s 0 1-s];
    B=A*G;
    
    % 2*2 GAUSS QUADRATURE
    
    B1=subs(B,{s,t},{-0.57735,-0.57735});
    B2=subs(B,{s,t},{-0.57735,0.57735});
    B3=subs(B,{s,t},{0.57735,-0.57735});
    B4=subs(B,{s,t},{0.57735,0.57735});
    
    k1=(B1.'*C*B1)*detJ;
    k2=B2.'*C*B2*detJ;
    k3=B3.'*C*B3*detJ;
    k4=B4.'*C*B4*detJ;
    
    K4=h*(k1+k2+k3+k4);
            
    
end