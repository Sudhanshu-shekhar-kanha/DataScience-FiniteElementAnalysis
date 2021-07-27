function [K8,xxc8,yyc8] = K_8_Node(EN,NCA,CCORD,C,h)
    syms t s
    % SHAPE FUNCTION
    N1=-0.25*(-1+s)*(-1+t)*(1+s+t);
    N2=0.5*(-1+s^2)*(-1+t);
    N3=0.25*(-1+t)*(1-s^2+t+s*t);
    N4=-0.5*(1+s)*(-1+t^2);
    N5=0.25*(1+s)*(1+t)*(-1+s+t);
    N6=-0.5*(-1+s^2)*(1+t);
    N7=0.25*(-1+s)*(1+s-t)*(1+t);
    N8=0.5*(-1+s)*(-1+t^2);
    Ne=[N1 N2 N3 N4 N5 N6 N7 N8];
    
    x=[];
    y=[];
    NN=[];
    
    for j=2:9
        N=NCA(NCA(:,1)==EN,j);
        NN=[NN,N];
        xc=[CCORD(CCORD(:,1)==N,2)];
        x=[x,xc];
        yc=[CCORD(CCORD(:,1)==N,3)];
        y=[y,yc];
    end
    
    xxc8=mean(x);
    yyc8=mean(y);
    
    xs=diff(Ne,s)*x.';
    xt=diff(Ne,t)*x.';
    ys=diff(Ne,s)*y.';
    yt=diff(Ne,t)*y.';
    
    detJ=xs*yt-xt*ys;
    
    A=(1/detJ)*[yt -ys 0 0;
                0 0 -xt xs;
                -xt xs yt -ys];
    
    u=[N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0];
    v=[0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];
    
    us=diff(u,s);
    ut=diff(u,t);
    vs=diff(v,s);
    vt=diff(v,t);
    
    G=[us;ut;vs;vt];
    
    B=A*G;
    
    % 3*3 GAUSS QUADRATURE
    so=0.774597;
    w1=0.308642;
    w2=0.493827;
    w3=0.790123;
    
    B1=subs(B,{s,t},{-so,-so});
    detJ1=subs(detJ,{s,t},{-so,-so});
    
    B2=subs(B,{s,t},{-so,0});
    detJ2=subs(detJ,{s,t},{-so,0});
    
    B3=subs(B,{s,t},{-so,so});
    detJ3=subs(detJ,{s,t},{-so,so});
    
    B4=subs(B,{s,t},{0,-so});
    detJ4=subs(detJ,{s,t},{0,-so});
    
    B5=subs(B,{s,t},{0.00000000000000000001,0.000000000000000000001});
    detJ5=subs(detJ,{s,t},{0.0000000000000000001,0.000000000000000000000001});
    
    B6=subs(B,{s,t},{0,so});
    detJ6=subs(detJ,{s,t},{0,so});
    
    B7=subs(B,{s,t},{so,-so});
    detJ7=subs(detJ,{s,t},{so,-so});
    
    B8=subs(B,{s,t},{so,0});
    detJ8=subs(detJ,{s,t},{so,0});
    
    B9=subs(B,{s,t},{so,so});
    detJ9=subs(detJ,{s,t},{so,so});
    
    k1=w1*(B1.'*C*B1)*detJ1;
    k2=w2*B2.'*C*B2*detJ2;
    k3=w1*B3.'*C*B3*detJ3;
    k4=w2*B4.'*C*B4*detJ4;
    k5=w3*(B5.'*C*B5)*detJ5;
    k6=w2*B6.'*C*B6*detJ6;
    k7=w1*B7.'*C*B7*detJ7;
    k8=w2*B8.'*C*B8*detJ8;
    k9=w1*(B9.'*C*B9)*detJ9;
    
    K8=h*(k1+k2+k3+k4+k5+k6+k7+k8);
    
    
    
    
    
    
end