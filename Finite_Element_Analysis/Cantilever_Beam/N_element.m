function [Nxc4,Nxc8,Nyc4,Nyc8]=N_element(EN,L,XXC4,XXC8,YYC4,YYC8)
    a=0.5*(L/20);
    b=0.5*(50/4);
    
    xc4=XXC4(EN);
    yc4=YYC4(EN);
    xc8=XXC8(EN);
    yc8=YYC8(EN);
    
    syms t s x y ;
    
    % SHAPE FUNCTION 4 NODE
    
    N1_4=0.25*(1-s)*(1-t);
    N2_4=0.25*(s+1)*(1-t);
    N3_4=0.25*(s+1)*(1+t);
    N4_4=0.25*(1-s)*(1+t);
    
    Nx4=[N1_4 0 N2_4 0 N3_4 0 N4_4 0];
    Ny4=[0 N1_4 0 N2_4 0 N2_4 0 N4_4];
    
    Nxc4=subs(Nx4,{s,t},{x-xc4,y-yc4});
    Nyc4=subs(Ny4,{s,t},{x-xc4,y-yc4});
    
    % SHAPE FUNCTION
    N1=-0.25*(-1+s)*(-1+t)*(1+s+t);
    N2=0.5*(-1+s^2)*(-1+t);
    N3=0.25*(-1+t)*(1-s^2+t+s*t);
    N4=-0.5*(1+s)*(-1+t^2);
    N5=0.25*(1+s)*(1+t)*(-1+s+t);
    N6=-0.5*(-1+s^2)*(1+t);
    N7=0.25*(-1+s)*(1+s-t)*(1+t);
    N8=0.5*(-1+s)*(-1+t^2);
    
    Nex8=[N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0];
    Ney8=[0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];
    
    Nxc8=subs(Nex8,{s,t},{x-xc8,y-yc8});
    Nyc8=subs(Ney8,{s,t},{x-xc8,y-yc8});
    
    
    
end