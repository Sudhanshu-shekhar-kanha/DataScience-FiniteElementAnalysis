function [Nxc8,Nyc8]=N8_element(EN,XXC8,YYC8)
   
    xc8=XXC8(EN);
    yc8=YYC8(EN);
    
    syms t s x y ;
    
    % SHAPE FUNCTION 8 NODE
    
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