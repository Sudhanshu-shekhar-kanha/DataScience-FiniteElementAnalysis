function [ C ] = C_plain_stress( E,nu)
    C=zeros(3,3);
    C(1,1)=1;
    C(1,2)=nu;
    
    C(2,1)=nu;
    C(2,2)=1;
    
    C(3,3)=(1-nu)/2;
    C=(E/(1-nu^2))*C;
end