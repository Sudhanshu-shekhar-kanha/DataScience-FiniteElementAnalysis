function [d_uni,d_bi,d_tc]=D_tri(EN,U_uni,U_bi,U_tc,rn)
    
    a8=rn(EN,:);
    d_uni=zeros(6,1);
    d_bi=zeros(6,1);
    d_tc=zeros(6,1);

    
    for i=1:length(a8)
        d_uni(i)=U_uni(a8(i));
        d_bi(i)=U_bi(a8(i));
        d_tc(i)=U_tc(a8(i));
        
    end
        
    
end