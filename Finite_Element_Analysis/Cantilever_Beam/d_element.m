function [d4,d8]=d_element(EN,U4,U8,rn4,rn8)
    a4=rn4(EN,:);
    a8=rn8(EN,:);
    
    d4=zeros(8,1);
    d8=zeros(16,1);
    
    
    for i=1:length(a4)
        d4(i)=U4(a4(i));
        
    end
    
    for i=1:length(a8)
        d8(i)=U8(a8(i));
        
    end
        
    
end