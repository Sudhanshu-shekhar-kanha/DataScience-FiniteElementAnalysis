function [RN4,RN8] = ROW( EN,NCA4,NCA8 )
    
    RN4=[];
    RN8=[];
    
    for i=1:4
        R1=2*(NCA4(NCA4(:,1)==EN,i+1))-1;
        R2=2*(NCA4(NCA4(:,1)==EN,i+1));
        R=[R1,R2];
        RN4=[RN4,R];
    end
    
    for i=1:8
        
    
        r1=2*(NCA8(NCA8(:,1)==EN,i+1))-1;
        r2=2*(NCA8(NCA8(:,1)==EN,i+1));
        r=[r1,r2];
        RN8=[RN8,r];
    end
    
end 