function [RN8] = Row_8(EN,NCA)
    RN8=[];
    
    for i=1:8
        r1=2*(NCA(NCA(:,1)==EN,i+1))-1;
        r2=2*(NCA(NCA(:,1)==EN,i+1));
        r=[r1,r2];
        RN8=[RN8,r];
    end
    
end