function [ RN ] = Row_3Node( EN,NCA )
    
    RN=[];
    for i=1:3
        R1=2*(NCA(NCA(:,1)==EN,i+1))-1;
        R2=2*(NCA(NCA(:,1)==EN,i+1));
        R=[R1,R2];
        RN=[RN,R];
    end
end 