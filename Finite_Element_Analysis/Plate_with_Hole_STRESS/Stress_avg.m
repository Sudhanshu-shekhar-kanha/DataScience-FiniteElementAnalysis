function [sig_mean_x,sigb_mean_x,sigtc_mean_x,sig_mean_y,sigb_mean_y,sigtc_mean_y]=Stress_avg(NELEMENTS,C,Bn,RNn,U_uni,U_bi,U_tc)

    
  
    sig_avg=[];
  
    sigb_avg=[];
  
    sigtc_avg=[];
    
    
    for i=1:NELEMENTS
        
        Bea=Bn{i};
        Rea=RNn{i};
        deua=[];
        deba=[];
        detca=[];
            
        for j=1:length(Rea)
            Rexa=Rea(j);
            deunia=U_uni(Rexa);
            deua=[deua,deunia];
            debia=U_bi(Rexa);
            deba=[deba,debia];
            detcia=U_tc(Rexa);
            detca=[detca,detcia];
            
        end
       
        epse=Bea*deua.';
        sige=C*epse;
      
        
        sig_avg=[sig_avg,sige];

        
        epseb=Bea*deba.';
     
        
    
        
        sigeb=C*epseb;
  
        
        sigb_avg=[sigb_avg,sigeb];
        
        epsetc=Bea*detca.';
     
        
     
        
        sigetc=C*epsetc;
 
        
        sigtc_avg=[sigtc_avg,sigetc];
        
        
    end 
      
    sig_mean_x=mean(sig_avg(1,:));
    sigb_mean_x=mean(sigb_avg(1,:));
    sigtc_mean_x=mean(sigtc_avg(1,:));
    
    sig_mean_y=mean(sig_avg(2,:));
    sigb_mean_y=mean(sigb_avg(2,:));
    sigtc_mean_y=mean(sigtc_avg(2,:));
    
    
    
end