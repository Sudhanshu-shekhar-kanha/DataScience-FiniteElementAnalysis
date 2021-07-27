function [eps,sig,epsb,sigb,epstc,sigtc,sig_mean_x]=Stress_calc(Element_Hole,C,Bn,RNn,U_uni,U_bi,U_tc)
    eps=[];
    sig=[];
    epsb=[];
    sigb=[];
    epstc=[];
    sigtc=[];
    
 
  

    
    for i=1:length(Element_Hole)
        a=Element_Hole(i);
        Be=Bn{a};
        Re=RNn{a};
        deu=[];
        deb=[];
        detc=[];
            
        for j=1:length(Re)
            Rex=Re(j);
            deuni=U_uni(Rex);
            deu=[deu,deuni];
            debi=U_bi(Rex);
            deb=[deb,debi];
            detci=U_tc(Rex);
            detc=[detc,detci];
            
        end
        
        
        
        epse=Be*deu.';
        epse1=epse;
        
        eps=[eps,epse1];
        
        sige=C*epse;
        sige1=sige;
        
        sig=[sig,sige1];

        
        epseb=Be*deb.';
        epse1b=epseb;
        
        epsb=[epsb,epse1b];
        
        sigeb=C*epseb;
        sige1b=sigeb;
        
        sigb=[sigb,sige1b];
        
        epsetc=Be*detc.';
        epse1tc=epsetc;
        
        epstc=[epstc,epse1tc];
        
        sigetc=C*epsetc;
        sige1tc=sigetc;
        
        sigtc=[sigtc,sige1tc];
        
    end 
    
    
    
    
        
        
        
        
        
        
    
    
end