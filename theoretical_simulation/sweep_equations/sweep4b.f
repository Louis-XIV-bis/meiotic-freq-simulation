    program sweep4b
!   program to calculate sweep effect for mutations with intermediate dominance
!   autosomal inheritance with facultative sex; selection and recombination adjusted for frequency of meioses
!   integrations use Simpsons'rule
    
    CHARACTER*20 FOUT,FIN
    
    real :: a,b,h,e1,FI,rho2,gam,TSC1,TSC2,q1,q2,recp1,recp2,rect1,rect2,cprob1,cprob2,recorr1,recorr2
    real ::  recps1,recps2,phi,delpi10,delpi20,ANP,ra
    integer :: nsimp
    
    write (*,*) 'Output file?'
    read(*,*)  FOUT
    OPEN (1,FILE=FOUT)
    write (*,*) 'Input file?'
    read(*,*)  FIN
    OPEN (2,FILE=FIN)
    
    read (2,*) nsimp
!   Number of points for Simpsons Rule

    write(1,*) 'Program for single sweep with facultative sex (no inbreeding)'
    write(1,*) 'Uses recombination rates calculated for Haldane mapping function & rec. rate per Mb'
    write(1,*) 'Diversity is scaled relative to the neutral value'
    write(1,*) ''
    write(1,*) 'Number of points for Simpsons rule= ',nsimp
    write(*,*) 'Number of points for Simpsons rule= ',nsimp
    read(2,*) s
    write(1,*) ''
    write(1,*) 'Selection coefficient= ',s
    write(*,*) 'Selection coefficient= ',s
    read(2,*) h
    write(1,*) 'Dominance coefficient= ',h
    write(*,*) 'Dominance coefficient= ',h
    write(1,*) ''
    write(*,*) 'Frequency of meioses?'
    read(*,*) alpha
    write(1,*) 'Frequency of meioses=',alpha
    write(*,*) 'Frequency of meioses=',alpha
    write(1,*) ''
    FI=0
    write(1,*) 'Fixation index= ',FI
    
    read(2,*) ra
    write(*,*) 'Recombination rate per Mb x alpha=',ra
    write(1,*) 'Recombination rate per Mb x alpha=',ra
    write(1,*) ''
    read(2,*) NPOP
    write(1,*) 'Population size= ',NPOP
    write(*,*) 'Population size= ',NPOP
    
    ANP=2*NPOP/(1+FI)
    write(1,*) 'Neutral mean coalescent time= ',ANP
    gam0=ANP*s
    write(1,*) 'Scaled selection coefficient= ',gam0
    write(1,*) ''
    
    a=h
    b=alpha*(1-2.0*h)
!   selection parameters

    q1=1.0/(2*a*gam0)
    p1=1-q1
    p2=1.0/(2*(a+b)*gam0)
    q2=1-p2
        
    write(*,*) 'Approximate conditional initial and final frequencies'
    write(*,*) 'q1= ',q1,' q2= ',q2
    write(*,*) 'p1= ',p1,' p2= ',p2
    write(1,*) 'Approximate conditional initial and final frequencies of beneficial mutation'
    write(1,*) 'q1= ',q1,' q2= ',q2
    write(1,*) 'p1= ',p1,' p2= ',p2
    write(1,*) ''
    
    c=a+b*q1
    d=a+b
    e1=a+b*q2
    
    G=(log(q2)+log(c)-log(q1)-log(e1))/a
    G=G+(log(p1)+log(e1)-log(p2)-log(c))/d
    
    A1=((1.0/q1)-1.0)+(b/a)*(log(q1)+log(e1)-log(c)-log(q2))
    A1=((A1/a)+G)/gam0
    PNC=exp(0.0-A1)
!   probability of no coalescence during sweep
        
    t1=(log(q2*c)-log(q1*e1))/a
    t2=(log(p1*e1)-log(p2*c))/d
    Td=(t1+t2)/gam0
!   deterministic fixation time on coalescent timescale
    write(1,*) 'Deterministic time to fixation (coalescent time scale)= ',Td
    write(1,*) ''
    write(1,*) 'Prob. of no coalescent event during sweep (with no rec.)= ',PNC
    write(1,*) ''

    ra=ra/alpha
!   absolute recombination rate per mB

    write(1,*) 'Unlinked loci'
    write(1,*)''
    rho=0.5*ANP*alpha
!   adjusts free recombination rate by frequency of meioses
    go to 30


10  read(2,*) z
    write(1,*) 'Distance from sweep (Mb)= ',z
    rho=0.5*ANP*alpha*(1-exp(0.0-2*z*ra))
    if(z.le.0.00001) go to 20
    write(1,*)''
    write(1,*) 'Scaled recombination rate (random mating autosomal value)= ',rho
    write(1,*) ''
    
30  rfac=rho/(2.0*NPOP)
    S1=2*FI/(1.0+FI)
!   selfing rate

    phi=S1*(2-S-2*(1.0-rfac)*rfac*(2.0-3.0*S1))
    phi=phi/((2.0-S1)*(2.0-(1.0-2.0*rfac*(1.0-rfac))*S1))
!   2-locus IBD probability

    rho1=rho*(1.0-2*FI+phi)/(1.0+FI)
    ros=rho1/gam0
    
    write(1,*) 'Unscaled recombination rate= ',rfac,' phi= ',phi
    write(1,*) 'Scaled recombination rate (true value)= ',rho1,' r/s = ',ros
    
    rho2=rho/(1.0+FI)
!   rho scaled by corrected Ne, but not corrected for inbreeding

   call integ1(cprob1,cprob2,TSC1,TSC2,recp1,recp2,rect1,rect2,recps1,recps2,delpi10,q1,q2,gam0,a,b,rho2,FI,phi,nsimp,ANP)

!   cprob is coalescent probability during sweep
    TS31=TSC1/cprob1
    TS32=TSC1/cprob2
!   TS3 is coalescence time during sweep, conditioned on coalescence
!   recp is probability of recombination during sweep
!   index 1 is the value corrected for selection effect on F; 2 is the uncorrected value

    if(rho1.le.0.001) then
    rect1=0.0
    rect2=0.0
    recp1=0.0
    recp2=0.0
    else
    rect1=rect1/recp1
    rect2=rect2/recp2
    end if
!   rect is mean time of recombination event, conditioned on recombination

    if(rho1.le.0.001) then
    PNR1=1.0
    PNR2=1.0
    else
    call recterm1(a,b,c,FI,gam0,q1,q2,recorr1,recorr2,phi,ANP)
    B1=2*rho2*recorr1
    B2=2*rho2*recorr2
    PNR1=exp(0.0-B1)
    PNR2=exp(0.0-B2)
!   write(*,*) 'recorr1,recorr2= ',recorr1,recorr2
!   write(*,*) 'PNR1,PNR2= ',PNR1,PNR2
    end if
!   PNR is probability of no recombination during sweep

    TS11=PNC*PNR1*Td
    TS12=PNC*PNR2*Td
!   coalescent time with no rec or coal during sweep

    TS21=TSC1+TS11
    TS22=TSC2+TS12
!   net coalescence time for sweep
    
    delpi1=1-recp1-TS21
    delpi2=1-recp2-TS22
!   net reduction in diversity at end of sweep-uncorrected
    delpi11=delpi1-recps1*Td
    delpi12=delpi2-recps2*Td
!   net reduction in diversity at end of sweep with Td correction
    delpi21=delpi1-(recps1*Td+(recp1-recps1)*rect1)
    delpi22=delpi2-(recps2*Td+(recp2-recps2)*rect2)
!   net reduction in diversity at end of sweep with Td and Tw correction

    delpi10=delpi10**2

    write(1,*) ''
!   write(1,*) 'Prob. of no recombination event during sweep (approx, corr. for within-gen seln.)= ',PNR1
!   write(1,*) 'Prob. of no recombination event during sweep (exact, uncorr.)= ',PNR2
!   write(1,*) 'Prob. of coalescence during sweep (approx, corr.)= ',cprob1
!   write(1,*) 'Prob. of coalescence during sweep (exact, uncorr.)= ',cprob2
!   write(1,*) ''
!   write(1,*) 'Mean coalescent time during sweep (approx, corr.)= ',TSC1
!   write(1,*) 'Mean coalescent time during sweep (exact, uncorr.)= ',TSC2
!   write(1,*) 'Mean coalescent time during sweep, conditioned on coal. (approx, corr.) = ',TS31
!   write(1,*) 'Mean coalescent time during sweep, conditioned on coal. (exact, uncorr.) = ',TS32
!   write(1,*) 'Net time to coalescence for sweep (approx, corr.)= ',TS21
    write(1,*) 'Net time to coalescence for sweep (exact, uncorr.)= ',TS22
    write(1,*) ''
!   write(1,*) 'Net prob. of recombination (approx, corr.).= ',recp1
!   write(1,*) 'Net prob. of recombination (exact, uncorr).= ',recp2
!   write(1,*) ''
!   write(1,*) 'recorr1= ',recorr1,' recorr2= ',recorr2
!   write(1,*) 'Net prob. of single recombination event (approx, corr.).= ',recps1
!   write(1,*) 'Net prob. of single recombination event (exact, uncorr.).= ',recps2
!   write(1,*) 'Mean time to recombination during sweep, conditioned on rec. (approx, corr.)= ',rect1
!   write(1,*) 'Mean time to recombination during sweep, conditioned on rec. (exact, uncorr.)= ',rect2
!   write(1,*) ''
!   write(1,*) 'Expected reduction in diversity at end of sweep (approx, corr.)= ',delpi1
!   write(1,*) 'Expected reduction in diversity at end of sweep (exact, uncorr.)= ',delpi2
!   write(1,*) 'Expected red. in diversity at end of sweep (approx, corr.): Td corr. = ',delpi11
    delpi12=1.0-delpi12
    write(1,*) 'Expected relative diversity at end of sweep (exact, uncorr. for inbreed): Td corr. = ',delpi12
    write(1,*) ''
    delpi21=1.0-delpi21
!   write(1,*) 'Expected reduction in diversity at end of sweep (approx, corr.): Td & Tw corr. = ',delpi21
    write(1,*) 'Expected relative diversity at end of sweep (exact, uncorr.): Td & Tw corr. = ',delpi21
    write(1,*) ''
!   write(1,*) 'Expected red. in diversity at end of sweep (direct calculation) = ',delpi10
    
    if(rho.le.0.001) then
    pnreca=1
    else
    pnreca=exp(0.0-rho1*Td)
    end if
    pnreca=1.0-pnreca
    write(1,*) 'Barton approximation for rel. diversity= ',pnreca
    write(1,*) ''
    write(1,*) ''
    go to 10
    
20  end program sweep4b


        subroutine recterm1(a,b,c,FI1,gam,x1,x2,recorr1,recorr2,phi,ANP)
!       calculates recombination factor A2>A1 for given A2 frequency x

        real :: a,b,c,e1,FI1,gam,x1,x2,recorr1,recorr2,phi,ANP
                
        an=1.0/ANP
        e1=a+b*x2
        f1=log(x2)+log(c)-log(x1)-log(e1)
        f2=log(e1)-log(c)
        
        recorr1=0.0
        recorr2=0.0
                
        if(abs(b).le.0.001) then
        recorr1=(1.0-FI1)*f1/gam
        recorr1=(recorr1-an*FI1*(log(x2)-log(x1)+2*(x1-x2)))/a
        else
        recorr1=(((1.0-FI1)/gam)-an*FI1)*f1/a
        recorr1=recorr1-an*f2*(1.0+FI1-2*FI1/b)
        recorr1=recorr1-an*a*(1.0+FI1)*f2/b
        recorr1=recorr1+an*(x2-x1)*(1.0+FI)
        end if
        
        recorr1=recorr1+an*FI1*(log(x2)-log(x1)+2*(x1-x2))
        
        recorr2=(1.0-2*FI1+phi)*f1/(a*gam)
!       uncorrected exact recombination factor
!       write (*,*) 'recorr1= ',recorr1,' recorr2 = ',recorr2
    
        end subroutine recterm1
        
        subroutine recterm2(a,b,FI1,gam,x,q1,recorr1,recorr2,phi,ANP)
!        calulates recombination factor A1>A2 for given A2 frequency x
        
        real :: a,b,FI1,gam,x,q1,recorr1,recorr2,phi,ANP
                
        an=1.0/ANP
        y=1.0-x
        p1=1.0-q1
        c=a+b*x
        f=a+b*q1
        d=a+b
        e2=abs(a-b)
        
        f1=log(p1)+log(c)-log(y)-log(f)
        recorr1=0.0
        recorr2=0.0
        
        if(abs(b).le.0.001) then
        recorr1=(1.0-FI1)*f1/(a*gam)
        recorr1=recorr1-an*FI*(2*(q-q1)+log(y)-log(p1))/a
        else
            if(e2.le.0.001) then
            recorr1=(1.0-FI1)*f1/(d*gam)
            recorr1=recorr1-an*(1+FI1)*((q-q1)+(a/b)*(log(c)-log(f)))
            recorr1=recorr1-an*FI1*(2/b)*(log(c)-log(f))
            recorr1=recorr1-an*FI1*(0.3333*(y**3-p1**3)/b)
        else
        recorr1=(1.0-FI1)*f1/(d*gam)
        recorr1=recorr1-an*(1+FI1)*((q-q1)+(a/b)*(log(c)-log(f)))
        recorr1=recorr1-an*FI1*(2/b)*(log(c)-log(f))
        recorr1=recorr1-an*FI1*(log(y)+log(f)-log(p1)-log(c))/(a-b)
            end if
        end if
        recorr1=recorr1+an*FI1*(2*(q-q1)+log(p)-log(p1))


!       recorr2=(1.0-2*FI+phi)*f1/(d*gam)
!       uncorrected exact recombination factor

    end subroutine recterm2
    
        subroutine integ1(ai11,ai12,ai21,ai22,ai31,ai32,ai41,ai42,ai51,ai52,ai61,q1,q2,gam,a,b,rho,FI,phi,nsimp,ANP)
!       uses Simpson's rule
!       integrals for coal. prob., coal. time, rec. prob. and rec. time

        real :: gam,ai11,ai12,ai21,ai22,ai31,ai32,ai41,ai42,ai51,ai52,ai61,q1,q2,rho,FI,phi
        real :: a,b,c,e1,f,recorr11,recorr21,recorr12,recorr22,ANP
        integer :: nsimp
        
        p1=1.0-q1
        p2=1.0-q2
        
        d=a+b
        e1=a+b*q2
        f=a+b*q1
        
        rhoc=rho*(1.0-2*FI+phi)
!       rho corrected for inbreeding
                
        ai11=0.0
        ai12=0.0
        ai21=0.0
        ai22=0.0
        ai31=0.0
        ai32=0.0
        ai41=0.0
        ai42=0.0
        ai51=0.0
        ai52=0.0
        
        ans=nsimp
        DX=(q2-q1)/ans
        j=0
        F11=0.0
        F12=0.0
        F31=0.0
        F32=0.0
        N1=nsimp+1
        
        do 50 i=1,N1
        
        x=q1+(i-1)*DX
        y=1.0-x
!       x is frequency of favourable allele
        c=a+b*x
                
        delxi=1.0/(x*c)
!       inverse of selection equation, without gamma or y factors
        
        if(i.eq.N1) then
        G=0.0
        PNC=1.0
        PNR=1.0
        go to 40
        else
        G=(log(q2)+log(c)-log(x)-log(e1))/a
        G=G+(log(y)+log(e1)-log(p2)-log(c))/d
        G=G/gam
!       backwards time to freq. x, on coalescent timescale

        G1=(x*(a+b*q1))/(q1*c)
        G2=(p1*c)/(y*(a+b*q1))
!       functions for direct calculation integral

        A1=((1.0/x)-(1.0/q2))+(b/a)*(log(x)+log(e1)-log(c)-log(q2))
        A1=A1/a
        A1=(A1/gam)+G
        PNC=exp(0.0-A1)
!       probability of no coalescence by frequency x
        end if

        if(rho.le.0.001) then
        PNR11=1.0
        PNR12=1.0
        PNR31=1.0
        PNR32=1.0
        else
        call recterm1(a,b,c,FI,gam,x,q2,recorr11,recorr21,phi,ANP)
        B11=2*rho*recorr11
        B12=2*rho*recorr21
!       index 1 is the value corrected for selection effect on F; 2 is the uncorrected value
        
        PNR11=exp(0.0-B11)
        PNR12=exp(0.0-B12)
!       probability of no recombination by frequency x

        call recterm2(a,b,FI,gam,x,q1,recorr12,recorr22,phi,ANP)
        B21=rho*recorr12
        B22=rho*recorr22
!       factor for prob. of no extra A1>A2 recombination

        call recterm1(a,b,f,FI,gam,q1,x,recorr11,recorr21,phi,ANP)
        B31=rho*recorr11
        B32=rho*recorr21
!       factor for prob. of no extra A2>A1 recombination

        PNR31=exp(0.0-B21-B31)
        PNR32=exp(0.0-B22-B32)
        end if
!       probability of no additional recombination by frequency q1

40      F11=delxi*PNC*PNR11/(x*y*gam)
        F12=delxi*PNC*PNR12/(x*y*gam)
!       integrands for prob. of coalescence during sweep
        
        recfac2=1.0-2*FI+phi
        recfac1=1.0-FI-gam*(x*y*(1+FI)*b+FI*(1.0-2*x))/ANP
        recfac1=recfac1+gam*FI*((a+b*x)*(1.0-2*x))/ANP
!       corrected and uncorrected recombination terms

        F31=2*rho*recfac1*delxi*PNC*PNR11/gam
        F32=2*rho*recfac2*delxi*PNC*PNR12/gam
!       integrands for prob. of  recombination during sweep

        F51=F31*PNR31
        F52=F32*PNR32
!       integrands for prob. of single recombination during sweep
            
        G11=G1**(0.0-rhoc/(a*gam))
        G12=G2**(0.0-rhoc/(d*gam))
        
        F61=G11*G12
!       integrand for direct calculation of reduction in diversity
        
        if(i.eq.N1) then
        F21=0.0
        F22=0.0
        F41=0.0
        F42=0.0
        else
        F21=F11*G
        F22=F12*G
!       integrands for mean time to coalescence, on coalescent timescale
        F41=F31*G
        F42=F32*G
!       integrand for mean time to recombination, on coalescent timescale
        end if
        
        if(i.eq.1) then
        ai11=F11
        ai12=F12
        ai21=F21
        ai22=F22
        ai31=F31
        ai32=F32
        ai41=F41
        ai42=F42
        ai51=F51
        ai52=F52
        ai61=F61
        go to 50
        end if
            if(j.eq.0) then
            ai11=ai11+4*F11
            ai12=ai12+4*F12
            ai21=ai21+4*F21
            ai22=ai22+4*F22
            ai31=ai31+4*F31
            ai32=ai32+4*F32
            ai41=ai41+4*F41
            ai42=ai42+4*F42
            ai51=ai51+4*F51
            ai52=ai52+4*F52
            ai61=ai61+4*F61
            
            j=1
            else
            ai11=ai11+2*F11
            ai12=ai12+2*F12
            ai21=ai21+2*F21
            ai22=ai22+2*F22
            ai31=ai31+2*F31
            ai32=ai32+2*F32
            ai41=ai41+2*F41
            ai42=ai42+2*F42
            ai51=ai51+2*F51
            ai52=ai52+2*F52
            ai61=ai61+2*F61
            j=0
            end if
50        continue
        
        ai11=DX*ai11/3.0
        ai12=DX*ai12/3.0
        ai21=DX*ai21/3.0
        ai22=DX*ai22/3.0
        ai31=DX*ai31/3.0
        ai32=DX*ai32/3.0
        ai41=DX*ai41/3.0
        ai42=DX*ai42/3.0
        ai51=DX*ai51/3.0
        ai52=DX*ai52/3.0
        ai61=DX*ai61/3.0
        ai62=DX*ai62/3.0

    end subroutine integ1

