    program Lowmeiosis1
!   program for sweeps of diploids with one meiosis per sweep

    CHARACTER*20 FOUT
    CHARACTER*20 FIN
    
    integer :: N
    real :: AN,r,alph,s,h,Pr,R1,R2,P10,Q10

    write(*,*) 'Input file?'
    read(*,*)  FIN
    write (*,*) 'Output file?'
    read  (*,*)  FOUT
    OPEN (1,FILE=FOUT)
    OPEN (2,FILE=FIN)

    read(2,*) s,h
    read(2,*) alph
    read(2,*) N
    read(2,*) ra
    
    AN=2*N
    write(1,*) 'Program for sweep with facultative sex and 1 meiosis per sweep'
    write(1,*) 'Selection coefficient= ',s,' Dominance coefficient= ',h
    write(1,*) 'Frequency of meioses= ', alph
    write(1,*) 'Recombination rate per Mb x alpha=',ra
    write(1,*) 'Population size= ',N
    write(1,*) ''

    write(*,*) 'Program for sweep with facultative sex and 1 meiosis per sweep'
    write(*,*) 'Selection coefficient= ',s,' Dominance coefficient= ',h
    write(*,*) 'Frequency of meioses= ', alph
    write(*,*) 'Recombination rate per Mb x alpha=',ra
    write(*,*) 'Population size= ',N
    write(*,*) ''

    write(*,*) 'Continue?'
    read (*,*) CONT
    if(CONT.eq.1) go to 400

    hs=h*s
    tm=0.5/alph
!   expected time to 1st meiosis

    gamm=AN*hs
!   scaled selection coefficient for sweep of A2A1 to fixation
    Tsc1=2*(log(gamm)+0.5772)/gamm
!   scaled sweep duration
    ts1=Tsc1*AN
!   unscaled sweep duration
    write(1,*) 'Scaled selection coefficient for sweep of A2A1= ',gamm
    write(1,*) 'Scaled sweep duration for A2A1= ',Tsc1
    write(1,*) 'Sweep duration for A2A1= ',ts1
    Pm=alph*ts1
    if(Pm.ge.1) then
    write(1,*) 'Condition for single meiosis violated'
    go to 400
    end if

    write(1,*) ''
    write(1,*) 'Prob. of  meiosis during sweep of A2A1= ',Pm
    
    tm=0.5/alph
!   expected time to meiosis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   No meiosis during sweep of A2A1
    tmc=((1/alph)-0.5*alph*ts1**2)/(1-Pm)
!   expected time to a meiosis after sweep of A2A1, if no meiosis occurs during a sweep
    Tsmc=tmc/AN
!   scaled expected time to a meiosis after sweep of A2A1

    p2=1/(AN*(1-h)*s)
!   frequency of A1A2 among A2 carrying genotypes at start of 2nd stochastic phase of 2nd sweep

    Tsc2=log(2*AN*(1-h)*s)*p2
!   scaled duration of sweep of A2A2 to fixation
    write(1,*) 'Expected scaled time to a meiosis after sweep of A2A1 = ',Tsmc
    write(1,*) 'Expected scaled time to fixation of A2A2 after meiosis post-fixation of A2A1  = ',Tsc2
    write(1,*) ''
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Meiosis occurs during sweep of A2A1

    Tscm2=log(AN*(1-h)*s)+log(AN*hs)-tm*hs
    Tscm2=Tscm2*p2
    write(1,*) 'Expected scaled time to a meiosis during sweep of A2A1= ',Tscm2
!   It is assumed that only the relative frequencies of A2A2 and A2A1 need to be considered
!   after the meiotic event
    
    pidur=(1.0/(AN*alph))+Pm*Tscm2+(1-Pm)*(Tsc1+Tsc2)
    write(1,*)''
    write(1,*)'Net increase in relative diversity during the sweep= ',pidur
    write(1,*)''
    write(1,*) 'Unlinked loci'
    write(1,*)''
    r=0.5
!   free recombination rate
    go to 30

10  read(2,*) z
    write(1,*)''
    write(1,*) 'Distance from sweep (Mb)= ',z
    r=0.5*(1-exp(0.0-2*z/ra))
!   recombination rate between neutral and selected loci
    if(z.le.0.00001) go to 400
    write(1,*)''
30  write(1,*) 'Recombination rate=',r
    Pr=2*r*(1-r)
    write(1,*) 'Post-sweep probability of non-i.b.d.= ',Pr
    pirel=Pr+pidur
    write(1,*) 'Net post-sweep relative diversity= ',pirel

    go to 10

400 end program Lowmeiosis1

    
