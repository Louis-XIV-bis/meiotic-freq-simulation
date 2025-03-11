    program lowmeiosis2
!   program for sweeps of diploids with two meioses per sweep

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
    write(1,*) 'Program for sweep with facultative sex and 2 meioses per sweep'
    write(1,*) 'Selection coefficient= ',s,' Dominance coefficient= ',h
    write(1,*) 'Frequency of meioses= ', alph
    write(1,*) 'Recombination rate per Mb x alpha=',ra
    write(1,*) 'Population size= ',N
    write(1,*) ''

    write(*,*) 'Program for sweep with facultative sex and 2 meioses per sweep'
    write(*,*) 'Selection coefficient= ',s,' Dominance coefficient= ',h
    write(*,*) 'Frequency of meioses= ', alph
    write(*,*) 'Recombination rate per Mb x alpha=',ra
    write(*,*) 'Population size= ',N
    write(*,*) ''

    write(*,*) 'Continue?'
    read (*,*) CONT
    if(CONT.eq.1) go to 400

    hs=h*s
    t1=0.5/alph
!   time to 1st meiosis
    q0=1/(AN*hs)
!   expected frequency of A2A1 at end of 1st stochastic phase of initial sweep.
    qm=(q0/(1-q0))*(1+hs)**t1
    qm=qm/(1+qm)
!   frequency of A2A1 at time of 1st meiosis
    pm=1-qm
!   frequency of A1A1 at time of 1st meiosis

    write(1,*) 't1= ',t1,' qm= ',qm
    
    t2=1.0/alph
!   time between 1st and 2nd meioses

    P10=2*pm*qm/(1-pm**2)
    Q10=(qm**2)/(1-pm**2)
!   conditional frequencies of A2B2/A1B1xA1B1/A1B1 and A2B2/A1B1xA2B2/A1B1

    P1=(1-0.5*qm)/(1-0.25*qm)
    Q1=1-P1
!   conditional frequencies of A2A1 and A2A2 individuals among progeny of 1st meiosis
!   conditioning is on the presence of A2
    write(1,*) '1st meiosis statistics'
    write(1,*) 'P10= ',P10,' Q10= ',Q10
    write(1,*) 'P1= ',P1,' Q1= ',Q1
    write(1,*) ''

    R1=(Q1/P1)*((1+s)/(1+hs))**t2
    Q2=R1/(1+R1)
    P2=1-Q2
!   conditional frequencies of A2A2 and A2A1 individuals at time of 2nd meiosis
    write(1,*) '2nd meiosis statistics'
    write(1,*) 't2= ',t2
    write(1,*) 'P2= ',P2,' Q2= ',Q2
    write(1,*) ''

    P3=(1+hs)/(AN*(1+s))
!   conditional frequency of A2A1 at start of 2nd stochastic phase of 3rd sweep phase
    Q3=1-P3
    t3=(log(Q3*P2)-log(Q2*P3))/((1-h)*s)
!   duration of 3rd sweep phase

    if(t3.le.0) then
    write(1,*) 'A2/A2 is fixed by time of 2nd meiosis'
    write(1,*) 'No meaningful solution for 2 meioses model'
    go to 400
    end if

    write(1,*) '3rd sweep phase statistics'
    write(1,*) 't3= ',t3
    write(1,*) 'P3= ',P3,' Q3= ',Q3
    
    ts1=t1+t2+t3
    Ts=(t1+t2+t3)/AN
!   scaled total sweep duration
    write(1,*) ''
    write(1,*) 'Total sweep duration= ',ts1
    write(1,*) 'Total increase in diversity during sweep= ',Ts

    at=alph*ts1
    if(at.ge.2) then
    write(1,*) 'More than 2 meioses per sweep; no meaningful solution'
    go to 400
    end if
    write(1,*)''
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
    if(z.le.0.00001) go to 400
    write(1,*)''

30  write(1,*) 'Recombination rate=',r
    call meiosis2(P10,Q10,R1,r,Pr)

    pirel=Ts+Pr
!   net post-sweep relative diversity
    write(1,*) 'Post-sweep probability of non-i.b.d.= ',Pr
    write(1,*) 'Net post-sweep relative diversity= ',pirel

    go to 10

400 end program lowmeiosis2

    subroutine meiosis2(P1,Q11,R1,r,Pr)
!   case of 2nd meiosis; calculation of haplotype frequencies accor
    real :: P1,Q1,Q11,R1,X(45),y(3),A(3,45),r
!   P1, Q1 and Q2 are the respective relative frequencies of category 1 A1A2, category 2 A2A1, and category 2 A2A2 individuals
    
    P1=P1
    Q1=Q11
    Q2=Q11*R1
!   R1 is ratio of A2A2 to A2A1 at time of 2nd meiosis
!   P1 is proportional to relative frequency of category 1 A1A2 offspring at time of meiosis 2
!   Q1 is proportional to relative frequency of category 2 A1A2 offspring at time of meiosis 2
!   Q2 is proportional to relative frequency of category 2 A2A2 offspring at time of meiosis 2

    X(1)=0.25*(P1*(1-r))**2
    X(2)=0.5*r*(1-r)*(P1)**2
    X(3)=0.25*(P1*r)**2
    X(4)=0.25*P1*Q2*(1-r)**3
    X(5)=0.25*P1*Q2*(1-r)*r**2
    X(6)=0.5*P1*Q1*(1-r)**3
    X(7)=0.5*P1*Q2*r*(1-r)**2
    X(8)=0.5*P1*Q1*r*(1-r)**2
    X(9)=X(8)
    X(10)=0.5*P1*Q1*(1-r)*r**2
    X(11)=0.25*P1*Q2*r*(1-r)**2
    X(12)=0.25*P1*Q2*r**3
    X(13)=0.5*P1*Q1*r*(1-r)**2
    X(14)=0.5*P1*Q2*(1-r)*r**2
    X(15)=0.5*P1*Q1*(1-r)*r**2
    X(16)=X(15)
    X(17)=0.5*P1*Q1*r**3
    X(18)=0.0625*(Q2**2)*(1-r)**4
    X(19)=0.125*(Q2*r*(1-r))**2
    X(20)=0.25*Q1*Q2*(1-r)**4
    X(21)=0.25*(Q2**2)*r*(1-r)**3
    X(22)=0.25*Q1*Q2*r*(1-r)**3
    X(23)=X(22)
    X(24)=0.25*Q1*Q2*(r*(1-r))**2
    X(25)=0.0625*(Q2**2)*r**4
    X(26)=0.25*Q1*Q2*(r*(1-r))**2
    X(27)=0.25*(Q2**2)*(1-r)*r**3
    X(28)=0.25*Q1*Q2*(1-r)*r**3
    X(29)=X(28)
    X(30)=0.25*Q1*Q2*r**4
    X(31)=0.25*(Q1**2)*(1-r)**4
    X(32)=0.5*Q1*Q2*r*(1-r)**3
    X(33)=0.5*(Q1**2)*r*(1-r)**3
    X(34)=X(33)
    X(35)=0.5*(Q1*r*(1-r))**2
    X(36)=0.25*(Q2*r*(1-r))**2
    X(37)=0.25*Q1*Q2*(r*(1-r))**2
    X(38)=X(37)
    X(39)=0.25*Q1*Q2*(1-r)*r**3
    X(40)=0.25*(Q1*r*(1-r))**2
    X(41)=X(4)
    X(42)=0.5*(Q1**2)*(1-r)*r**3
    X(43)=X(40)
    X(44)=0.5*(Q1**2)*(1-r)*r**3
    X(45)=0.25*(Q1**2)*r**4
!   vector of mating type frequencies for 2nd meiosis

    do 10 i=1,3
    do 20 j=1,45
    A(i,j)=0.0
20  continue
10  continue
!   initializes matrix for calculating haplotype frequencies

    A(1,1)=0.25*(1-r)
    A(3,1)=0.25*r
    A(1,2)=0.125*(1-r)
    A(2,2)=0.125*(1-r)
    A(3,2)=0.25*r
    A(2,3)=A(1,1)
    A(3,3)=A(3,1)
    A(1,4)=0.5*(1-0.5*r)
    A(3,4)=A(3,1)
    A(1,5)=0.25*(1-r)
    A(2,5)=0.25
    A(3,5)=0.25*r
    A(1,6)=0.25*(1-r)
    A(2,6)=0.125*r
    A(3,6)=0.125*r
    A(1,7)=0.125*(3-2*r)
    A(2,7)=0.125
    A(3,7)=0.25*r
    A(1,8)=0.125*(2-r)
    A(3,8)=0.125*r
    A(1,9)=0.125*(1-r)
    A(2,9)=0.125
    A(3,9)=0.125*r
    A(1,10)=0.125
    A(2,10)=0.125*(1-r)
    A(3,10)=0.125*r
    A(1,11)=0.25
    A(2,11)=0.25*(1-r)
    A(3,11)=0.25*r
    A(2,12)=0.25*(2-r)
    A(3,12)=0.25*r
    A(1,13)=0.125*(1-r)
    A(2,13)=0.125
    A(3,13)=0.125*r
    A(1,14)=0.125*(3-2*r)
    A(2,14)=0.125
    A(3,14)=0.25*r
    A(1,15)=0.125
    A(2,15)=0.125*(1-r)
    A(3,15)=0.125*r
    A(2,16)=0.125*(2-r)
    A(3,16)=0.125*r
    A(1,17)=0.125*r
    A(2,17)=0.25*(1-r)
    A(3,17)=A(1,17)
    A(1,18)=1
    A(1,19)=0.5
    A(2,19)=0.5
    A(1,20)=0.25*(2-r)
    A(2,20)=0.25*r
    A(1,21)=0.75
    A(2,21)=0.25
    A(1,22)=0.5
    A(1,23)=0.25
    A(2,23)=0.25
    A(1,24)=0.25*(1+r)
    A(2,24)=0.25*(1-r)
    A(2,25)=1
    A(1,26)=0.25*(1-r)
    A(2,26)=0.25*(1+r)
    A(1,27)=0.25
    A(2,27)=0.75
    A(1,28)=0.25
    A(2,28)=0.25
    A(2,29)=0.5
    A(1,30)=0.25*r
    A(2,30)=0.25*(2-r)
    A(2,31)=0.25*(1-r)
    A(2,31)=0.25*r
    A(1,32)=0.125*(3-2*r)
    A(2,34)=0.125*(1+2*r)
    A(1,33)=0.125*(2-r)
    A(2,33)=0.125*r
    A(1,34)=0.125*(1-r)
    A(2,34)=0.125*(1+r)
    A(1,35)=0.125
    A(2,35)=0.125
    A(1,36)=0.5
    A(2,36)=0.5
    A(1,37)=0.25
    A(2,37)=0.25
    A(1,38)=0.125
    A(2,38)=0.375
    A(1,39)=0.125*(1+2*r)
    A(2,39)=0.125*(2-2*r)
    A(1,40)=0.25
    A(1,41)=0.125
    A(2,41)=0.125
    A(1,42)=0.125*(1+r)
    A(1,42)=0.125*(1-r)
    A(2,43)=0.25
    A(1,44)=0.125*r
    A(2,44)=0.125*(2-r)
    A(1,45)=0.25*r
    A(1,45)=0.25*(1-r)
!   matrix of conditional haplotype frequencies for each mating
    
    yt=0
    do 30 i=1,3
    y(i)=0
    do 40 j=1,45
    y(i)=y(i)+A(i,j)*X(j)
40  continue
    yt=yt+y(i)
30  continue
!   calculates the non-normalized frequencies of the 3 haplotypes
    y(1)=y(1)/yt
    y(2)=y(2)/yt
    write(1,*) ''
    write(1,*) 'A2B2 and A2B1 haplotype freqs= ',y(1),y(2)
    Pr=1-(y(1)**2)-(y(2)**2)
!   probability of non-i.b.d. using normalized haplotype frequencies
    write(1,*) ''
    end subroutine meiosis2
    


