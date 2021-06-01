!
! Trapezoidal Rule
!
program Trapezoid

implicit none


real*8,allocatable,dimension(:) :: Shom,fg,phi
real*8,allocatable,dimension(:,:) :: Sfin
real*8,dimension(:) :: f(100),x(100),g(100),h(100)
real*8 :: L0,kupper,dk,k,pi,dx,dxs,xs,stdev,Sum1,Sum2,var,fh,fx
integer :: nsub,nfreq,ii,jj,imcs,mesh,ij,maxn,err,jk,kk
!
!write(*,*) "Number of Monte Carlo Simulations, nsim="
!read *,  nsim
!write(*,*) "Function f(x) discretization (no. of subintervals), nsub="
!read *,  nsub
!write(*,*) "Horizontal axis discretization discretization "
!read *,  mesh
!write(*,*) "Specify f(x) length, L="
!read *, L
!write(*,*) "Define upper value of frequency, kupper="
!read *,  kupper
write(*,*) "Specify frequency discretization, nfreq="
read *, nfreq
!
!
!hstep=L/nsub
pi=3.14159265

open(unit=20,file="functionf.txt",status="unknown")
open(unit=30,file="PSDF.out",status="unknown")
open(unit=40,file="Shom.out",status="unknown")
open(unit=50,file="f.out",status="unknown")
open(unit=60,file="mean.out",status="unknown")
open(unit=70,file="results.txt",status="unknown")


allocate(Shom(nfreq+1))
allocate(phi(nfreq+1))


! Loop for the Monte Carlo Simulations

h = 0.
imcs = 1
do
        write(70,*) 'CNT-'
        x=0.
        f=0.
        if (imcs==1)then 
            jj = 0
        else
            jj=1
        end if
        do    
            jj = jj + 1
            read(20,*,iostat=err) x(jj),f(jj)
            if (jj.eq.1) cycle
            if ((x(jj).eq.0..and.f(jj).eq.0.)) then 
                if (imcs.eq.1) then 
                    L0 = x(jj-1)
                    nsub = jj-1
                    dxs = L0 / dble(nsub-1)
                    kupper = pi / dxs
                    dk = kupper /nfreq
                end if
                exit
            end if
        end do
        
        do ij = 1,jj-1
            x(ij) = x(ij) * L0 / x(jj-1)
            if (ij==jj-1) then
               write(70,*) x(ij),0.
            else
               write(70,*) x(ij),f(ij)
            end if
        end do
        
        xs = 0.
        do ii = 1,nsub
            fh = 0.
            do jk = 2,jj-1
                if (xs.le.x(jk)) then
                    fh = f(jk-1) + (f(jk)-f(jk-1))*(xs-x(jk-1))/(x(jk)-x(jk-1))
                    if (abs(xs-x(jk)).le.(1.e-05)) fh = f(jk)
                    exit
                end if
            end do
            h(ii) = h(ii) + fh
            
            if ( (ii==1).or.(ii==nsub))then
                h(ii)=0.
            end if
            
            xs = xs + dxs    
        end do

if (err ==-1) exit
imcs = imcs + 1
end do

h = h/dble(imcs-1)

write (*,*) (h(ii), ii=1,nsub)

rewind(20)
allocate(fg(nsub))
Shom=0.
g = 0.
dk = kupper / nfreq
imcs = 1
do
        !write(70,*) 'CNT-'
        x=0.
        f=0.
        if (imcs==1)then 
            jj = 0
        else
            jj=1
        end if
        do    
            jj = jj + 1
            read(20,*,iostat=err) x(jj),f(jj)  
            if (jj.eq.1) cycle
            if ((x(jj).eq.0..and.f(jj).eq.0.)) exit
        end do
        
        do ij = 1,jj-1
            x(ij) = x(ij) * L0 / x(jj-1)
        end do
        
        xs = 0.
        fg = 0.
        do ii = 1,nsub
            do jk = 2,jj-1
                if (xs.le.x(jk)) then
                    fg(ii) = f(jk-1) + (f(jk)-f(jk-1))*(xs-x(jk-1))/(x(jk)-x(jk-1))
                    !write(70,*) xs,fg(ii)
                    fg(ii) = fg(ii) - h(ii)
                    exit
                end if
            end do
            g(ii) = g(ii) + (fg(ii))**2
            
            if ( (ii==1).or.(ii==nsub))then
                g(ii)=0.
            end if
            
            xs = xs + dxs    
        end do

!*******************************************************************
    k=0.
    do ii = 1,nfreq+1
        Sum1=0.
        Sum2=0.
        xs = 0.
        do ij = 2,nsub
            Sum1=Sum1+(fg(ij-1)*cos(k*(xs-dxs))+fg(ij)*cos(k*xs))*dxs/2.
            Sum2=Sum2+(fg(ij-1)*sin(k*(xs-dxs))+fg(ij)*sin(k*xs))*dxs/2.
            xs = xs + dxs    
        enddo
        Shom(ii) = Shom(ii) + ((1/(2*pi*L0))*(Sum1**2+Sum2**2))
        k = k + dk
    end do
!*****************************************************************
if (err ==-1) exit
imcs = imcs + 1
end do

Shom = Shom / dble(imcs-1)
g = g / dble(imcs-1)
var = 0.
do ii =1,nfreq-1
    Shom(ii+1) =(Shom(ii) + Shom(ii+1) + Shom(ii+2))/3.
    var = var +(Shom(ii+1)+Shom(ii))*dk/2.
end do

var = var + (Shom(nfreq)+Shom(nfreq+1))*dk/2.
var=2.*var !Symmetric part of integral
g = g /var

allocate(Sfin(nsub,nfreq+1))

do jj =1,nfreq+1
    do ii = 1,nsub
        Sfin(ii,jj) =Shom(jj)*g(ii)
    end do
end do 

do kk = 1,10
xs = 0.
call random_number(phi)
do ii = 1,nsub
    fx = 0.
    k = 0.
    do jj = 1,nfreq+1
        fx = fx + sqrt(2.) * sqrt(2. * Sfin(ii,jj) * dk) * cos(k*xs + phi(jj)*2.*pi)
        k = k + dk
    end do
    xs = xs + dxs
    fx = fx + h(ii)
    write (50,999) fx
end do
end do

write (30,888) ((Sfin(ii,jj), ii=1,nsub),jj=1,nfreq+1)	
!write (*,888) ((Sfin(ii,jj), ii=1,nsub),jj=1,nfreq+1)	
write (40,999) (Shom(ii), ii = 1,nfreq+1)
write (60,999) (h(ii), ii=1,nsub)


close(20)
close(30)
close(40)
close(50)

deallocate(Shom)
deallocate(phi)
deallocate(fg)
deallocate(Sfin)

 888  format (47E23.15)	  
 999  format (E23.15)	  

pause
end program