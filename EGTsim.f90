program EGT_sim

! ============================================================================!
! DESCRIPTION OF CODE

! Generates numerous ellipsoidal growth tessellation realizations for the
! purpose perfoming a dimension reduction. The first attempt (besides)
! PCA will be locally linear embedding. The distance metric is a local
! kernel average (sweeps across domain) followed by an power spectral density
! computation, and then the 2 norm.

! Since many samples are needed, I will generate each sample on 1 processor
! and this will be done in parallel

  
! END DESCRIPTION OF CODE  
! ============================================================================!
  
use random

implicit none


include 'mpif.h'

logical first
character buffer*150,filname1*150,idstr*150
integer Nx(3), Ngamma, Ntsample,Nsite,Ngrain,Ntot,ieorr, myid, j1, j2, j3,jmcs, &
        numprocs, Nsim,iseed1,irg,SIZE
integer, allocatable :: cellid(:),PUT(:),iseed2(:)
real*8 dx(3),d5bound(5), hard_core, lambdarate
real*8, allocatable :: Gammainv(:),param_mat_CDF(:,:),XEG(:,:)


! INITIALIZE MPI
call MPI_INIT(ieorr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ieorr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ieorr)


Nsim = 10
lambdarate = .0339
Nx = (/64,64,64/)
Ntot = Nx(1)*Nx(2)*Nx(3)
dx = (/.25, .25, .25/)

first = .TRUE.
! initialize rand(0) to be unique wrt myid 
iseed1 = myid*1000+100000
call srand(iseed1)
! initialize rand_poisson and random_normal to be unique wrt myid 
call random_seed(SIZE=j1)
allocate(iseed2(j1),PUT(j1))

do j2 = 1,j1
iseed2(j2) =myid*1050+200310 + j2*35000
enddo 
call random_seed(PUT = iseed2)
deallocate(PUT)

Ngamma = 81*5 ! 9x9 covariance, 5 conditional joint distributions
Ntsample = 1024*5 !chose 1024 elements in matlab file
allocate(param_mat_CDF(Ntsample,18),Gammainv(Ngamma))
if (myid .eq. 0 ) then
open(unit=1,file='/Users/kteferra/Documents/research/projects/ICME/codes/EGT/Target_Gammainv.txt',status='old')
read(1,*) buffer
do j1 = 1,Ngamma
read(1,*) Gammainv(j1)
enddo
close(1)
open(unit=1,file='/Users/kteferra/Documents/research/projects/ICME/codes/EGT/Target_CDFs.txt',status='old')
read(1,*) buffer
do j1 = 1,5
read(1,*) d5bound(j1)
do j2 = 1,1024
read(1,*) (param_mat_CDF(1024*(j1-1)+j2,j3),j3=1,18)
enddo ! j2 = 1,1024
enddo ! j1 = 1,5
close(1)
endif ! if (myid .eq. 0) then

! END READ IN TARGET MARGINAL CDFs AND INVERSE OF SQRT OF COVARIANCE MATRICES
! broadcast to all
do j1 = 1,18
call MPI_BCAST(param_mat_CDF(:,j1),Ntsample,MPI_DOUBLE_PRECISION,&
               0,MPI_COMM_WORLD,ieorr)
enddo ! j1 = 1,18
call MPI_BCAST(d5bound,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ieorr)
call MPI_BCAST(Gammainv,Ngamma,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ieorr)
hard_core = .3


! initialize random module so dont have to recompute certain parameters
! (note the change of first from true to false)
Nsite = random_poisson(real(lambdarate*Ntot*dx(1)*dx(2)*dx(3),4),first)
first = .FALSE.

! start MCS here
do jmcs = 1,Nsim/numprocs
Nsite = random_poisson(real(lambdarate*Ntot*dx(1)*dx(2)*dx(3),4),first)
if (Nsite .lt. 1 ) Nsite = 1

allocate(cellid(Ntot),XEG(Nsite,18))
!! simulate grain parameters to get XEG
irg = 1 ! = 0 for read, = 1 for generate simulation
call Get_XEG1(Nx(1),Nx(2),Nx(3),dx,Nsite, Ngamma,Ntsample,hard_core, &
     XEG,param_mat_CDF,Gammainv,d5bound,irg)



call Get_Morph1(Nx(1),Nx(2),Nx(3),Ntot,Nsite,dx,XEG,cellid,Ngrain)



! write microstructure to VTK. Send all processors to myid=0
!! ---------------------------------------------------------------------------!!
! LATER YOU'LL PROBABLY WANT TO WRITE ALL CELLID'S TO 1 FILE
!! ---------------------------------------------------------------------------!!

if (myid .eq. 0) then
write(idstr,'(i5)') ( (numprocs)*(jmcs-1) + 1 )
idstr = adjustl(idstr)
filname1='Microstructure_test'//idstr(1:len_trim(idstr))//'.vti'
call Write_Microstructure_VTI(Nx(1),Nx(2),Nx(3),dx,cellid,filname1)
do j1 = 1,(numprocs-1)

call MPI_RECV(cellid,Ntot,MPI_INTEGER,j1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ieorr)

write(idstr,'(i5)') ( (numprocs)*(jmcs-1)+j1+1 )
idstr = adjustl(idstr)
filname1='Microstructure_test'//idstr(1:len_trim(idstr))//'.vti'
call Write_Microstructure_VTI(Nx(1),Nx(2),Nx(3),dx,cellid,filname1)


enddo ! do j1 = 1,(numprocs-1)

else ! If (myid .eq. 0) then

call MPI_SEND(cellid,Ntot,MPI_INTEGER, 0, 0, MPI_COMM_WORLD,ieorr)


endif ! if (myid .eq. 0) then
deallocate(cellid,XEG)
enddo ! do jmcs = 1,Nsim/numprocs

call MPI_FINALIZE(ieorr)

end program !EGT_sim
