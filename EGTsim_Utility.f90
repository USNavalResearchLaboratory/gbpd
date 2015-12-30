! ============================================================================!
subroutine Get_Morph1(N1,N2,N3,Ntot,Nsite,dx,XEG,cellid,Ngrain)

! this subroutine generates an EGT microstructure and returns 
! vector cellid, and Gsize (which is equiv diameter) 

implicit none
integer N1,N2,N3,Ntot,j1,j2,j3, Ngrain,gclen,z, blen,L1,L2,Lv1, &
        Nslip,Nsite,c,j,cellid(Ntot),Nltmp
integer, allocatable :: gcompete(:),gc(:), potneigh(:), &
 itmp1(:),itmp2(:), badgrain(:,:),voxelid(:),baccounted(:), gneighbor(:,:), &
 gneighborL(:),itmp3(:)
real*8 dx(3), XEG(Nsite,18), node(3)
real*8, allocatable :: x1(:),x2(:),x3(:),xp1(:), xp2(:),xp3(:),th(:),ph(:), &
                       vang(:),d(:),t1(:)

allocate(gcompete(Nsite) )
gcompete = (/ (j1, j1=1,Nsite) /)

cellid(:) = 0
allocate(x1(Nsite),x2(Nsite),x3(Nsite),xp1(Nsite),xp2(Nsite),xp3(Nsite), & 
         th(Nsite),ph(Nsite),vang(Nsite), d(Nsite),t1(Nsite) )
do j1 = 1,N1
 do j2 = 1,N2
  do j3 = 1,N3
node = (/ ((real(j1)-.5)*dx(1)), ((real(j2)-.5)*dx(2)), ((real(j3)-.5)*dx(3)) /)
x1 = node(1) - XEG(gcompete,1)
x2 = node(2) - XEG(gcompete,2)
x3 = node(3) - XEG(gcompete,3)
d = sqrt( x1**2 + x2**2 + x3**2 )
xp1 = XEG(gcompete,7)*x1 + XEG(gcompete,8)*x2 + XEG(gcompete,9)*x3
xp2 = XEG(gcompete,10)*x1 + XEG(gcompete,11)*x2 + XEG(gcompete,12)*x3
xp3 = XEG(gcompete,13)*x1 + XEG(gcompete,14)*x2 + XEG(gcompete,15)*x3
th = atan( (abs(xp2)/abs(xp1)) )
ph = atan( ( sqrt(xp1**2 + xp2**2)/abs(xp3) ) )
vang = sqrt(XEG(gcompete,4)**2*XEG(gcompete,5)**2*XEG(gcompete,6)**2/ & 
(XEG(gcompete,5)**2*XEG(gcompete,6)**2*(cos(th))**2*(sin(ph))**2 + &
XEG(gcompete,4)**2*XEG(gcompete,6)**2*(sin(th))**2*(sin(ph))**2 + &
XEG(gcompete,4)**2*XEG(gcompete,5)**2*(cos(ph))**2 ) )
t1 = d/vang
cellid( N1*N2*(j3-1) + N1*(j2-1) + j1 ) = gcompete(minloc(t1,1))
  enddo !j3
 enddo !j2
enddo !j1
deallocate(gcompete )

! find disconnected grains
c = 0
Lv1 = 1
allocate(badgrain(N1*N2*N3,15) )
badgrain(:,:) = 0

do while (Lv1 .ne. 0)
 Lv1 = 0
 allocate(voxelid(N1*N2*N3))
 voxelid = 0
 c = c+1
 do j = 1,Nsite
  !get the indices of cellid that =1
  allocate(gc(N1*N2*N3))
  gclen = 0
  do j1 = 1,(N1*N2*N3)
   if (j > cellid(j1)-.2 .and. j < cellid(j1)+.2) then
    gclen = gclen + 1
    gc(gclen)  = j1
   endif
  enddo ! j1
  j1 = nint( XEG(j,1)/dx(1) )
  j2 = nint( XEG(j,2)/dx(2) )
  j3 = nint( XEG(j,3)/dx(3) )
  allocate(baccounted(N1*N2*N3))
  blen = 1
  baccounted(blen) = N1*N2*(j3-1)+N1*(j2-1)+j1
  z = blen+1
  do while ( z .ne. blen)
   z = blen
   allocate(potneigh(blen*6))
   do j1 = 1,blen
   potneigh(6*(j1-1)+1:6*j1) = (/ baccounted(j1)+1,baccounted(j1)+N1, &
    baccounted(j1)-1,baccounted(j1)-N1,baccounted(j1)+N1*N2, &
    baccounted(j1)-N1*N2 /)
   enddo ! j1 = 1,blen
   ! get intersection of gc and potneigh
   Nltmp = gclen + 6*blen
   allocate(itmp1(Nltmp))
   call Get_Intersect(gc(1:gclen),gclen,potneigh,6*blen,itmp1,Nltmp,L1)
   deallocate(potneigh)
   ! merge intersect with baccounted
   allocate(itmp2(L1+blen))
   itmp2(1:L1) = itmp1(1:L1)
   itmp2((L1+1):(L1+blen)) = baccounted(1:blen)
   deallocate(itmp1,baccounted)
   ! find unique values of itmp2
   allocate(baccounted(L1+blen))
   call Get_Unique(itmp2,(L1+blen),baccounted,blen)
   deallocate(itmp2)
  enddo ! do while ( z .ne. blen)
  if (gclen .gt. 0 ) then 
  allocate(itmp1(gclen))
  call Get_Difference(gc(1:gclen),gclen,baccounted(1:blen),blen,itmp1,L1)
  if (L1 .gt. 0 ) then
  badgrain(itmp1(1:L1),c) = j
  voxelid((Lv1+1):(L1+Lv1)) = itmp1(1:L1)
  endif 
  Lv1 = L1 + Lv1
  deallocate(gc,itmp1,baccounted)
  else
     deallocate(gc,baccounted)
  endif !  if (gclen .gt. 0 ) then 
 enddo !j
 ! let grains compete for open voxels

 allocate(gc(Nsite))
 gc = (/ (L1, L1=1,Nsite) /)
 do j = 1,Lv1
  j3 = floor( real(voxelid(j)-1)/real(N1*N2) ) +1
  j2 = floor( real(voxelid(j)-1-N1*N2*(j3-1))/real(N1)) + 1 
  j1 = voxelid(j) - N1*N2*(j3-1) - N1*(j2-1)
  allocate(itmp1(Nsite))
 
  call Get_Difference(gc,Nsite,badgrain(voxelid(j),1:c),c,itmp1,L1)
  node = (/ ((real(j1)-.5)*dx(1)), ((real(j2)-.5)*dx(2)), &
            ((real(j3)-.5)*dx(3)) /)
  x1(1:L1) = node(1) - XEG(itmp1(1:L1),1)
  x2(1:L1) = node(2) - XEG(itmp1(1:L1),2)
  x3(1:L1) = node(3) - XEG(itmp1(1:L1),3)
  d(1:L1) = sqrt( x1(1:L1)**2 + x2(1:L1)**2 + x3(1:L1)**2 )
  xp1(1:L1) = XEG(itmp1(1:L1),7)*x1(1:L1) + XEG(itmp1(1:L1),8)*x2(1:L1) + &
              XEG(itmp1(1:L1),9)*x3(1:L1)
  xp2(1:L1) = XEG(itmp1(1:L1),10)*x1(1:L1) + XEG(itmp1(1:L1),11)*x2(1:L1) + &
              XEG(itmp1(1:L1),12)*x3(1:L1)
  xp3(1:L1) = XEG(itmp1(1:L1),13)*x1(1:L1) + XEG(itmp1(1:L1),14)*x2(1:L1) + &
              XEG(itmp1(1:L1),15)*x3(1:L1)
  th(1:L1) = atan( (abs(xp2(1:L1))/abs(xp1(1:L1))) )
  ph(1:L1) = atan( ( sqrt(xp1(1:L1)**2 + xp2(1:L1)**2)/abs(xp3(1:L1)) ) )
  vang(1:L1) = sqrt(XEG(itmp1(1:L1),4)**2*XEG(itmp1(1:L1),5)**2* & 
   XEG(itmp1(1:L1),6)**2/(XEG(itmp1(1:L1),5)**2*XEG(itmp1(1:L1),6)**2* & 
   (cos(th(1:L1)))**2*(sin(ph(1:L1)))**2 + XEG(itmp1(1:L1),4)**2* & 
   XEG(itmp1(1:L1),6)**2*(sin(th(1:L1)))**2*(sin(ph(1:L1)))**2 + &
   XEG(itmp1(1:L1),4)**2*XEG(itmp1(1:L1),5)**2*(cos(ph(1:L1)))**2 ) )
  t1(1:L1) = d(1:L1)/vang(1:L1)
  cellid( N1*N2*(j3-1) + N1*(j2-1) + j1 ) = itmp1(minloc(t1(1:L1),1))
  deallocate(itmp1)
 enddo ! do j = 1,Lv1
 deallocate( gc,voxelid)
enddo ! do while (Lv1 .ne. 0)
deallocate(badgrain,x1,x2,x3,xp1,xp2,xp3, th,ph, vang,d )

! eliminate grains with only one neighbor and give voxels to neighbor
Ngrain = Nsite
if (Nsite .gt. 100) then ! it doesn't make sense to eliminate small domains
allocate(gneighbor(Nsite,Ntot),gneighborL(Nsite))
gneighbor(:,:) = 0
gneighborL(:) = 0
do j1 = 1,(N1-1)
 do j2 = 1,(N2-1)
  do j3 = 1,(N3-1)
   if (cellid(N1*N2*(j3-1)+N1*(j2-1)+j1) .ne. &
       cellid(N1*N2*(j3-1)+N1*(j2-1)+j1+1) ) then
    gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1)) = &
     gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1)) + 1
    gneighbor(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1),gneighborL( & 
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1))) = &
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1+1)  
    gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1+1)) = &
     gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1+1)) + 1
    gneighbor(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1+1),gneighborL( & 
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1+1))) = &
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1) 
   endif
   if (cellid(N1*N2*(j3-1)+N1*(j2-1)+j1) .ne. &
       cellid(N1*N2*(j3-1)+N1*(j2)+j1) ) then
    gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1)) = &
     gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1)) + 1
    gneighbor(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1),gneighborL( & 
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1))) = &
     cellid(N1*N2*(j3-1)+N1*(j2)+j1) 
    gneighborL(cellid(N1*N2*(j3-1)+N1*(j2)+j1)) = &
     gneighborL(cellid(N1*N2*(j3-1)+N1*(j2)+j1)) + 1
    gneighbor(cellid(N1*N2*(j3-1)+N1*(j2)+j1),gneighborL( & 
     cellid(N1*N2*(j3-1)+N1*(j2)+j1))) = &
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1) 
   endif
   if (cellid(N1*N2*(j3-1)+N1*(j2-1)+j1) .ne.&
       cellid(N1*N2*(j3)+N1*(j2-1)+j1) ) then
    gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1)) = &
     gneighborL(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1)) + 1
    gneighbor(cellid(N1*N2*(j3-1)+N1*(j2-1)+j1),gneighborL( & 
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1))) = &
     cellid(N1*N2*(j3)+N1*(j2-1)+j1) 
    gneighborL(cellid(N1*N2*(j3)+N1*(j2-1)+j1)) = &
     gneighborL(cellid(N1*N2*(j3)+N1*(j2-1)+j1)) + 1
    gneighbor(cellid(N1*N2*(j3)+N1*(j2-1)+j1),gneighborL( & 
     cellid(N1*N2*(j3)+N1*(j2-1)+j1))) = &
     cellid(N1*N2*(j3-1)+N1*(j2-1)+j1) 
   endif

  enddo !j3 = 1,(N3-1)
 enddo ! j2 = 1,(N2-1)
enddo ! j1 = 1,(N1-1)

allocate(gc(Nsite),itmp3(1))
gc = (/ (L1, L1=1,Nsite) /)
gclen = Nsite
do j = 1,Nsite
 allocate(itmp1(Ntot))
 call Get_Unique(gneighbor(j,:),Ntot,itmp1,L1)
 if (L1 .eq. 2) then ! use two b/c 0 will be first value in itmp1 
  do j1 = 1,(N1*N2*N3)
   if (j > cellid(j1)-.2 .and. j < cellid(j1)+.2) then
    cellid(j1) = itmp1(L1)    
   endif !   if (j > cellid(j1)-.2 .and. j < cellid(j1)+.2) then
  enddo ! j1
  allocate(itmp2(Nsite))  
  itmp3(1) = j
  call Get_Difference(gc,gclen,itmp3,1,itmp2,L1)
  deallocate(gc)
  allocate(gc(L1))
  gc = itmp2(1:L1)
  gclen = L1
  deallocate(itmp2)
 endif ! (L1 .eq. 2)
  deallocate(itmp1)
enddo ! j = 1,Nsite
Ngrain = gclen
deallocate(gc,gneighbor,gneighborL)
endif ! if (Nsite .gt. 100) then
return
end subroutine ! Get_Morph1
! ============================================================================!
subroutine Get_Difference(A,lenA,B,lenB,AoB,lenAoB)
implicit none
integer j1,j2,lenA,lenB,A(lenA), B(lenB),lenAoB,AoB(lenA),lenAu,lenBu
integer, allocatable :: itmp1(:), itmp2(:)

allocate( itmp1(lenA), itmp2(lenB) )

call Get_Unique(A,lenA,itmp1,lenAu)
call Get_Unique(B,lenB,itmp2,lenBu)

j1 = 1
j2 = 1
lenAoB = 0
do while (j1 .le. lenAu .and. j2 .le. lenBu)
 if (itmp1(j1) .lt. itmp2(j2)) then
   lenAoB = lenAoB + 1
   AoB(lenAoB) = itmp1(j1)
   j1 = j1 + 1
 elseif (itmp2(j2) .lt. itmp1(j1)) then
   j2 = j2+1
 else
  j1 = j1+1
  j2 = j2+1
 endif ! (itmp1(j1) .lt. itmp2(j2)) 
enddo ! while (j1 .lt. gclen .and. j2 .lt. 6*blen) 

if (j2 .ge. lenBu) then
do j2 = j1,lenAu
 lenAoB = lenAoB+1
 AoB(lenAoB) = itmp1(j2)
enddo

endif 



deallocate(itmp1,itmp2)

return
end subroutine
! ============================================================================!
subroutine Get_Unique(A,lenA,Au,lenAu)
implicit none
integer j1,lenA,A(lenA),lenAu,Au(lenA)
integer, allocatable :: itmp1(:)

allocate( itmp1(lenA) )
itmp1 = A
call sort(itmp1,lenA)

lenAu = 0
j1 = 2
do while (j1 .le. lenA)
do while ( itmp1(j1) .eq. itmp1(j1-1) .and. j1 .lt. lenA) 
 j1 = j1 + 1
enddo !while
lenAu = lenAu+1
Au(lenAu) = itmp1(j1-1)
j1 = j1 + 1
enddo ! while

if (lenA .gt. 1) then
if (itmp1(lenA) .ne. itmp1(lenA-1)) then
 lenAu = lenAu+1
 Au(lenAu) = itmp1(lenA)
endif
endif


deallocate(itmp1)

return
end subroutine
! ============================================================================!
subroutine Get_XEG1(N1,N2,N3,dx,Nsite,Ngamma,Ntsample,R,XEG,param_mat_CDF, &
           Gammainv,d5bound,irg)

use random

!                 XEG parameters
! XEG(:,1), XEG(:,2), XEG(:,3) ..... grain nucleation sites
! XEG(:,4), XEG(:,5), XEG(:,6) ..... grain velocities
! XEG(:,7), XEG(:,8), XEG(:,9) ..... map to local x1 coor
! XEG(:,10), XEG(:,11), XEG(:,12)... map to local x2 coor
! XEG(:,13), XEG(:,14), XEG(:,15)... map to local x3 coor
! XEG(:,16), XEG(:,17), XEG(:,18)... crystallographic orientation
implicit none

integer irg, j,Nsite,N1,N2,N3,Nsite2,j1,Ngamma,Ntsample,status, &
        iloc,j2,minN5
integer, allocatable :: iz(:)
real*8 XEG(Nsite,18),dx(3),R, varphi1,Phi,varphi2,Gammainv(Ngamma), &
       param_mat_CDF(Ntsample,18),d5bound(5),bound,q
real*8, allocatable :: Xtmp(:,:),Xtmp2(:),d5(:), &
       Ynor(:)

if (irg .eq. 0 ) then ! read in from file

!open(unit=10,file='EGT_Parameters_Simulation.txt',status='old')

open(unit=10,file='XEG_test.txt',status='old')

do j=1,Nsite
read(10,*) XEG(j,1),XEG(j,2),XEG(j,3),XEG(j,4),XEG(j,5),XEG(j,6), &
          XEG(j,7),XEG(j,8),XEG(j,9),XEG(j,10),XEG(j,11),XEG(j,12), &
          XEG(j,13),XEG(j,14),XEG(j,15),XEG(j,16),XEG(j,17),XEG(j,18)
enddo ! j = 1,Nsite

close(10)

else  ! generate realizations

minN5 = minval( (/ Nsite, 6 /) )
! this version of the code generates a matern hard core Poisson process
! I am not thinning, I am just looping until the hard core condition is
! satisfied, which is easy since hard core is small

allocate(Xtmp(Nsite,3),Xtmp2(Nsite-1),iz(Nsite-1),d5(Nsite**2))
Nsite2 = 0
do while (Nsite2 .ne. Nsite)
Nsite2 = 0
do j = 1,Nsite
Xtmp(j,:) = (/ rand(0)*real(N1-1)*dx(1),rand(0)*real(N2-1)*dx(2), &
               rand(0)*real(N3-1)*dx(3) /)
enddo ! j = 1,Nsite
do j = 1,Nsite
iz = (/ (j1, j1=1,(j-1)), (j1,j1=(j+1),Nsite) /)
Xtmp2(:) = sqrt( (Xtmp(j,1)-Xtmp(iz,1))**2 + (Xtmp(j,2)-Xtmp(iz,2))**2 + &
    (Xtmp(j,3)-Xtmp(iz,3))**2)
d5( (/ (j1,j1=(Nsite*(j-1)+1),(Nsite*j))/)) = sqrt( (Xtmp(j,1)-Xtmp(:,1))**2 + &
    (Xtmp(j,2)-Xtmp(:,2))**2 + (Xtmp(j,3)-Xtmp(:,3))**2) 
if (minval(Xtmp2) .gt. R ) then
Nsite2 = Nsite2 + 1
endif ! if (minval(Xtmp2) .gt. R ) then
enddo ! j = 1,Nsite
enddo ! do while (Nsite2 .ne. Nsite)

XEG(:,1:3) = Xtmp
deallocate(Xtmp,Xtmp2,iz)
! determine conditional group each grain belongs in
allocate(Xtmp2(Nsite))
do j = 1,Nsite
 call sort_real(d5( Nsite*(j-1)+1:Nsite*j),Nsite)
 Xtmp2(j) = &
sum(d5( (/ (j1,j1=(Nsite*(j-1)+1),(Nsite*(j-1)+MinN5))/)))/real(MinN5)
enddo
deallocate(d5)
allocate(d5(Nsite))
do j = 1,Nsite
 if (Xtmp2(j) .le. d5bound(2)) d5(j) = 1 
 if (Xtmp2(j) .le. d5bound(3) .and. Xtmp2(j) .gt. d5bound(2)) d5(j) = 2 
 if (Xtmp2(j) .le. d5bound(4) .and. Xtmp2(j) .gt. d5bound(3)) d5(j) = 3 
 if (Xtmp2(j) .le. d5bound(5) .and. Xtmp2(j) .gt. d5bound(4)) d5(j) = 4 
 if (Xtmp2(j) .gt. d5bound(5)) d5(j) = 5 
enddo ! do j = 1,Nsite
deallocate(Xtmp2)
! end determine conditional group each grain belongs in
allocate(Xtmp2(9),Ynor(9))

! simulate marks according to distributions obtained from IN100 data set
do j = 1,Nsite
do j1 = 1,9
 Xtmp2(j1) = random_normal()
enddo ! do j1 = 1,9
Ynor(:) = 0.0
do j1 = 1,9
do j2 = 1,9
Ynor(j1) = Xtmp2(j2)*Gammainv( 81*(int(d5(j))-1)+ 9*(j2-1)+j1)+Ynor(j1)
enddo ! j2 = 1,9
enddo ! j1 = 1,9
do j1 = 1,9
call cdfnor(1,Xtmp2(j1),q,Ynor(j1),0.0d0,1.0d0,status,bound)
 iloc = minloc( abs(param_mat_CDF( &
  (/ (j2,j2=1024*(int(d5(j))-1)+1,1024*int(d5(j)))/),2*j1)-Xtmp2(j1)),1)
Xtmp2(j1) =  param_mat_CDF(iloc,2*j1-1)
enddo ! do j1 = 1,9

XEG(j,4) = Xtmp2(1)
XEG(j,5) = Xtmp2(2)
XEG(j,6) = Xtmp2(3)
varphi1 = Xtmp2(4)
Phi = Xtmp2(5)
varphi2 = Xtmp2(6)
XEG(j,7) = cos(varphi2)*cos(varphi1) -  sin(varphi2)*cos(Phi)*sin(varphi1)
XEG(j,8) = -sin(varphi2)*cos(varphi1) - cos(varphi2)*cos(Phi)*sin(varphi1)
XEG(j,9) = sin(Phi)*sin(varphi1)
XEG(j,10) = cos(varphi2)*sin(varphi1)+sin(varphi2)*cos(Phi)*cos(varphi1)
XEG(j,11) = -sin(varphi2)*sin(varphi1)+cos(varphi2)*cos(Phi)*cos(varphi1)
XEG(j,12) = -sin(Phi)*cos(varphi1)
XEG(j,13) = sin(varphi2)*sin(Phi)
XEG(j,14) = cos(varphi2)*sin(Phi)
XEG(j,15) = cos(Phi)
XEG(j,16) = Xtmp2(7)
XEG(j,17) = Xtmp2(8)
XEG(j,18) = Xtmp2(9)

enddo ! j = 1,Nsite

deallocate(Xtmp2,d5,Ynor)
endif ! if (irg .eq. 0 ) ! read in from file



return
end subroutine ! Get_XEG1
! ============================================================================!
subroutine Get_Intersect(A,lenA,B,lenB,AnB,lenAB,lenAnB)
implicit none
integer j1,j2,lenA,lenAB,lenB,A(lenA), B(lenB),AnB(lenAB),lenAnB
integer, allocatable :: itmp1(:), itmp2(:)

allocate( itmp1(lenA), itmp2(lenB) )
itmp1 = A
itmp2 = B
call sort(itmp1,lenA)
call sort(itmp2,lenB)
j1 = 1
j2 = 1
lenAnB = 0
do while (j1 .le. lenA .and. j2 .le. lenB)
 if (itmp1(j1) .lt. itmp2(j2)) then
   j1 = j1 + 1
 elseif (itmp1(j1) .gt. itmp2(j2)) then
   j2 = j2+1
 else
  lenAnB = lenAnB + 1
  AnB(lenAnB) = itmp1(j1)
  j1 = j1+1
  j2 = j2+1
!write(*,*) lenAnB
 endif ! (itmp1(j1) .lt. itmp2(j2)) 
enddo ! while (j1 .lt. gclen .and. j2 .lt. 6*blen) 

deallocate(itmp1,itmp2)

return
end subroutine
! ============================================================================!
! THE FOLLOWING ROUTINES ARE TO SORT AN ARRAY AND TAKEN FROM
! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      INTEGER, INTENT(INOUT) :: a, b
      INTEGER                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

SUBROUTINE sort(arr,n)
INTEGER n,M,NSTACK,arr(n)
PARAMETER (M=7,NSTACK=50)
! Sorts an array arr(1:n) into ascending numerical order using the Quicksort
! algorithm. n is input; arr is replaced on output by its sorted rearrangement.
! Parameters: M is the size of subarrays sorted by straight insertion and 
! NSTACK is the required auxiliary storage.
INTEGER i,ir,j,jstack,k,l,istack(NSTACK),a,temp
jstack=0
l=1
ir=n
1 if(ir-l.lt.M)then ! Insertion sort when subarray small enough.
do j=l+1,ir
a=arr(j)
do i=j-1,l,-1
if(arr(i).le.a) goto 2
arr(i+1)=arr(i)
enddo 
i=l-1
2 arr(i+1)=a
enddo 
if(jstack.eq.0) return
ir=istack(jstack) !Pop stack and begin a new round of partitioning.
l=istack(jstack-1)
jstack=jstack-2
else
k=(l+ir)/2 
temp=arr(k)
arr(k)=arr(l+1)
arr(l+1)=temp
if(arr(l).gt.arr(ir))then
temp=arr(l)
arr(l)=arr(ir)
arr(ir)=temp
endif
if(arr(l+1).gt.arr(ir))then
temp=arr(l+1)
arr(l+1)=arr(ir)
arr(ir)=temp
endif
if(arr(l).gt.arr(l+1))then
temp=arr(l)
arr(l)=arr(l+1)
arr(l+1)=temp
endif
i=l+1 !Initialize pointers for partitioning.
j=ir
a=arr(l+1) !Partitioning element.
3 continue !Beginning of innermost loop.
i=i+1 !Scan up to nd element > a.
if(arr(i).lt.a)goto 3
4 continue
j=j-1 !Scan down to nd element < a.
if(arr(j).gt.a)goto 4
if(j.lt.i)goto 5 !Pointers crossed. Exit with partitioning complete.
temp=arr(i) !Exchange elements.
arr(i)=arr(j)
arr(j)=temp
goto 3 !End of innermost loop.
5 arr(l+1)=arr(j) !Insert partitioning element.
arr(j)=a
jstack=jstack+2
!if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
if(ir-i+1.ge.j-l)then
istack(jstack)=ir
istack(jstack-1)=i
ir=j-1
else
istack(jstack)=j-1
istack(jstack-1)=l
l=i
endif
endif
goto 1
return
end subroutine

! END
! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
! ============================================================================!
! ============================================================================!
SUBROUTINE sort_real(arr,n)
INTEGER n,M,NSTACK
real*8 arr(n)
PARAMETER (M=7,NSTACK=50)
! Sorts an array arr(1:n) into ascending numerical order using the Quicksort
! algorithm. n is input; arr is replaced on output by its sorted rearrangement.
! Parameters: M is the size of subarrays sorted by straight insertion and 
! NSTACK is the required auxiliary storage.
INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
real*8 a,temp
jstack=0
l=1
ir=n
1 if(ir-l.lt.M)then ! Insertion sort when subarray small enough.
do j=l+1,ir
a=arr(j)
do i=j-1,l,-1
if(arr(i).le.a) goto 2
arr(i+1)=arr(i)
enddo 
i=l-1
2 arr(i+1)=a
enddo 
if(jstack.eq.0) return
ir=istack(jstack) !Pop stack and begin a new round of partitioning.
l=istack(jstack-1)
jstack=jstack-2
else
k=(l+ir)/2 
temp=arr(k)
arr(k)=arr(l+1)
arr(l+1)=temp
if(arr(l).gt.arr(ir))then
temp=arr(l)
arr(l)=arr(ir)
arr(ir)=temp
endif
if(arr(l+1).gt.arr(ir))then
temp=arr(l+1)
arr(l+1)=arr(ir)
arr(ir)=temp
endif
if(arr(l).gt.arr(l+1))then
temp=arr(l)
arr(l)=arr(l+1)
arr(l+1)=temp
endif
i=l+1 !Initialize pointers for partitioning.
j=ir
a=arr(l+1) !Partitioning element.
3 continue !Beginning of innermost loop.
i=i+1 !Scan up to nd element > a.
if(arr(i).lt.a)goto 3
4 continue
j=j-1 !Scan down to nd element < a.
if(arr(j).gt.a)goto 4
if(j.lt.i)goto 5 !Pointers crossed. Exit with partitioning complete.
temp=arr(i) !Exchange elements.
arr(i)=arr(j)
arr(j)=temp
goto 3 !End of innermost loop.
5 arr(l+1)=arr(j) !Insert partitioning element.
arr(j)=a
jstack=jstack+2
!if(jstack.gt.NSTACK)pause 'NSTACK too small in sort'
if(ir-i+1.ge.j-l)then
istack(jstack)=ir
istack(jstack-1)=i
ir=j-1
else
istack(jstack)=j-1
istack(jstack-1)=l
l=i
endif
endif
goto 1
return
end subroutine ! sort_real

! END 
! http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap08/sorting.f90
! ============================================================================!
subroutine Write_Microstructure_VTI(N1,N2,N3,dx,cellid,filname)
implicit none

character filname*150
integer j1,cellid(N1*N2*N3),Ntot,N1,N2,N3
real*8 dx(3)

Ntot = N1*N2*N3
open(unit=99,file=filname,status='replace')

write(99,1) 
write(99,2) N1,N2,N3,dx(1),dx(2),dx(3)

write(99,3) N1,N2,N3
write(99,4) 
write(99,5) 
write(99,6) 
write(99,7) 

do j1 = 1,Ntot
write(99,'(I3)') cellid(j1)
enddo

write(99,8) 
write(99,9) 
write(99,10) 
write(99,11) 
write(99,12) 

close(99)


1 format('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">')
2 format('<ImageData WholeExtent="0 ',i3' 0 ',i3,' 0 ',i3,&
          '" Origin="0 0 0" Spacing="',2(f5.3,1x),f5.3,'">')
3 format('<Piece Extent="0 ',i3,' 0 ',i3,' 0 ',i3'">')
4 format('<PointData>')
5 format('</PointData>')
6 format('<CellData>')
7 format('<DataArray type="Float32" Name="Region" format="ascii">')
8 format('</DataArray>')
9 format('</CellData>')
10 format('</Piece>')
11 format('</ImageData>')
12 format('</VTKFile>')

return
end subroutine
! ============================================================================!
