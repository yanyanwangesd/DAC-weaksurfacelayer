subroutine find_hack (geometry,network,stack,params)


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (parm) params
  integer i,j,k
  double precision l
  double precision,dimension(:),allocatable::contributing_area, max_length, horton_length
  
  !!!!!!@ Y.Wang, add to output every desired timesteps
  integer icount
  character cs*8
  icount=params%istep
  write(cs,'(I8)') icount
  if (icount.lt.10)      cs(1:7)='0000000'
  if (icount.lt.100)     cs(1:6)='000000'
  if (icount.lt.1000)    cs(1:5)='00000'
  if (icount.lt.10000)   cs(1:4)='0000'
  if (icount.lt.100000)  cs(1:3)='000'
  if (icount.lt.1000000)  cs(1:2)='00'
  if (icount.lt.10000000)  cs(1:1)='0'
  !!!!!!@ Y.Wang, add to output every desired timesteps
  
  

  call time_in ('find_hack')

  allocate(contributing_area(geometry%nnode), max_length(geometry%nnode), horton_length(geometry%nnode))

  contributing_area=geometry%surface
  max_length=0.
  do i=stack%nnode,1,-1
     j=stack%order(i)
     k=network%receiver(j)
     if (k.ne.0) then
        contributing_area(k) = contributing_area(k) + contributing_area(j)
        l = dsqrt((geometry%x(k)-geometry%x(j))**2.d0 + (geometry%y(k)-geometry%y(j))**2.d0)
        max_length(k) = max(max_length(k),max_length(j)+l)
     endif
  enddo

  !!!!!!@ Y.Wang, modified to store in ASCII folder and name files with timesteps
  open (40,file='ASCII/HackData'//cs,status='unknown')
  !!!!!!@ Y.Wang, modified to store in ASCII folder and name files with timesteps
  do i=1,geometry%nnode
     write(40,'(2i7,2f16.3,i3)') i, geometry%catchment(i), contributing_area(i), max_length(i), geometry%strahler(i)
  enddo
  close (40,err=1232)
  
  !!!!!!@ Y.Wang, modified to store in ASCII folder and name files with timesteps
  open (41,file='ASCII/HortonData'//cs,status='unknown')
  !!!!!!@ Y.Wang, modified to store in ASCII folder and name files with timesteps
  horton_length = 0.
  do i=stack%nnode,1,-1
      j=stack%order(i)
      if (network%receiver(j) .eq. 0) then
         if (horton_length(j).gt.0.d0) then
            write(41,'(2i7,2f16.3,i3)') j, geometry%catchment(j), contributing_area(j), horton_length(j), geometry%strahler(j)
         endif
      else
         k=network%receiver(j)
         l = dsqrt((geometry%x(k)-geometry%x(j))**2.d0 + (geometry%y(k)-geometry%y(j))**2.d0)
         if (geometry%strahler(j).eq.geometry%strahler(k)) then    
            horton_length(k) = horton_length(j) + l
         else
            horton_length(j) =  horton_length(j) + l
            write(41,'(2i7,2f16.3,i3)') j, geometry%catchment(j), contributing_area(j), horton_length(j), geometry%strahler(j)
         endif
      endif
   enddo
  close (41,err=1232)
  deallocate(contributing_area, max_length,horton_length)
  call time_out ('find_hack')

  return
1232 STOP 'problem closing HackData'
end subroutine find_hack
