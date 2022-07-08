subroutine output_z_tau (params,geometry,network,stack)


  use definitions

  implicit none

  type (geom) geometry
  type (netw) network
  type (stck) stack
  type (parm) params
  integer i,j,k,icount,ii,node1,node2
  double precision l,ks_K, slope,Anot, xc 
  double precision xdiv,ydiv,xnode1,ynode1,xnode2,ynode2,xhead1,yhead1,xhead2,yhead2,flux_head1,flux_head2,chi_head1,chi_head2
  character cs*8
  double precision,dimension(:),allocatable::contributing_area, tau, faketau,chi
  
  !!!!!!add to write channel head elevation @Y.Wang, Aug. 2021 !!!!!!!!!!!!!!!
  double precision zhead1, zhead2, tau2, ave_erosion, fact1,fact2, kappa
  integer num_samples, num_full_bin,sample_in_last_bin, h
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  


  call time_in ('output_z_tau')

  !Anot = 1.d7
  Anot = 1.
  icount=params%istep

  write(cs,'(I8)') icount
  if (icount.lt.10)      cs(1:7)='0000000'
  if (icount.lt.100)     cs(1:6)='000000'
  if (icount.lt.1000)    cs(1:5)='00000'
  if (icount.lt.10000)   cs(1:4)='0000'
  if (icount.lt.100000)  cs(1:3)='000'
  if (icount.lt.1000000)  cs(1:2)='00'
  if (icount.lt.10000000)  cs(1:1)='0'


! !!! @YW, the chi is not geometry%chi, the chi has uplift and discharge in there, more like a chi* in Eqn.(5) of Willett et al 2014



  allocate(contributing_area(geometry%nnode), tau(geometry%nnode),faketau(geometry%nnode),chi(geometry%nnode))

  contributing_area=geometry%surface
  do i=stack%nnode,1,-1
     j=stack%order(i)
     k=network%receiver(j)
     if (k.ne.0) then
        contributing_area(k) = contributing_area(k) + contributing_area(j)
     endif
  enddo
  ks_K=params%k_scalar1*params%rainfall_height**params%m
  tau=0.
  faketau=0.
  do i=1,stack%nnode
     j=stack%order(i)
     if (network%receiver(j).eq.0) then
        tau(j)=0.
        faketau(j)=0.
        chi(j) = 0.
        geometry%chi(j) = 0.
     else                 
        k=network%receiver(j)
        l = dsqrt((geometry%x(k)-geometry%x(j))**2.d0 + (geometry%y(k)-geometry%y(j))**2.d0)
        slope = (geometry%z(j)-geometry%z(k))/l
        tau(j)=tau(k)+(1./(ks_K*contributing_area(j)**params%m*(slope**(params%n - 1.))))*l
        faketau(j)=faketau(k)+(1./(ks_K*contributing_area(j)**params%m))*l
        chi(j) = chi(k) + l*(geometry%w(j)**(1/params%n))*&
             &(Anot/geometry%discharge(j))**(params%m/params%n)
             
        geometry%chi(j) = geometry%chi(k)+l*(Anot/(contributing_area(j)))**params%m
     endif
  enddo
  open (75,file='ASCII/tau_z'//cs,status='unknown')
  do i=1,geometry%nnode
     write(75,'(f8.3,f15.3,f15.4,f15.3,i5,f15.3)') geometry%z(i), tau(i),  chi(i), geometry%chi(i), geometry%catchment(i), contributing_area(i)
  enddo
  close (75)

  open (76,file='ASCII/channel_head'//cs,status='unknown')
  open (77,file='ASCII/delta_tau'//cs,status='unknown')
  open (78,file='ASCII/no_channel_connection'//cs,status='unknown')
  open (79,file='ASCII/divides'//cs,status='unknown')
  xc = params%xc
  do i =1,geometry%ndivide
     xdiv = geometry%xdiv(i)
     ydiv = geometry%ydiv(i)
     !first node of the divide
     node1 = geometry%nndivnode(i,1)
     xnode1 = geometry%x(node1)
     ynode1 = geometry%y(node1)
     l = dsqrt((xnode1-xdiv)**2.d0 + (ynode1-ydiv)**2.d0) ! distance between bounding node and divide
     if (l.gt.xc) then
     
         l = l-xc ! distance between bounding node and channel head
         
        !!!!!!add to write channel head elevation @Y.Wang, Aug. 2021 !!!!!!!!!!!!!!!
         if (params%transient_divide) then
             
             kappa = (geometry%k(node1)*(geometry%precipitation(node1)*params%ka)**(params%m))
             if (params%h*params%m.eq.1) then
                 tau2=log((l+xc)/params%xc)/kappa
             else
                 tau2= ((l+xc)**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m)) /&
                     &kappa/(1.d0-params%h*params%m)
             endif
             num_samples=INT(min(tau2,params%time)/params%deltat)
             if (num_samples.eq.0) then
                 fact1=(geometry%erosion_rate(node1)/geometry%k(node1)/ &
                     &(geometry%precipitation(node1)*params%ka)**params%m)**(1.d0/params%n)
             else
                 if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                     num_full_bin=params%num_bins
                     sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                 else
                     num_full_bin = INT(num_samples/params%sample_per_bin)+1
                     sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                 endif
                 ave_erosion=0.d0
                 if (num_full_bin.eq.1) then
                     ave_erosion=geometry%erosion_rate_history(node1,1)
                 else
                     ave_erosion=0.d0
                     do h=1,num_full_bin-1
                         ave_erosion=ave_erosion +geometry%erosion_rate_history(node1,h)*params%sample_per_bin
                     enddo
                     ave_erosion=ave_erosion +geometry%erosion_rate_history(node1,num_full_bin)*sample_in_last_bin
                     ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                 endif
                 fact1=(ave_erosion/kappa)**(1.d0/params%n)
             endif
             
         else
             fact1=(geometry%erosion_rate(node1)/kappa)**(1.d0/params%n)
         endif
         
         
         if (params%hmn.eq.0.d0) then
             fact2=log((l+xc)/xc)
         else
             fact2=(1.d0/params%hmn)*((l+xc)**(params%hmn)-xc**(params%hmn))
         endif
         
         ! channel head elevation z = fact1*fact2 according to equation (7) of Goren et al.2014
         zhead1 = fact1*fact2+geometry%z(node1)
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         
         if (xnode1.eq.xdiv) then
             xhead1 = xnode1
             if (ynode1.gt.ydiv) then
                 yhead1 = ydiv + xc
             else
                 yhead1 = ydiv - xc
             endif
         elseif (ynode1.eq.ydiv) then
             yhead1 = ynode1
             if (xnode1.gt.xdiv) then
                 xhead1 = xdiv+xc
             else
                 xhead1 = xdiv-xc
             endif
         else
             slope = (ydiv-ynode1)/(xdiv-xnode1)
             if (xnode1.gt.xdiv) then 
                 xhead1 = xdiv + xc/dsqrt(1+slope**2)
             else
                 xhead1 = xdiv - xc/dsqrt(1+slope**2)
             endif
             yhead1 = slope*(xhead1-xdiv) + ydiv
         endif
         
        flux_head1 = (params%ka*xc**params%h)*geometry%precipitation(node1)
        chi_head1 = chi(node1) +   l*(geometry%w(node1)**(1/params%n))*&
             &(Anot/flux_head1)**(params%m/params%n)
        write(76,'(f15.3,f15.3,f15.3,f15.3,f15.3)') xhead1, yhead1, zhead1, flux_head1, chi_head1
        
     else ! this is a short divide, 
         !!!!!!add to write channel head elevation @Y.Wang, Aug. 2021 !!!!!!!!!!!!!!!
         zhead1 = geometry%z(node1)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         chi_head1 = chi(node1)
         write(76,'(f15.3,f15.3,f15.3,f15.3,f15.3)') xnode1, ynode1, zhead1, geometry%discharge(node1), chi_head1
         
     endif


     !second node of the divide
     node2 = geometry%nndivnode(i,2)
     xnode2 = geometry%x(node2)
     ynode2 = geometry%y(node2)    
     l = dsqrt((xnode2-xdiv)**2.d0 + (ynode2-ydiv)**2.d0)
     if (l.gt.xc) then
        l = l-xc
        
        !!!!!!add to write channel head elevation @Y.Wang, Aug. 2021 !!!!!!!!!!!!!!!
         if (params%transient_divide) then
             kappa = geometry%k(node2)*(geometry%precipitation(node2)*params%ka)**(params%m)
             if (params%h*params%m.eq.1) then
                 tau2=log((l+xc)/xc)/kappa
             else
                 tau2= ((l+xc)**(1.d0-params%h*params%m) - params%xc**(1.d0-params%h*params%m))/kappa/(1.d0-params%h*params%m)
             endif
             num_samples=INT(min(tau2,params%time)/params%deltat)
             if (num_samples.eq.0) then
                 fact1=(geometry%erosion_rate(node2)/kappa)**(1.d0/params%n)
             else
                 if ((params%num_bins-1)*params%sample_per_bin.lt.num_samples) then
                     num_full_bin=params%num_bins
                     sample_in_last_bin=num_samples-(params%num_bins-1)*params%sample_per_bin
                 else
                     num_full_bin = INT(num_samples/params%sample_per_bin)+1
                     sample_in_last_bin=num_samples-(num_full_bin-1)*params%sample_per_bin
                 endif
                 ave_erosion=0.d0
                 if (num_full_bin.eq.1) then
                     ave_erosion=geometry%erosion_rate_history(node2,1)
                 else
                     ave_erosion=0.d0
                     do h=1,num_full_bin-1
                         ave_erosion=ave_erosion +geometry%erosion_rate_history(node2,h)*params%sample_per_bin
                     enddo
                     ave_erosion=ave_erosion +geometry%erosion_rate_history(node2,num_full_bin)*sample_in_last_bin
                     ave_erosion=ave_erosion/((num_full_bin-1)*params%sample_per_bin + sample_in_last_bin)
                 endif
                 fact1=(ave_erosion/kappa)**(1.d0/params%n)
             endif
             
         else
             fact1=(geometry%erosion_rate(node2)/kappa)**(1.d0/params%n)
         endif
         
         
         if (params%hmn.eq.0.d0) then
             fact2=log((l+xc)/xc)
         else
             fact2=(1.d0/params%hmn)*((l+xc)**(params%hmn)-(params%xc)**(params%hmn))
         endif
         
         ! channel head elevation z = fact1*fact2 according to equation (7) of Goren et al.2014
         zhead2 = fact1*fact2+geometry%z(node2)
         
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
         
        if (xnode2.eq.xdiv) then
           xhead2 = xnode2
           if (ynode2.gt.ydiv) then
              yhead2 = ydiv + xc
           else
              yhead2 = ydiv-xc
           endif
        elseif (ynode2.eq.ydiv) then
           yhead2 = ynode2
           if (xnode2.gt.xdiv) then
              xhead2 = xdiv+xc
           else
              xhead2 = xdiv-xc
           endif
        else
           slope = (ydiv-ynode2)/(xdiv-xnode2)
           if (xnode2.gt.xdiv) then
              xhead2 = xdiv + xc/dsqrt(1+slope**2)
           else
              xhead2 = xdiv - xc/dsqrt(1+slope**2)
           endif
           yhead2 = slope*(xhead2-xdiv) + ydiv
        endif
        flux_head2 = (params%ka*xc**params%h)*geometry%precipitation(node2)
        chi_head2 = chi(node2) +   l*(geometry%w(node2)**(1/params%n))*&
             &(Anot/flux_head1)**(params%m/params%n)
        write(76,'(f15.3,f15.3,f15.3,f15.3,f15.3)') xhead2, yhead2, zhead2, flux_head2, chi_head2
     else ! this is a short divide
     
         !!!!!!add to write channel head elevation @Y.Wang, Aug. 2021 !!!!!!!!!!!!!!!
        zhead2 = geometry%z(node2)
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         
        chi_head2 = chi(node2)
        
        
        write(76,'(f15.3,f15.3,f15.3,f15.3,f15.3)') xnode2, ynode2, zhead2, geometry%discharge(node2), chi_head2
        
     endif
     write(77,'(f15.3)') abs(chi_head2-chi_head1)
     write(78,'(i10,i10)') node1, node2
     write(79,'(f15.3,f15.3,f15.3)') xdiv,ydiv,geometry%zdiv(i)
  enddo

  close (76)
  close (77)
  close (78)
  close (79)

  deallocate(contributing_area, tau, faketau,chi)


  if(params%istep.eq.0)geometry%chi=0

  call time_out ('output_z_tau')

  return
end subroutine output_z_tau
