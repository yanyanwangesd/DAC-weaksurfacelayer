subroutine flexure2D (geometry, params)

! Routine to calculate flexure uplift due to erosion
! Algorithm is from Beaumont (1978), solving the zero-order modified Bessel function 
! of the second kind. The Bessel function is solved by call ZBESK (see the document of zbesh.f90)

! Load of erosion is taken as the product of erosion rate, and the surface area of river nodes

! Added in Oct. 2021 by Yanyan Wang

!-----------------------------------------------------------------------------
  use definitions
  implicit none
  type (geom) geometry
  type (parm) params

  double precision poisson, g, young, d, alpha
  double precision rhom, rhoa, rhoc, pi
  double precision rhor, kr
  double precision pres

  integer i, j
  integer nnode
  double precision delz, deltal
  double precision, dimension(:), allocatable:: p, w

  integer IERR,  KODE, N,  NZ
  double precision  ZI, ZR
  double precision  FNU
  double precision  CYI, CYR
  double precision kei
  
  
  
  ! constants used in flexure calcualtion
  rhom = 3150 !density mantle, kg/m^3
  rhoa = 1 ! density air, kg/m^3
  rhoc = 2800 ! density of the crust, unweathered rock, kg/m^3

  poisson = 0.25 ! Poisson ratio
  g = 9.81 ! gravity, m/s^2
  young = 1.d11 ! Young's modulus, Pa
  pi = 3.14159265359
    
  d = young*(params%flexure_te**3)/12/(1-poisson**2)
  alpha = (4*d/(rhom-rhoa)/g)**0.25; 

  
  ! density of regolith layer scales with the erodibility contrast used in find_erodibility.f90
  ! minimum density of regolity layer is set to be 1000 kg/m^3
  kr = maxval(geometry%k)

  if ((1-kr/params%k_scalar1/100) < 1000/rhoc) then
      rhor = rhoc*(1-kr/params%k_scalar1/100)
  else
      rhor = 1000
  endif
  
  
  ! constants used for solving the BesselK function (the modified bessel function of the second kind) calcualtion
  ! don't change the three parameters!!!
  FNU = 0
  KODE = 1
  N = 1
  
  
  
  nnode = geometry%nnode
  allocate(p(nnode), w(nnode)) ! weight of erosion and deflection
  
call time_in ('flexure2D')

  if (params%time >= params%flexure_time ) then
  
      ! calcualte  the load of each node
      do i = 1,nnode
          delz = geometry%erosion_rate(i)*params%deltat
          ! note here I used half of the area*reosion rate as load weight. negative means it is the eroded weight, kg @YWang
          if (geometry%dregolith(i) .gt. 0) then
              if (delz .lt. geometry%dregolith(i)) then
                  p(i) = -geometry%surface(i)*delz*rhor*g/2 ! only regolith layer is eroded away
              else
                  p(i) = -(geometry%surface(i)*geometry%dregolith(i)*rhor+(delz-geometry%dregolith(i))*rhoc)*g/2 !regolith layer and bedrock are eroded
              endif
          else
              p(i) = -geometry%surface(i)*delz*rhoc*g/2 !only bedrock is eroded
          endif

      enddo
  
      ! calculate deflection from erosion
      do i = 1,nnode
          
          w(i) = 0
          do j = 1,nnode
              deltal = dsqrt((geometry%x(i)-geometry%x(j))**2+(geometry%y(i)-geometry%y(j))**2)
              ! calculate BesselK function
              ZR =  deltal/alpha
              ZI =  deltal/alpha
              call ZBESK(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
              kei = -CYR
              pres = p(j)/(pi*(rhom-rhoa)*g*alpha**2)
              w(i) = w(i)+pres*kei
          enddo
          ! update z
          geometry%z(i) = geometry%z(i)+w(i)
          ! save flexure uplift rate
          geometry%flexure(i) = w(i)/params%deltat
          
      enddo


  endif ! end of time loop
  
  
  
deallocate(p, w)


call time_out ('flexure2D')

return

end subroutine flexure2D


