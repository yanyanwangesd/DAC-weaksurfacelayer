subroutine set_params_as_function (geometry,params)

  ! Routine to ...

  use definitions

  implicit none

  type (geom) geometry
  type (parm) params

  integer i,j,k,l, nset, icond, inbound
  double precision  d1,d2,d3, fmin,fmax, fwk, foffs, tmin,tmax
  double precision,pointer:: fdep(:), depvar(:)
  double precision tinterval

!   ! legacy line of code that is needed so that find_precipitation does not need to be called  - set anyway where dependent variable is not P
!   geometry%precipitation=params%rainfall_height
!   ! and the resetting part from uplift_and_advect
!   geometry%u = 0.
!   geometry%v = 0.


  if(.NOT.params%f_varies_with_xyz)then
  
      geometry%w=params%uplift_scalar1
      
      else   
          ! add pieces of precheck whether uplift is set in the setting blocks, if uplift is set in block, 
          ! but there is gap of time intervals, the gap of time interval is set to uplift_scaler1 as default uplift @Y.Wang on May 25, 2021
          if(params%time==0) then
              print*, '  '
              print*, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - '
              print*, 'Precheck of uplift specification in setting blocks:'                                      
              tinterval = 0
              do nset=1, params%f_num_sets             
                 if(params%f_variable_determined(nset)==2) then                
                     tinterval = tinterval+params%f_timebound(nset,2)-params%f_timebound(nset,1)
                 endif          
              enddo
              if(tinterval < params%tfinal)then
                  print*, '  Warning!!! There is GAP of time interval for uplift specification.' 
                  print*, '  uplift_scalar1 is used for gaps until gaps are rewrite with its previous uplift block'                   
                  geometry%w=params%uplift_scalar1
                                          
              else if(tinterval > params%tfinal)then
                  print*, '  Warning!!! There is OVERLAP between time intervals for uplift specification.'   
              
              else
                  print*, '  Great!!! Uplift uplift specification covers all time intervals and no overlap.'                    
              endif  
              
              print*, '- - - - - - - - - - - - - - - - - - - - - - - - - - - - '
              print*, '  '
              print*, '  '
                          
          endif     
      
          
      
      
          do nset=1, params%f_num_sets
              ! independent variable, select, and set pointer to
              select case(params%f_depends_on(nset))
                case(0)
                    fdep=>geometry%x
                case(1)
                    fdep=>geometry%y
                case(2)
                    fdep=>geometry%z
              end select
        
              ! dependent variable, select, and set pointer to
              select case(params%f_variable_determined(nset))
                case(0)
                    depvar=>geometry%u
                case(1)
                    depvar=>geometry%v
                case(2)
                    depvar=>geometry%w
                case(3)
                    depvar=>geometry%precipitation
                case(4)
                    depvar=>geometry%k
                case(5)
                    depvar=>geometry%z
              end select
        
                tmin = params%f_timebound(nset,1)
                tmax = params%f_timebound(nset,2)
                if(tmin.eq.0 .and. tmax.eq.0) tmax = 1.d18
                
                if(  params%time.ge.tmin  .and. params%time.le.tmax  )then
                        do i = 1, geometry%nnode
                            ! which condition segment (line) applies
                            inbound=0
                            do icond=1,params%f_polyseg(nset)
                                fmin=params%poly(nset,icond,1)
                                fmax=params%poly(nset,icond,2)
                                foffs=params%poly(nset,icond,3)
                                if(fdep(i).ge.fmin  .and. fdep(i).le.fmax)then
                                    inbound=1
                                    exit
                                endif
                            enddo
                            if(inbound.eq.0)then                            !!RECONSIDER - simply not setting when out of bounds
                                if(fdep(i).lt.params%poly(nset,1,1))then
                                      icond=1
                                      !print*, "node below lower bound of condition coverage"
                                else if(fdep(i).gt.params%poly(nset,params%f_polyseg(nset),2))then
                                      icond=params%f_polyseg(nset)
                                      !print*, "node above upper bound of condition coverage"
                                else
                                      stop "fatally screwed up polynome input - this can be the case if defined segments do not connect"
                                endif
                            endif
                            fwk=fdep(i)-foffs
        
                            d3=0
                            do j=1,params%pn(nset,icond)  ! loop and sum up contributions of all polynomial degrees
                                d3=d3  + params%poly(nset,icond,j+3) * (fwk)**(j-1)
                            enddo
                            select case(params%f_superpose(nset))
                                case(0)
                                    depvar(i)=d3
                                case(1)
                                    depvar(i)=depvar(i)+d3
                                case(2)
                                    depvar(i)=depvar(i)*d3
                            end select
                        enddo
                endif
        
        !!! @YW debug for multiple setting blocks when restart
                !if (params%read_restart)then
                !   print*,
                !   print*,
                !   print*,  'settting the ', nset, 'block'
                    
                !endif
        
        !!! @YW debug for multiple setting blocks when restart !!!!!!!!!!!!!!!!!!
        
            ! dissociate pointers for next condition
            nullify(fdep)
            nullify(depvar)
          enddo
  endif
end subroutine set_params_as_function





