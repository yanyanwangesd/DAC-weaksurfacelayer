subroutine find_erodibility (geometry,params)
! update nodal erodibility at regolith regime, specifically for the escarpment topogrpahy
! philosophy is the slow-eroding plateau has a finite-thickness surface layer of regolith
! the fast-eroding escarpment side is harder bedrock. @YWANG, Nov2021
! 


use definitions

implicit none

type (geom) geometry
type (parm) params


double precision ave_erosion, erosion_threshold
integer i
double precision twindow
integer num_samples, num_full_bin,sample_in_last_bin
integer cindex

call time_in ('find_erodibility')

! define values of thresholds

! Tstart = 90.d6 ! time to start make erodibility vary with erosion rate


! parameters to make K varies with erosion rate, not fully implemented yet
twindow = 5.d5 ! characteristic time scale of weathering
erosion_threshold = 5.d-6 ! characteristic threshold erosion rate to weathering-dominated

! Only start when model run longer than Tstart
if (params%time >= params%regolith_time ) then

    ! loop through nodes
    do i=1,geometry%nnode
        cindex = geometry%catchment(i)
        
        !! side boundary basins, for plateau rivers that also flow to the side basins
        !if(geometry%fix(cindex)==3 .or. geometry%fix(cindex)==4) then ! initialize regolith layer
        !     if(params%time .eq. Tstart)then
        !       geometry%k(i) = k_ratio*params%k_scalar1
        !        geometry%dregolith(i) = regolith_thickness
        !    else
        !        if(geometry%dregolith(i).gt.0.d0)then
        !           geometry%k(i) = k_ratio*params%k_scalar1
        !            geometry%dregolith(i) = geometry%dregolith(i)-geometry%erosion_rate(i)*params%deltat
        !            geometry%precipitation(i) = params%rainfall_height/(k_ratio**(1/params%m)) !scale rainfall to keep erosion rate unchanged
        !        else
        !            geometry%k(i) = params%k_scalar1
        !            geometry%dregolith(i) = -1.d0
        !        endif
        !    endif
        ! endif
        
        ! side boundary basins, for plateau basins that only flow to the right boundary
        if(geometry%fix(cindex)==3 .or. geometry%fix(cindex)==4) then
            geometry%k(i) = params%k_scalar1
            geometry%dregolith(i) = 0.d0
        endif
        
        ! plateau basins, only update K for catchment but don't change erosion rate for basins flow out of the right boundary
        ! do this by scale precipitation (discharge) accordingly
        if(geometry%fix(cindex)==2) then
            geometry%k(i) = params%k_ratio*params%k_scalar1
            geometry%dregolith(i) = params%regolith_thickness
            geometry%precipitation(i) = params%rainfall_height/(params%k_ratio**(1/params%m)) !scale rainfall to keep erosion rate unchanged
        endif
        
        ! escarpment basins, update K when regoligh layer is fully eroded away
        if(geometry%fix(cindex)==1) then
            if (geometry%dregolith(i)> 0) then
                geometry%dregolith(i) = geometry%dregolith(i)-geometry%erosion_rate(i)*params%deltat
                if (geometry%dregolith(i)>0) then
                    geometry%k(i) = params%k_ratio*params%k_scalar1
                else
                    geometry%k(i) = params%k_scalar1
                    geometry%dregolith(i) = 0
                endif
            endif
        endif
    enddo ! end lood of nodes
    
    !another option is to make K varies with the erosion rate history. Setting a critial erosion rate, 
    ! if the average erosion rate over a characteristic time window is smaller than the critical erosion rate,
    ! a layer of regolith is formed and it increases thickness until the limit, not fully implemented yet.
    
    
    
    
    
    
endif


call time_out ('find_erodibility')

return

end subroutine find_erodibility
