!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
!
! United States and French Government Sponsorship Acknowledged.

  subroutine finalize_simulation()

  use adios_manager_mod
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use pml_par
  use gravity_perturbation, only: gravity_output, GRAVITY_SIMULATION

  implicit none

  integer :: ier
  !1.lucas 
  ! write gravity perturbations
  if (GRAVITY_SIMULATION) call gravity_output()

  ! save last frame
  !2.lucas,
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call save_forward_arrays_adios()
    else
      open(unit=IOUT,file=prname(1:len_trim(prname))//'save_forward_arrays.bin', &
            status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'error: opening save_forward_arrays.bin'
        print *,'path: ',prname(1:len_trim(prname))//'save_forward_arrays.bin'
        call exit_mpi(myrank,'error opening file save_forward_arrays.bin')
      endif

      if (ACOUSTIC_SIMULATION) then
        write(IOUT) potential_acoustic
        write(IOUT) potential_dot_acoustic
        write(IOUT) potential_dot_dot_acoustic
      endif

      if (ELASTIC_SIMULATION) then
        write(IOUT) displ
        write(IOUT) veloc
        write(IOUT) accel

        if (ATTENUATION) then
          write(IOUT) R_trace
          write(IOUT) R_xx
          write(IOUT) R_yy
          write(IOUT) R_xy
          write(IOUT) R_xz
          write(IOUT) R_yz
          write(IOUT) epsilondev_trace
          write(IOUT) epsilondev_xx
          write(IOUT) epsilondev_yy
          write(IOUT) epsilondev_xy
          write(IOUT) epsilondev_xz
          write(IOUT) epsilondev_yz
        endif
      endif

      if (POROELASTIC_SIMULATION) then
        write(IOUT) displs_poroelastic
        write(IOUT) velocs_poroelastic
        write(IOUT) accels_poroelastic
        write(IOUT) displw_poroelastic
        write(IOUT) velocw_poroelastic
        write(IOUT) accelw_poroelastic
      endif

      close(IOUT)
      if(CTD_SEM) then !-------------------------------------------m2
      open(unit=IOUT+100,file=prname(1:len_trim(prname))//'save_forward_arrays_m2.bin', &
            status='unknown',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'error: opening save_forward_arrays_m2.bin'
        print *,'path: ',prname(1:len_trim(prname))//'save_forward_arrays_m2.bin'
        call exit_mpi(myrank,'error opening file save_forward_arrays_m2.bin')
      endif
      if (ELASTIC_SIMULATION) then
        write(IOUT+100) displ_m2
        write(IOUT+100) veloc_m2
        write(IOUT+100) accel_m2

        if (ATTENUATION) then
          write(IOUT+100) R_trace_m2
          write(IOUT+100) R_xx_m2
          write(IOUT+100) R_yy_m2
          write(IOUT+100) R_xy_m2
          write(IOUT+100) R_xz_m2
          write(IOUT+100) R_yz_m2
          write(IOUT+100) epsilondev_trace_m2
          write(IOUT+100) epsilondev_xx_m2
          write(IOUT+100) epsilondev_yy_m2
          write(IOUT+100) epsilondev_xy_m2
          write(IOUT+100) epsilondev_xz_m2
          write(IOUT+100) epsilondev_yz_m2
        endif
      endif
      close(IOUT+100)
      endif !------------------------------------------------------m2


    endif
  endif
  
  !3.lucas 
  ! adjoint simulations
  if (SIMULATION_TYPE == 3) then
    ! adjoint kernels
    call save_adjoint_kernels() !CTD_SEM
  endif
  
  !4.lucas 
  ! seismograms and source parameter gradients for (pure type=2) adjoint simulation runs
  if (SIMULATION_TYPE == 2) then
    if (nrec_local > 0) then
      ! seismograms (strain)
      call write_adj_seismograms2_to_file(myrank,seismograms_eps,number_receiver_global,nrec_local,it,DT,NSTEP,t0)
      ! source gradients  (for sources in elastic domains)
      call save_kernels_source_derivatives()
    endif
  endif

  !5.lucas 
  ! stacey absorbing fields will be reconstructed for adjoint simulations
  ! using snapshot files of wavefields
  if (STACEY_ABSORBING_CONDITIONS) then
    ! closes absorbing wavefield saved/to-be-saved by forward simulations
    if (num_abs_boundary_faces > 0 .and. (SIMULATION_TYPE == 3 .or. &
          (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

      if (ELASTIC_SIMULATION) call close_file_abs(IOABS)
      if (ELASTIC_SIMULATION .and. CTD_SEM) call close_file_abs(IOABS_m2) ! CTD_SEM
       
      if (ACOUSTIC_SIMULATION) call close_file_abs(IOABS_AC)

    endif
  endif
  
  !6.lucas
  ! mass matrices
  if (ELASTIC_SIMULATION) then
    deallocate(rmassx)
    deallocate(rmassy)
    deallocate(rmassz)
    if(CTD_SEM) then !------------------------------------------m2
    deallocate(rmassx_m2)
    deallocate(rmassy_m2)
    deallocate(rmassz_m2)
    if(compute_approx_Hessian) then !----m1
    deallocate(rmassx_m1)
    deallocate(rmassy_m1)
    deallocate(rmassz_m1)
    endif !------------------------------m1
    endif !-----------------------------------------------------m2
  endif
  if (ACOUSTIC_SIMULATION) then
    deallocate(rmass_acoustic)
  endif
  
  !7.lucas 
  ! C-PML absorbing boundary conditions
  if (PML_CONDITIONS) then
    ! outputs informations about C-PML elements in VTK-file format
    if (NSPEC_CPML > 0) call pml_output_VTKs()
    ! deallocates C_PML arrays
    call pml_cleanup()
  endif
  
  !8.lucas 
  ! boundary surfaces
  deallocate(ibelm_xmin)
  deallocate(ibelm_xmax)
  deallocate(ibelm_ymin)
  deallocate(ibelm_ymax)
  deallocate(ibelm_bottom)
  deallocate(ibelm_top)
  
  !9.lucas 
  ! ADIOS file i/o
  if (ADIOS_ENABLED) then
    call adios_cleanup()
  endif
  
  !10.lucas 
  ! asdf finalizes
  if ((SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) .and. READ_ADJSRC_ASDF) then
    call asdf_cleanup()
  endif
  
  !11.lucas 
  ! close the main output file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of the simulation'
    write(IMAIN,*)
    close(IMAIN)
  endif

  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  end subroutine finalize_simulation
