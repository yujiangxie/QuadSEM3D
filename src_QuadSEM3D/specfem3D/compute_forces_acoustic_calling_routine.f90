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

! acoustic solver

! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
!
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then:
!     u = grad(Chi) / rho
! Velocity is then:
!     v = grad(Chi_dot) / rho
! (Chi_dot being the time derivative of Chi)
! and pressure is:
!     p = - Chi_dot_dot
! (Chi_dot_dot being the time second derivative of Chi).
!
! The source in an acoustic element is a pressure source.
!
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore Chi_dot_dot
! is continuous, therefore Chi is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u = grad(Chi) would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions. ! lucas, to be exactly undstood.

subroutine compute_forces_acoustic_calling()! lucas, should have simulation_type=2? 

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  use pml_par, only: is_CPML,spec_to_CPML,nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential, &
                     PML_potential_acoustic_old,PML_potential_acoustic_new

  implicit none

  ! local parameters
  integer:: iface,ispec,iglob,igll,i,j,k,ispec_CPML
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements (on MPI interface),
  !         iphase = 2 is for computing inner elements
  integer :: iphase
  !1.lucas
  !enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        BOTTOM_FREE_SURFACE,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)!lucas, to set output: potential_acoustic(iglob)=0,potential_dot_acoustic(iglob)=0,potential_dot_dot_acoustic(iglob)=0, 
                                                                 !lucas, where iglob = ibool(i,j,k,ispec), and ispec = free_surface_ispec(iface).

  if (USE_LDDRK) then
    call acoustic_enforce_free_surface_lddrk(NSPEC_AB,NGLOB_AB_LDDRK,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        BOTTOM_FREE_SURFACE,potential_acoustic_lddrk,potential_dot_acoustic_lddrk, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)! lucas, only set output: potential_acoustic_lddrk(iglob)=0, and potential_dot_acoustic_lddrk(iglob)=0. and not involving the dot_dot term.
  endif


  ! distinguishes two runs: for elements in contact with MPI interfaces, and elements within the partitions
  do iphase = 1,2 ! lucas ----------------------------iphase--------------------------------------------------------

    ! acoustic pressure term
    ! no need to test the two others because NGLLX == NGLLY = NGLLZ in unstructured meshes
    !2.lucas
    if (NGLLX == 5 .or. NGLLX == 6 .or. NGLLX == 7) then
     ! call compute_forces_acoustic_fast_Deville(iphase, potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,.false.) ! lucas, in:iphase, backward_simulation=false
       call compute_forces_acoustic_NGLL5_fast(iphase, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        .false.) ! lucas, move from specfem3D_original.
    else
      call compute_forces_acoustic_generic_slow(iphase, &
                        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        .false.)
    endif





    ! computes additional contributions
    if (iphase == 1) then ! lucas ---------------------------iphase =1-----//////////////////////////////////////////////////////////////////////////
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
        call compute_stacey_acoustic(NSPEC_AB,NGLOB_AB, &
                          potential_dot_dot_acoustic,potential_dot_acoustic, &
                          ibool,iphase, &
                          abs_boundary_jacobian2Dw,abs_boundary_ijk,abs_boundary_ispec, &
                          num_abs_boundary_faces,rhostore,kappastore,ispec_is_acoustic, &
                          SIMULATION_TYPE,SAVE_FORWARD,it,b_reclen_potential, &
                          b_absorb_potential,b_num_abs_boundary_faces)!lucas, to update potential_dot_dot_acoustic, which is inout, 
                                                                      !lucas, to get b_absorb_potential, which is out, 
                                                                      !lucas, to write b_absorb_potential in case of adjoint simulations.

                                                                      !lucas, the other variables are input.
      endif

      ! elastic coupling, lucas -----------------------------------------------------------------
      if (ELASTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          if (SIMULATION_TYPE == 1) then
            ! forward definition: \bfs=\frac{1}{\rho}\bfnabla\phi
            call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                                ibool,displ,potential_dot_dot_acoustic, &
                                num_coupling_ac_el_faces, &
                                coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                coupling_ac_el_normal, &
                                coupling_ac_el_jacobian2Dw, &
                                iphase, &
                                PML_CONDITIONS,SIMULATION_TYPE,.false.)!lucas, to update(inout) potential_dot_dot_acoustic, the others all are input

          else !lucas, SIMULATION_TYPE /= 1
            ! handles adjoint runs coupling between adjoint potential and adjoint elastic wavefield
            ! adjoint definition: \partial_t^2 \bfs^\dagger=-\frac{1}{\rho}\bfnabla\phi^\dagger
            call compute_coupling_acoustic_el(NSPEC_AB,NGLOB_AB, &
                                ibool,accel,potential_dot_dot_acoustic, &
                                num_coupling_ac_el_faces, &
                                coupling_ac_el_ispec,coupling_ac_el_ijk, &
                                coupling_ac_el_normal, &
                                coupling_ac_el_jacobian2Dw, &
                                iphase, &
                                PML_CONDITIONS,SIMULATION_TYPE,.false.)! lucas, here is accel, the bove is displ. 
                                                                       ! lucas, in this adjoint case, just set the displ_x = -displ(1,iglob)[!where displ=accel here], also for the y and z.
                                                                       ! lucas, need to unstand the differences between forward definition and adjoint definition, which is interesting. todo, 04/2019
          endif
        endif
      endif ! lucas--------------------------------------------------------------------------------

      ! poroelastic coupling, lucas--------------------------
      if (POROELASTIC_SIMULATION ) then
        if (num_coupling_ac_po_faces > 0) then
          if (SIMULATION_TYPE == 1) then
            call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                          ibool,displs_poroelastic,displw_poroelastic, &
                          potential_dot_dot_acoustic, &
                          num_coupling_ac_po_faces, &
                          coupling_ac_po_ispec,coupling_ac_po_ijk, &
                          coupling_ac_po_normal, &
                          coupling_ac_po_jacobian2Dw, &
                          iphase)!lucas, to update(inout) potential_dot_dot_acoustic, the others are input
          else
            stop 'not implemented yet'
          endif
        endif
      endif !lucas -------------------------------------------

      ! adds sources
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      call compute_add_sources_acoustic()! lucas, to update potential_dot_dot_acoustic

    endif ! iphase---------------------------------end iphase1-----//////////////////////////////////////////////////////////////////////////////////




    ! assemble all the contributions between slices using MPI
    !3.lucas
    if (iphase == 1) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      !lucas, sum the 'potential_dot_dot_acoustic' of each ibool of each iinterface togather,
      call assemble_MPI_scalar_async_send(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        my_neighbors_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)!lucas, array to send is potential_dot_dot_acoustic
    else !lucas, iphase=2

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_async_recv(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)
    endif
  enddo! lucas ----------------------------------------------------------end iphase i=1 and 2----------------------------------------


  !4.lucas, update potential_dot_dot_acoustic determined above, 
  !divides pressure with mass matrix
  potential_dot_dot_acoustic(:) = potential_dot_dot_acoustic(:) * rmass_acoustic(:) ! lucas, rmass_acoustic is obtained in create_mass_matrices.f90

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

!! DK DK this paragraph seems to be from Zhinan or from ChangHua:
! However, enforcing explicitly potential_dot_dot_acoustic, potential_dot_acoustic, potential_acoustic
! to be zero on outer boundary of PML help to improve the accuracy of absorbing low-frequency wave components
! in case of long-time simulation.
  !5.lucas, set p, dp/dt, and (dp)^2/dt^2 of the boundary equal to 0.
  !impose Dirichlet conditions for the potential (i.e. Neumann for displacement) on the outer edges of the C-PML layers
  if (PML_CONDITIONS .and. SET_NEUMANN_RATHER_THAN_DIRICHLET_FOR_FLUID_PMLs) then
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
!!! It is better to move this into do iphase=1,2 loop
        if (ispec_is_acoustic(ispec) .and. is_CPML(ispec)) then
          ! reference GLL points on boundary face
          ispec_CPML = spec_to_CPML(ispec)
          do igll = 1,NGLLSQUARE
            ! gets local indices for GLL point
            i = abs_boundary_ijk(1,igll,iface)
            j = abs_boundary_ijk(2,igll,iface)
            k = abs_boundary_ijk(3,igll,iface)

            iglob=ibool(i,j,k,ispec)

            potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
            potential_dot_acoustic(iglob) = 0._CUSTOM_REAL
            potential_acoustic(iglob) = 0._CUSTOM_REAL
            if (ELASTIC_SIMULATION) then
              PML_potential_acoustic_old(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
              PML_potential_acoustic_new(i,j,k,ispec_CPML) = 0._CUSTOM_REAL
            endif
          enddo
        endif ! ispec_is_acoustic
!!!   endif
    enddo
  endif

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 DELTA_T CHI_DOT_DOT( T + DELTA_T )
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the chi_dot term which requires chi_dot_dot(t+delta)
  ! corrector
  !6.lucas
  if (USE_LDDRK) then
    call update_potential_dot_acoustic_lddrk() !lucas, ! updates acceleration, velocity and displacement in acoustic region (outer core), i.e., p, dp/dt, (dp)^2/dt^2
  else
    potential_dot_acoustic(:) = potential_dot_acoustic(:) + deltatover2*potential_dot_dot_acoustic(:)
  endif

  !7.enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        BOTTOM_FREE_SURFACE,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)

  if (USE_LDDRK) then
    call acoustic_enforce_free_surface_lddrk(NSPEC_AB,NGLOB_AB_LDDRK,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        BOTTOM_FREE_SURFACE,potential_acoustic_lddrk,potential_dot_acoustic_lddrk, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)
  endif
  !8.lucas, to get potential_acoustic_adj_coupling for simulation_type/=1.
  if (SIMULATION_TYPE /= 1) then
    potential_acoustic_adj_coupling(:) = potential_acoustic(:) &
                            + deltat * potential_dot_acoustic(:) &
                            + deltatsqover2 * potential_dot_dot_acoustic(:)
  endif
  !8.lucas, save p, dp/dt, and (dp)^2/dt^2 for the each GLL points of the interfaces.
  if (PML_CONDITIONS) then
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
      if (nglob_interface_PML_acoustic > 0) then
        !lucas, save b_PML_potential(1,iglob)=potential_acoustic, b_PML_potential(2,iglob)=potential_dot_acoustic, b_PML_potential(3,iglob)=potential_dot_dot_acoustic, for backforward simulation
        !lucas, where iglob = 1, nglob_interface_PML_acoustic
        call save_potential_on_pml_interface(potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                             nglob_interface_PML_acoustic,b_PML_potential,b_reclen_PML_potential)!lucas,   out: b_PML_potential, the others are input. to be check the potential_acoustic here =0 or not in case of free surface, it should be 0 based on my knowledge right now, 04/2019.???????????????????????????
      endif
    endif
  endif

end subroutine compute_forces_acoustic_calling

!
!=====================================================================

! acoustic solver

! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
!
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then:
!     u = grad(Chi) / rho
! Velocity is then:
!     v = grad(Chi_dot) / rho
! (Chi_dot being the time derivative of Chi)
! and pressure is:
!     p = - Chi_dot_dot
! (Chi_dot_dot being the time second derivative of Chi).
!
! The source in an acoustic element is a pressure source.
!
! First-order acoustic-acoustic discontinuities are also handled automatically
! because pressure is continuous at such an interface, therefore Chi_dot_dot
! is continuous, therefore Chi is also continuous, which is consistent with
! the spectral-element basis functions and with the assembling process.
! This is the reason why a simple displacement potential u = grad(Chi) would
! not work because it would be discontinuous at such an interface and would
! therefore not be consistent with the basis functions.


subroutine compute_forces_acoustic_backward_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: iphase

  ! checks
  if (SIMULATION_TYPE /= 3) &
    call exit_MPI(myrank,'error calling compute_forces_acoustic_backward() with wrong SIMULATION_TYPE')

  ! adjoint simulations
  !1.lucas
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                      BOTTOM_FREE_SURFACE,b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      ibool,free_surface_ijk,free_surface_ispec, &
                      num_free_surface_faces,ispec_is_acoustic)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! adjoint simulations
    !2.lucas
    ! no need to test the two others because NGLLX == NGLLY = NGLLZ in unstructured meshes
    if (NGLLX == 5 .or. NGLLX == 6 .or. NGLLX == 7) then
      call compute_forces_acoustic_fast_Deville(iphase, &
                      b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      .true.)
    else
      call compute_forces_acoustic_generic_slow(iphase, &
                      b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      .true.)
    endif

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
         call compute_stacey_acoustic_backward(NSPEC_AB, &
                              ibool,iphase, &
                              abs_boundary_ijk,abs_boundary_ispec, &
                              num_abs_boundary_faces,ispec_is_acoustic, &
                              SIMULATION_TYPE,NSTEP,it,NGLOB_ADJOINT, &
                              b_potential_dot_dot_acoustic,b_reclen_potential, &
                              b_absorb_potential,b_num_abs_boundary_faces)
      endif

      ! elastic coupling
      if (ELASTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          ! adjoint/kernel simulations
          call compute_coupling_acoustic_el(NSPEC_ADJOINT,NGLOB_ADJOINT, &
                            ibool,b_displ,b_potential_dot_dot_acoustic, &
                            num_coupling_ac_el_faces, &
                            coupling_ac_el_ispec,coupling_ac_el_ijk, &
                            coupling_ac_el_normal, &
                            coupling_ac_el_jacobian2Dw, &
                            iphase, &
                            PML_CONDITIONS,SIMULATION_TYPE,.true.)
        endif
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION ) then
        if (num_coupling_ac_po_faces > 0) then
            stop 'coupling acoustic-poroelastic domains not implemented yet...'
        endif
      endif

      ! adds sources
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      call compute_add_sources_acoustic_backward()

    endif ! iphase
    !3.lucas
    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends b_potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      ! adjoint simulations
      call assemble_MPI_scalar_async_send(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                      b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                      num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                      my_neighbors_ext_mesh, &
                      b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

    else
      ! adjoint simulations
      call assemble_MPI_scalar_async_recv(NPROC,NGLOB_ADJOINT,b_potential_dot_dot_acoustic, &
                      b_buffer_recv_scalar_ext_mesh,num_interfaces_ext_mesh, &
                      max_nibool_interfaces_ext_mesh, &
                      nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                      b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)
    endif
  enddo

  ! divides pressure with mass matrix
  ! adjoint simulations
  ! 4.lucas
  b_potential_dot_dot_acoustic(:) = b_potential_dot_dot_acoustic(:) * rmass_acoustic(:)

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 DELTA_T CHI_DOT_DOT( T + DELTA_T )
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
!   updates the chi_dot term which requires chi_dot_dot(t+delta)

  ! corrector
  ! adjoint simulations
  !5.lucas
  b_potential_dot_acoustic(:) = b_potential_dot_acoustic(:) + b_deltatover2*b_potential_dot_dot_acoustic(:)

! enforces free surface (zeroes potentials at free surface)
  ! adjoint simulations
  !6.lucas
  call acoustic_enforce_free_surface(NSPEC_AB,NGLOB_ADJOINT,STACEY_INSTEAD_OF_FREE_SURFACE, &
                      BOTTOM_FREE_SURFACE,b_potential_acoustic,b_potential_dot_acoustic,b_potential_dot_dot_acoustic, &
                      ibool,free_surface_ijk,free_surface_ispec, &
                      num_free_surface_faces,ispec_is_acoustic)

end subroutine compute_forces_acoustic_backward_calling

!
!-------------------------------------------------------------------------------------------------
!

subroutine compute_forces_acoustic_GPU_calling()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  ! local parameters
  integer:: iphase

  ! check
  if (PML_CONDITIONS) call exit_MPI(myrank,'PML conditions not yet implemented on GPUs')

  ! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE)

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! acoustic pressure term
    ! includes code for SIMULATION_TYPE==3
    call compute_forces_acoustic_cuda(Mesh_pointer, iphase, &
                                      nspec_outer_acoustic, nspec_inner_acoustic)

    ! computes additional contributions
    if (iphase == 1) then
      ! Stacey absorbing boundary conditions
      if (STACEY_ABSORBING_CONDITIONS) then
         call compute_stacey_acoustic_GPU(iphase,num_abs_boundary_faces, &
                              SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
                              b_reclen_potential,b_absorb_potential, &
                              b_num_abs_boundary_faces,Mesh_pointer)
      endif

      ! elastic coupling
      if (ELASTIC_SIMULATION) then
        if (num_coupling_ac_el_faces > 0) then
          ! on GPU
          call compute_coupling_ac_el_cuda(Mesh_pointer,iphase, &
                                           num_coupling_ac_el_faces)
        endif
      endif

      ! poroelastic coupling
      if (POROELASTIC_SIMULATION ) then
        if (num_coupling_ac_po_faces > 0) then
          if (SIMULATION_TYPE == 1) then
            call compute_coupling_acoustic_po(NSPEC_AB,NGLOB_AB, &
                          ibool,displs_poroelastic,displw_poroelastic, &
                          potential_dot_dot_acoustic, &
                          num_coupling_ac_po_faces, &
                          coupling_ac_po_ispec,coupling_ac_po_ijk, &
                          coupling_ac_po_normal, &
                          coupling_ac_po_jacobian2Dw, &
                          iphase)
          else
            stop 'not implemented yet'
          endif
          if (SIMULATION_TYPE == 3) &
            stop 'not implemented yet'
        endif
      endif

      ! adds sources
      ! note: we will add all source contributions in the first pass, when iphase == 1
      !       to avoid calling the same routine twice and to check if the source element is an inner/outer element
      !
      call compute_add_sources_acoustic_GPU()

    endif ! iphase

    ! assemble all the contributions between slices using MPI
    if (iphase == 1) then
      ! sends potential_dot_dot_acoustic values to corresponding MPI interface neighbors (non-blocking)
      call transfer_boun_pot_from_device(Mesh_pointer, &
                                         potential_dot_dot_acoustic, &
                                         buffer_send_scalar_ext_mesh, &
                                         1) ! -- 1 == fwd accel
      call assemble_MPI_scalar_send_cuda(NPROC, &
                        buffer_send_scalar_ext_mesh,buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                        nibool_interfaces_ext_mesh, &
                        my_neighbors_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call transfer_boun_pot_from_device(Mesh_pointer, &
                                           b_potential_dot_dot_acoustic, &
                                           b_buffer_send_scalar_ext_mesh, &
                                           3) ! -- 3 == adjoint b_accel

        call assemble_MPI_scalar_send_cuda(NPROC, &
                          b_buffer_send_scalar_ext_mesh,b_buffer_recv_scalar_ext_mesh, &
                          num_interfaces_ext_mesh,max_nibool_interfaces_ext_mesh, &
                          nibool_interfaces_ext_mesh, &
                          my_neighbors_ext_mesh, &
                          b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh)

      endif

    else

      ! waits for send/receive requests to be completed and assembles values
      call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,potential_dot_dot_acoustic, &
                        Mesh_pointer, &
                        buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        request_send_scalar_ext_mesh,request_recv_scalar_ext_mesh, &
                        1)

      ! adjoint simulations
      if (SIMULATION_TYPE == 3) then
        call assemble_MPI_scalar_write_cuda(NPROC,NGLOB_AB,b_potential_dot_dot_acoustic, &
                        Mesh_pointer, &
                        b_buffer_recv_scalar_ext_mesh, &
                        num_interfaces_ext_mesh, &
                        max_nibool_interfaces_ext_mesh, &
!                       nibool_interfaces_ext_mesh,ibool_interfaces_ext_mesh, &
                        b_request_send_scalar_ext_mesh,b_request_recv_scalar_ext_mesh, &
                        3)
      endif
    endif
  enddo

! update velocity
! note: Newmark finite-difference time scheme with acoustic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 DELTA_T
! CHI_DOT_DOT( T + DELTA_T )
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic
! u(t+delta_t) + f(t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with
!   respect to time)
!   M is mass matrix, K stiffness matrix and B boundary term
!   f denotes a source term
!
! corrector:
! updates the chi_dot term which requires chi_dot_dot(t+delta)

  call kernel_3_acoustic_cuda(Mesh_pointer,deltatover2,b_deltatover2)

! enforces free surface (zeroes potentials at free surface)
  call acoustic_enforce_free_surf_cuda(Mesh_pointer,STACEY_INSTEAD_OF_FREE_SURFACE)

end subroutine compute_forces_acoustic_GPU_calling

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface(NSPEC_AB,NGLOB_AB,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        BOTTOM_FREE_SURFACE,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)! lucas, set potential_acoustic=0,potential_dot_acoustic=0,potential_dot_dot_acoustic=0

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB
  logical :: STACEY_INSTEAD_OF_FREE_SURFACE, BOTTOM_FREE_SURFACE

! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB) :: &
        potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! free surface
  integer :: num_free_surface_faces
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! checks if free surface became an absorbing boundary
  if (STACEY_INSTEAD_OF_FREE_SURFACE .and. .not. BOTTOM_FREE_SURFACE) return

! enforce potentials to be zero at surface
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)

        ! sets potentials to zero
        potential_acoustic(iglob)         = 0._CUSTOM_REAL
        potential_dot_acoustic(iglob)     = 0._CUSTOM_REAL
        potential_dot_dot_acoustic(iglob) = 0._CUSTOM_REAL
      enddo
    endif

  enddo

end subroutine acoustic_enforce_free_surface

!
!-------------------------------------------------------------------------------------------------
!

subroutine acoustic_enforce_free_surface_lddrk(NSPEC_AB,NGLOB_AB_LDDRK,STACEY_INSTEAD_OF_FREE_SURFACE, &
                        BOTTOM_FREE_SURFACE, potential_acoustic_lddrk,potential_dot_acoustic_lddrk, &
                        ibool,free_surface_ijk,free_surface_ispec, &
                        num_free_surface_faces,ispec_is_acoustic)! lucas, potential_acoustic_lddrk,potential_dot_acoustic_lddrk

  use constants

  implicit none

  integer :: NSPEC_AB,NGLOB_AB_LDDRK
  logical :: STACEY_INSTEAD_OF_FREE_SURFACE,BOTTOM_FREE_SURFACE

! acoustic potentials
  real(kind=CUSTOM_REAL), dimension(NGLOB_AB_LDDRK) :: &
            potential_acoustic_lddrk,potential_dot_acoustic_lddrk

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_AB) :: ibool

! free surface
  integer :: num_free_surface_faces
  integer :: free_surface_ijk(3,NGLLSQUARE,num_free_surface_faces)
  integer :: free_surface_ispec(num_free_surface_faces)

  logical, dimension(NSPEC_AB) :: ispec_is_acoustic

! local parameters
  integer :: iface,igll,i,j,k,ispec,iglob

  ! checks if free surface became an absorbing boundary
  if (STACEY_INSTEAD_OF_FREE_SURFACE .and. .not. BOTTOM_FREE_SURFACE) return

! enforce potentials to be zero at surface
  do iface = 1, num_free_surface_faces

    ispec = free_surface_ispec(iface)

    if (ispec_is_acoustic(ispec)) then

      do igll = 1, NGLLSQUARE
        i = free_surface_ijk(1,igll,iface)
        j = free_surface_ijk(2,igll,iface)
        k = free_surface_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)

        ! sets potentials to zero
        potential_acoustic_lddrk(iglob)         = 0._CUSTOM_REAL
        potential_dot_acoustic_lddrk(iglob)     = 0._CUSTOM_REAL
      enddo
    endif

  enddo

end subroutine acoustic_enforce_free_surface_lddrk

