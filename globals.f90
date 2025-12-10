module globals
        real(kind=8), allocatable :: r(:,:),v(:,:),F(:,:)
        real(kind=8) :: L,sigma=1.0,eps=1.0, E_tot, r_cut, v_cut, dt, m=1.0, T, rho, lgv_gam=0.5,T_inst,p_inst
        integer :: N
end module globals

