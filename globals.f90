module globals
        real(kind=8), allocatable :: r(:,:),v(:,:),F(:,:)
        real(kind=8) :: L,sigma=1.0,eps=1.0, E_tot, r_cut, v_cut, dt, m=1.0, T, rho, lgv_gam=0.5,T_inst,p_inst
        integer :: N

        
        real(kind=8) :: sigma6,c24eps,L2_4,frac_L,Ldiv2,rcut2,ex_cut,F_cut
end module globals

