module MD

use globals
use ziggurat

contains

subroutine init_coords()
        integer :: i,j
        
        do i=1,N
                do j=1,3
                        r(j,i) = uni()*L
                end do
        end do
end subroutine init_coords

subroutine update_E_and_F() 
    implicit none 
    integer :: i,j,k,exss
    real(kind=8) :: r_ij(3), norm2, ex, ex2, v_ij, p, F_esc
    real(kind=8) :: ex_cut, F_cut, delta, frac_L, Ldiv2
    real(kind=8) :: sigma6, c24eps, L2_4, rcut2, norm2_3

    ! Inicializaciones
    p = 0.0
    E_tot = 0.0
    F = 0.0

    ! Cálculos fuera del loop
    sigma6 = sigma**6
    c24eps = 24.0*eps
    L2_4   = L*L/4
    frac_L = 2.0/L
    Ldiv2 = L/2.0
    rcut2  = r_cut*r_cut
    ex_cut = (sigma/r_cut)**6
    F_cut  = c24eps*(ex_cut - 2.0*ex_cut*ex_cut)/rcut2

    !$OMP PARALLEL PRIVATE(i,j,k,exss,r_ij,norm2,ex,ex2,v_ij,delta,norm2_3,F_esc) REDUCTION(+:p, E_tot)
    ! Loop sobre pares de partículas
    !$OMP DO SCHEDULE(STATIC,2)
    do i=1,N-1
        do j=i+1,N
            ! Diferencia de posición
            !!r_ij = r(:,i) - r(:,j)
            !norm2 = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)

            ! Corrección de imágenes periódicas
            !if(norm2 > L2_4) then
            !    r_ij = r_ij - (L * int(2*r_ij/L))
            !    norm2 = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)
            !end if
            !!r_ij = r_ij - L * nint(2*r_ij/L)
            
            !se hace la correccion de imagen periodica elemento a elemento y luego se calcula norm2, para calcularlo una sola vez
            do k = 1,3
                delta = r(k,i) - r(k,j)
                ! Aplicar PBC solo si es necesario
                if(abs(delta) > Ldiv2) then
                    exss = nint(delta*frac_L)
                    delta = delta - exss*L
                end if
                r_ij(k) = delta
            end do
            norm2 = r_ij(1)*r_ij(1) + r_ij(2)*r_ij(2) + r_ij(3)*r_ij(3)

            ! Computar fuerzas y energía si dentro de r_cut
            if(norm2 < rcut2) then
                ! Evitar potencia usando multiplicaciones
                norm2_3 = norm2*norm2*norm2
                ex   = sigma6 / norm2_3
                ex2 = ex*ex
                v_ij = 4.0*eps*(-ex + ex2)
                v_ij = v_ij - v_cut
                E_tot = E_tot + v_ij

                F_esc = c24eps*(ex - 2.0*ex2)/norm2 - F_cut
                F(:,i) = F(:,i) - F_esc*r_ij
                F(:,j) = F(:,j) + F_esc*r_ij

                p = p  -F_esc * norm2
                !p = p + dot_product(r_ij, -F_esc*r_ij)
            end if
        end do
    end do
   !$OMP END DO
   !$OMP END PARALLEL
    p_inst = N*get_T()/(L**3) + 1.0/(3*L**3)*p
end subroutine update_E_and_F

subroutine update_lgv_F()
        real(kind=8) :: rand_dev

        rand_dev = sqrt(2.0*lgv_gam*m*T/dt)
        F = F - lgv_gam*v*m
        do i=1,N
                do j=1,3
                        F(j,i) = F(j,i) + rand_dev*rnor()
                end do
        end do
end subroutine update_lgv_F

subroutine pbc()
        r = r - L*floor(r/L)
end subroutine pbc       

function E_minimization(steps, stride) result(Es)
        integer, intent(in) :: steps, stride
        real(kind=8), allocatable :: Es(:)
        character :: fname(14)
        allocate(Es(int(steps/stride)))

        do i=1,steps
                r = r + 0.5*(F/m)*dt**2
                call pbc
                call update_E_and_F()

                if (MOD(i,stride) == 0.0) then
                        fname = "coords_min.xyz"
                        call save_coords("coords_min.xyz")
                        Es(int(i/stride)) = E_tot
                end if

        end do
end function E_minimization

subroutine save_coords(filename)
        character(len=*), intent(in) :: filename
        logical :: es
        
        inquire(file = filename, exist=es)
        if(es) then
                open(unit=1, file=filename, access="append", status="old", action="write")
        
        else
                open(unit=1, file=filename, status="new", action="write")
                write(1,*) N
        end if

        
        do i=1,N
                write(1,"(A1,F16.8,A1,F16.8,A1,F16.8)") "C       ",r(1,i),"      ",r(2,i),"      ",r(3,i)
        end do
        close(1)
end subroutine save_coords


subroutine integrate() 
        !asume que la fuerza del paso anterior ya fue calculada
        r = r + v*dt + 0.5*(F/m)*dt**2  !revisar si esto es paralelizable
        call pbc()
        v = v + 0.5*F*dt/m
        call update_E_and_F()
        call update_lgv_F()
        v = v + 0.5*F*dt/m
        T_inst = get_T()
end subroutine integrate

subroutine initiate_velocities()
        v = 0.0
        do i=1,N
                v(1,i) = rnor()*sqrt(T/m)
                v(2,i) = rnor()*sqrt(T/m)
                v(3,i) = rnor()*sqrt(T/m)
        end do
end subroutine initiate_velocities


function get_T() result(T)
        real(kind=8) :: T
        T = 1.0/(3.0*N-3.0) * sum(v*v)
end function get_T

function get_Ek() result(Ek)
        real(kind=8) :: Ek
        Ek = sum(0.5*(v*v))
end function get_Ek

end module MD


