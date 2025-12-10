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
        integer :: i,j
        real(kind=8) :: r_ij(3), norm, ex, v_ij, p
        p = 0.0 
        E_tot = 0.0
        F = 0.0
        do i=1,N-1
                do j=i+1,N
                        r_ij = r(:,i)-r(:,j)
                        norm = sqrt(sum(r_ij**2))

                        !corregimos distancia a la imagen mas cercana
                        if(norm > L/2) then
                                r_ij = r_ij - (L * int(2*r_ij/L))
                                norm = sqrt(sum(r_ij**2))
                        end if

                        !computamos fuerzas y energia solo si la distancia es menor a r_cut
                        if(norm < r_cut) then
                                ex = (sigma/norm)**6
                                v_ij = 4.0*eps*(-ex+(ex*ex)) - v_cut
                                E_tot = E_tot + v_ij
                                
                                ex_cut = (sigma/r_cut)**6
                                F_cut = 24.0*eps*(ex_cut-2.0*ex_cut*ex_cut)/(r_cut*r_cut)
                                F_esc = 24.0*eps*(ex-2.0*ex*ex)/(norm*norm) - F_cut
                                F(:,i) = F(:,i) - F_esc*r_ij
                                F(:,j) = F(:,j) + F_esc*r_ij
                                !F(:,i) = F(:,i)-(r(:,i)*(4.0*eps*(6*ex-12*ex*ex))/norm**2)
                                !F(:,j) = F(:,j)-(r(:,j)*(4.0*eps*(6*ex-12*ex*ex))/norm**2)
                                p = p + dot_product(r_ij,-F_esc*r_ij)
                        end if
                end do
        end do
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
        r = r + v*dt + 0.5*(F/m)*dt**2
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


