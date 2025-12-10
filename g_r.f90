program dist_radios
    
    implicit none

    integer :: N_frame,frame,N, i, j, k, bin, nbins
    real(kind(8)) :: L,dr,pi      
    real(kind(8)), allocatable :: r(:,:), r_dist(:), g_sum(:), g(:)
    integer, allocatable :: acumulado(:)
    real(kind(8)) :: dxyz(3), rij, volumen, rho, n_id,r_max
    character(len=2) :: sym
    character(len=5) :: f
    open(10, file='punto_d/produccion.xyz', status='old')
    N_frame=1100
    L=8.0
    dr=0.01
    pi = acos(-1.0)
    volumen = L**3
    r_max = L/2
    nbins = int(L/dr/2)
    read(10,*) N
    read(10,*) f 
    print*, N
    print*, f

    allocate(r(N,3))
    allocate(acumulado(nbins), r_dist(nbins), g(nbins), g_sum(nbins))
    g_sum=0 
    do k=1,N_frame
        rho   = N / volumen
        !frame=1
    print*, k*N+1 , (k+1)*N
        
        do i = 1, N
            read(10,*) sym, r(i,1), r(i,2), r(i,3)
            !print*, r(i,1), r(i,2), r(i,3)
        end do
    
        acumulado = 0
        do i = 1, N-1
            do j = i+1, N
                dxyz = r(i,:) - r(j,:)
                dxyz = dxyz - L * nint(dxyz / L)
                rij = sqrt(sum(dxyz**2))
                
                !print*, rij 
                if (rij < r_max) then
                    bin = int(rij/dr) + 1
    !pongo 2 porque hay que contar dos veces una dist para la particula i y otra para la j
                    if (bin>1 .and. bin <= nbins) then  
                        acumulado(bin) = acumulado(bin) + 2
                    end if
                end if 
            end do
        end do
    !print*, 'suma de acumulado', sum(acumulado/N)
        do i = 1, nbins
            r_dist(i) = (i-0.5)*dr
            n_id = 4.0/3.0*pi*rho*((r_dist(i)+dr)**3-r_dist(i)**3)
            g(i) = acumulado(i) / n_id/N
        end do
        
        g_sum = g_sum + g
        read(10,*) N
        read(10,*) f 
        print*, N
        print*, f
    end do
    g_sum = g_sum / N_frame
    !guardo
        open(20, file='punto_d/g_r.dat', status='replace', action='write')
        do i = 1, nbins
            write(20,*) r_dist(i), g_sum(i)
        end do
        close(20)
        close(10)

end program dist_radios
