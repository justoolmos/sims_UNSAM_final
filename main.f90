program main
        use globals
        use ziggurat
        use MD

        integer :: seed,i
        logical :: es
        real(kind=8), allocatable :: Es_min(:)
        real(kind=8) :: Ek, inst_T, p_in, p

        inquire(file='seed.dat',exist=es)
        if(es) then
            open(unit=10,file='seed.dat',status='old')
            read(10,*) seed
            close(10)
            print *,"  * Leyendo semilla de archivo seed.dat"
        else
            seed = 24583490
        end if

        call zigset(seed)

        open(unit=30, file="in.dat", status="old", action="read")
        read(30,*) T,rho,L,n_steps,dt
        close(30)
       
        N = rho*(L**3)
        print *, N
        r_cut = sigma*2.5
        v_cut = 4*eps*(-(sigma/r_cut)**6+(sigma/r_cut)**12)


        allocate(r(3,N))
        allocate(v(3,N))
        allocate(F(3,N))

        call init_coords()
        call update_E_and_F()
        call update_lgv_F()
        
        Es_min = E_minimization(5000,1)
        call initiate_velocities()
        T_inst = get_T()

        open(unit=10, file = "produccion.dat", action="write", status="replace", form="formatted")
        write(10,"(A34)") "step,E total,Ep,Ek,F max,v max,T,p"

        do i=1,n_steps
                call integrate()
                if(mod(i,100)==0) then
                        Ek = get_Ek()
                        write(10,"(I5,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8)") &
                        i,",",E_tot+Ek,",",E_tot,",",Ek,",",maxval(F),",",maxval(v),",",T_inst,",",p_inst
                        write(*,"(A6,I5)") "step: ", i
                        print *, E_tot + Ek 
                        call save_coords("produccion.xyz")
                end if
        end do
        close(10)
                
end program main







