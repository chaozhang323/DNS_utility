module MGRIDGEN2D
use omp_lib
use cm_gridgen
use modMetrics
use modTecbin
implicit none
real(8),parameter,private :: real8_tol = 1.e-14
real(8),private :: err_inner = 1.e-8, err_outer = 1.e-6, pi = 3.14159265
integer,private :: nsave = 50, maxiter_inner = 150, maxiter_outer = 100, max_change = 0, nrefine = 4, nblocks = 4
integer,dimension(4),private :: boundary_symmetric = (/0,0,0,0/)
integer,dimension(4),private :: boundary_fixed = (/1,1,0,0/)

contains
    subroutine InitOrthogonalMesh()
        character(100) :: fn_control
        logical :: ex
        integer :: i

        fn_control = 'orthogonal_control.inp'
        inquire(file=trim(fn_control),exist=ex)
        if(ex) then
            open(33,file=trim(fn_control),status='old')
            read(33,*)
            read(33,*) (boundary_fixed(i),i=1,4)
            read(33,*)
            read(33,*) (boundary_symmetric(i),i=1,4)
            read(33,*)
            read(33,*) nrefine, nsave, maxiter_outer, maxiter_inner, max_change
            read(33,*)
            read(33,*) err_inner, err_outer
        else
            print*,"Control file 'orthogonal_control.inp' NOT exist."
            print*,'Default control parameters are used. '
        endif
        print*,'------------ Control parameters for orthogonal solver ------------'
        print*,'boundary fixed (left, right, bot, top): '
        print*,boundary_fixed
        print*,'boundary symmetric (left, right, bot, top): '
        print*,boundary_symmetric
        print*,'nrefine: ',nrefine
        print*,'nstdout: ',nsave
        print*,'maxiter_outer: ',maxiter_outer
        print*,'maxiter_inner: ',maxiter_inner
        print*,'max_change: ',max_change
        print*,'tol_outer: ',err_outer
        print*,'tol_inner: ',err_inner

        pi = 4.d0*atan(1.d0)
    end subroutine InitOrthogonalMesh

    !--------------------------------------------------------------------------------------------
    ! Orthogonal mesh generator with distortion function (weak constraint method)
    ! See the subroutine "compute_grid_new" for details.
    !--------------------------------------------------------------------------------------------
    subroutine GenOrthogonalMesh(rleft,rright,rbot,rtop,x2d,y2d)
        real(8),dimension(:,:),intent(in) :: rleft,rright,rbot,rtop
        real(8),dimension(:,:),intent(inout) :: x2d,y2d
        integer :: n,i,j,nx,ny,nloop
        logical :: iloop = .True.
        ! orthogonal mesh generation
        integer num_choice,nx_tmp1,nx_tmp2
        integer,dimension(:),allocatable :: idx_left,idx_right,idx_bot,idx_top
        real(8),dimension(:,:),allocatable :: choice_left,choice_right,choice_bot,choice_top,choice_bot_tmp,choice_top_tmp
        real(8),dimension(:,:),allocatable :: x2d_new,z2d_new
        ! dimension
        nx = size(rbot,dim=1)
        ny = size(rleft,dim=1)
        if(size(x2d,dim=1).ne.nx .or. size(x2d,dim=2).ne.ny) then
            print*,'Inconsistent mesh dimension ...'
        endif
        nx_tmp1 = nx
        allocate(choice_bot(nx_tmp1,2),choice_top(nx_tmp1,2))
        ! assume uniform distribution
        choice_bot = rbot
        choice_top = rtop
        do n=1,nrefine
            ! first refine
            nx_tmp2 = 2*nx_tmp1-1
            allocate(choice_bot_tmp(nx_tmp2,2),choice_top_tmp(nx_tmp2,2))
            call refine_boundary(choice_bot,nx_tmp1,choice_bot_tmp)
            call refine_boundary(choice_top,nx_tmp1,choice_top_tmp)
            ! second refine
            nx_tmp1 = 2*nx_tmp2-1
            deallocate(choice_bot,choice_top)
            allocate(choice_bot(nx_tmp1,2),choice_top(nx_tmp1,2))
            call refine_boundary(choice_bot_tmp,nx_tmp2,choice_bot)
            call refine_boundary(choice_top_tmp,nx_tmp2,choice_top)
            deallocate(choice_bot_tmp,choice_top_tmp)
        enddo
        num_choice = size(choice_bot,dim=1)
        print*,'Final shape of choice_bot & choice_top: ',size(choice_bot,dim=1),size(choice_bot,dim=2)
        print*,'Number of choice: ',num_choice
        ! uniform left and right: xprof,zprof_top,zprof_bot
        allocate(choice_left(ny,2),choice_right(ny,2))
        ! left/right boundary
        choice_left = rleft
        choice_right = rright
        ! initialize idx_bot,idx_top,idx_left,idx_right
        allocate(idx_bot(nx),idx_top(nx),idx_left(ny),idx_right(ny))
        call init_index(choice_bot,num_choice,idx_bot,x2d(1:nx,1),y2d(1:nx,1),nx,1)
        call init_index(choice_top,num_choice,idx_top,x2d(1:nx,ny),y2d(1:nx,ny),nx,1)
        call init_index(choice_left,ny,idx_left,x2d(1,1:ny),y2d(1,1:ny),ny,1)
        call init_index(choice_right,ny,idx_right,x2d(nx,1:ny),y2d(nx,1:ny),ny,1)
        print*,'TFI Initialization ...'
        call init_TFI(x2d,y2d,x2d,y2d,nx,ny)
        ! orthogonal mesh generation
        print*,'Compute grid ...'
        call compute_grid_new(x2d,y2d,idx_left,idx_right,idx_bot,idx_top,nx,ny,&
                              choice_left,ny,choice_right,ny,choice_bot,num_choice,choice_top,num_choice,&
                              boundary_fixed,boundary_symmetric,err_inner,err_outer,maxiter_inner,maxiter_outer,max_change,nsave)
        deallocate(choice_left,choice_right,choice_top,choice_bot)

    end subroutine GenOrthogonalMesh


    !----------------------------------------------
    ! Old version of orthogonal mesh generator:
    ! Use the inverse-Laplace solver without distortion function,
    ! orthogonal mesh is only applicable to uniform mesh.
    ! Left & right boundary grid distributions are applied through
    ! a combination of TFI mesh and orthogonal uniform mesh.
    !----------------------------------------------
    subroutine GenOrthogonalMesh_old(rleft,rright,rbot,rtop,x2d,y2d)
        real(8),dimension(:,:),intent(in) :: rleft,rright,rbot,rtop
        real(8),dimension(:,:),intent(inout) :: x2d,y2d
        integer :: n,i,j,nx,ny,nloop
        logical :: iloop = .True.
        ! orthogonal mesh generation
        integer num_choice,nx_tmp1,nx_tmp2
        integer,dimension(:),allocatable :: idx_left,idx_right,idx_bot,idx_top
        real(8),dimension(:,:),allocatable :: choice_left,choice_right,choice_bot,choice_top,choice_bot_tmp,choice_top_tmp
        real(8),dimension(:,:),allocatable :: x2d_new,z2d_new

        ! dimension
        nx = size(rbot,dim=1)
        ny = size(rleft,dim=1)
        if(size(x2d,dim=1).ne.nx .or. size(x2d,dim=2).ne.ny) then
            print*,'Inconsistent mesh dimension ...'
        endif

        nx_tmp1 = nx
        allocate(choice_bot(nx_tmp1,2),choice_top(nx_tmp1,2))
        ! assume uniform distribution
        choice_bot = rbot
        choice_top = rtop
        do n=1,nrefine
            ! first refine
            nx_tmp2 = 2*nx_tmp1-1
            allocate(choice_bot_tmp(nx_tmp2,2),choice_top_tmp(nx_tmp2,2))
            call refine_boundary(choice_bot,nx_tmp1,choice_bot_tmp)
            call refine_boundary(choice_top,nx_tmp1,choice_top_tmp)
            ! second refine
            nx_tmp1 = 2*nx_tmp2-1
            deallocate(choice_bot,choice_top)
            allocate(choice_bot(nx_tmp1,2),choice_top(nx_tmp1,2))
            call refine_boundary(choice_bot_tmp,nx_tmp2,choice_bot)
            call refine_boundary(choice_top_tmp,nx_tmp2,choice_top)
            deallocate(choice_bot_tmp,choice_top_tmp)
        enddo
        num_choice = size(choice_bot,dim=1)
        print*,'Final shape of choice_bot & choice_top: ',size(choice_bot,dim=1),size(choice_bot,dim=2)
        print*,'Number of choice: ',num_choice
        ! uniform left and right: xprof,zprof_top,zprof_bot
        allocate(choice_left(ny,2),choice_right(ny,2))
        ! left/right boundary
        choice_left = rleft
        choice_right = rright
        call apply_uniform(choice_left(:,1),rleft(1,1),rleft(ny,1),ny)  ! x
        call apply_uniform(choice_left(:,2),rleft(1,2),rleft(ny,2),ny)  ! y
        call apply_uniform(choice_right(:,1),rright(1,1),rright(ny,1),ny)   ! x
        call apply_uniform(choice_right(:,2),rright(1,2),rright(ny,2),ny)   ! y
        ! initialize idx_bot,idx_top,idx_left,idx_right
        allocate(idx_bot(nx),idx_top(nx),idx_left(ny),idx_right(ny))
        call init_index(choice_bot,num_choice,idx_bot,x2d(1:nx,1),y2d(1:nx,1),nx,1)
        call init_index(choice_top,num_choice,idx_top,x2d(1:nx,ny),y2d(1:nx,ny),nx,1)
        call init_index(choice_left,ny,idx_left,x2d(1,1:ny),y2d(1,1:ny),ny,1)
        call init_index(choice_right,ny,idx_right,x2d(nx,1:ny),y2d(nx,1:ny),ny,1)
        print*,'TFI Initialization ...'
        call init_TFI(x2d,y2d,x2d,y2d,nx,ny)
        ! orthogonal mesh generation
        print*,'Compute grid ...'
        call compute_grid(x2d,y2d,x2d,y2d,idx_left,idx_right,idx_bot,idx_top,idx_left,idx_right,idx_bot,idx_top,nx,ny,&
                          choice_left,ny,choice_right,ny,choice_bot,num_choice,choice_top,num_choice,&
                          err_inner,max_change,maxiter_outer,maxiter_inner,nsave,boundary_fixed,boundary_symmetric)
!        nloop = 0
!        do while(iloop)
!            nloop = nloop + 1
!            if (iconverge.eq.1) then
!                print*,'Converged !'
!                iloop = .False.
!            endif
!            if(nloop.ge.nmax) then
!                print*,'Maximum iteration #: ',nmax*nsave,' !!!'
!                iloop = .False.
!            endif
!        enddo
        deallocate(choice_left,choice_right,choice_top,choice_bot)
        ! apply the left/right grid distribution
        print*,'Apply left/right grid distribution ...'
        call init_index(rleft,ny,idx_left,x2d(1,1:ny),y2d(1,1:ny),ny,1)
        call init_index(rright,ny,idx_right,x2d(nx,1:ny),y2d(nx,1:ny),ny,1)
        call gridgen_checker(x2d,y2d,x2d,y2d,nx,ny)
    end subroutine GenOrthogonalMesh_old

    !----------------------------------------------
    ! Extract several orthogonal i lines from orthogonal mesh,
    ! use TFI to generate 2D mesh within each block.
    ! NOTE: apply for symmetry geometry only (rbot_in = -rtop_in)
    subroutine MultiBlockTFI(x2d,y2d,nx,ny,rbot_in,rtop_in)
        integer,intent(in) :: nx,ny
        real(8),dimension(nx,2),intent(inout) :: rbot_in,rtop_in
        real(8),dimension(nx,ny),intent(inout) :: x2d,y2d
        real(8),dimension(ny,2) :: rleft, rright
        real(8),dimension(:,:),allocatable :: xtmp,ytmp
        integer :: i,j,n,idx,nx0,nx1
        real(8) :: p0(2),p1(2),p2(2),ratio

        print*,'compute multiblock TFI mesh ...'
        allocate(xtmp(nx,ny),ytmp(nx,ny))
        rleft(:,1) = x2d(1,:)
        rleft(:,2) = y2d(1,:)
!        xtmp(1,:) = rleft(:,1)
!        ytmp(1,:) = rleft(:,2)
        nx0 = 1
        do n = 1,nblocks-1
!            print*,'Zone #: ',n
            p1(1) = x2d(1,1)+dble(n)*(x2d(nx,1)-x2d(1,1))/dble(nblocks)
            call BisectionSearch(rbot_in(:,1),nx,p1(1),idx)
            nx1 = idx
            p1 = (/rbot_in(idx,1),rbot_in(idx,2)/)
            rright(1,:) = p1
            rright(ny,:) = (/rtop_in(idx,1),rtop_in(idx,2)/)
            call BisectionSearch(x2d(:,1),nx,p1(1),idx)
            p0 = (/x2d(idx,1),y2d(idx,1)/)
            p2 = (/x2d(idx+1,1),y2d(idx+1,1)/)
            ratio = (p1(1)-p0(1))/(p2(1)-p0(1))
!            call LinearInterp2P(p0(1),p2(1),p0(2),p2(2),p1(1),p1(2))
            do j = 2,ny-1
                ! get p1
                p0 = (/x2d(idx,j),y2d(idx,j)/)
                p2 = (/x2d(idx+1,j),y2d(idx+1,j)/)
                p1(1) = p0(1)+ratio*(p2(1)-p0(1))
                call LinearInterp2P(p0(1),p2(1),p0(2),p2(2),p1(1),p1(2))
                rright(j,:) = p1
            enddo
!            rright(ny,:) = (/rtop_in(nx0+nx1-1,1),rtop_in(nx0+nx1-1,2)/)
            call tfi2d_struc(rleft,rright,rbot_in(nx0:nx1,1:ny),rtop_in(nx0:nx1,1:ny),xtmp(nx0:nx1,1:ny),ytmp(nx0:nx1,1:ny))
            nx0 = nx1
            rleft = rright
        enddo
        ! the last block
        rright(:,1) = x2d(nx,:)
        rright(:,2) = y2d(nx,:)
        nx1 = nx
        call tfi2d_struc(rleft,rright,rbot_in(nx0:nx1,1:ny),rtop_in(nx0:nx1,1:ny),xtmp(nx0:nx1,1:ny),ytmp(nx0:nx1,1:ny))
        open(unit=20,file="Check_BlockTFI.dat",status='unknown')
        write(20,*) 'Variables = x,z'
        write(20,*) 'Zone I = ',nx,', J = ',ny,', F=POINT'
        do j=1,ny
        do i=1,nx
            write(20,*) xtmp(i,j),ytmp(i,j)
        enddo
        enddo
        close(20)
        x2d = xtmp
        y2d = ytmp
        deallocate(xtmp,ytmp)
    end subroutine MultiBlockTFI

    subroutine CheckOrthogonalMesh(x2d,y2d,nx,ny)
        integer,intent(in) :: nx,ny
        real(8),dimension(nx,ny),intent(in) :: x2d,y2d
        integer :: i,j
        real(8),dimension(2,2) :: mm
        real(8),dimension(2) :: p11,p21,p01,p12,p10
        real(8) :: g12,hksi,heta
        character(200) :: tecname

        !          1 2 3     4            5        6        7    8    9    10
        tecname = 'x z angle aspect_ratio smooth_i smooth_j dxdi dxdj dydi dydj'
        call InitTec(1,nx,ny,1,10,tecname,0)
        vartmp(:,:,1,1) = x2d
        vartmp(:,:,1,2) = y2d
!        print*,'Start writing Check_OrthogonalMesh2d.dat ...'
        do j=1,ny
        do i=1,nx
            mm = Calmm_2D(i,j,x2d,y2d)
            g12 = mm(1,1)*mm(1,2) + mm(2,1)*mm(2,2)
            hksi = sqrt((mm(1,1))**2+(mm(2,1))**2)
            heta = sqrt((mm(1,2))**2+(mm(2,2))**2)
            vartmp(i,j,1,3) = acos(g12/(hksi*heta+1.d-20))*180./pi
            vartmp(i,j,1,4) = max(hksi,heta)/(min(hksi,heta)+1.d-20)
            vartmp(i,j,1,7) = mm(1,1)
            vartmp(i,j,1,8) = mm(1,2)
            vartmp(i,j,1,9) = mm(2,1)
            vartmp(i,j,1,10) = mm(2,2)
        enddo
        enddo

        vartmp(:,:,:,5:6) = 1.d0
        do j=2,ny-1
        do i=2,nx-1
            p11 = (/x2d(i,j),y2d(i,j)/)
            p21 = (/x2d(i+1,j),y2d(i+1,j)/)
            p01 = (/x2d(i-1,j),y2d(i-1,j)/)
            p12 = (/x2d(i,j+1),y2d(i,j+1)/)
            p10 = (/x2d(i,j-1),y2d(i,j-1)/)
            vartmp(i,j,1,5) = get_smoothness(p01,p11,p21)
            vartmp(i,j,1,6) = get_smoothness(p10,p11,p12)
        enddo
        enddo

        call WriteTec('Check_OrthogonalMesh2d.plt')
!        print*,'Finish writting Check_OrthogonalMesh2d.dat'
    end subroutine

    !----------------------------------------------
    ! Symmetric 2D mesh with respect to top boundary
    !----------------------------------------------
    subroutine mesh_symmetric_top(x2d,y2d,nx,ny,x2d_sym,y2d_sym,itp)
        integer,intent(in) :: nx,ny
        real(8),dimension(nx,ny),intent(in) :: x2d,y2d
        real(8),dimension(nx,ny),intent(out) :: x2d_sym,y2d_sym
        integer,intent(in),optional :: itp
        integer :: itype,i,j,ny_sym,ny_cut
        real(8) :: y_sym
        real(8),dimension(:),allocatable :: dist_left,dist_right
        real(8),dimension(:,:),allocatable :: x2d_tmp,y2d_tmp,x2d_tmp2,y2d_tmp2,rleft,rright

!        if(any(y2d(:,ny).ne.y2d(1,ny))) then
        if(any(abs(y2d(:,ny)-y2d(1,ny)).gt.real8_tol)) then
            print*,'mesh_symmetric_top: Top boundary NOT a symmetric plane!'
            stop
        endif
        if(.not.present(itp)) then
            itype=0
        else
            itype = itp
        endif
        y_sym = y2d(1,ny)
!        allocate(x2d_tmp(nx,ny),y2d_tmp(nx,ny))
        x2d_sym = x2d
        do j=1,ny
        do i=1,nx
            y2d_sym(i,j) = 2.*y_sym-y2d(i,ny-j+1)
        enddo
        enddo
        if(itype.ge.1) then
            y2d_sym = y2d_sym-y_sym
            if(itype.eq.2) then ! offset from symmetric plane
                ! step 1: create 2D mesh by symmetry with respect to j=1
                allocate(x2d_tmp(nx,2*ny-1),y2d_tmp(nx,2*ny-1))
                x2d_tmp(:,ny:2*ny-1) = x2d_sym; y2d_tmp(:,ny:2*ny-1) = y2d_sym
                call mesh_symmetric_bot(x2d_sym,y2d_sym,nx,ny,x2d_tmp(:,1:ny),y2d_tmp(:,1:ny))
                ! step 2: create 1D left/right distributions
                allocate(dist_left(2*ny-1),dist_right(2*ny-1))
                call copy_symmetric(y2d_sym(1,:),ny,dist_left)
                call copy_symmetric(y2d_sym(nx,:),ny,dist_right)
                ny_sym = 2*ny
                ! left & right
                allocate(rleft(ny_sym,2),rright(ny_sym,2))
                rleft(:,1) = x2d_sym(1,1); rright(:,1) = x2d_sym(nx,1)
                call apply_distribution(rleft(:,2),-y2d_sym(1,ny),y2d_sym(1,ny),ny_sym,dist_left)
                call apply_distribution(rright(:,2),-y2d_sym(nx,ny),y2d_sym(nx,ny),ny_sym,dist_right)
                deallocate(dist_left,dist_right)
                ! top & bot
!                allocate(rbot(nx,2),rtop(nx,2))
!                rtop(:,1) = x2d_sym(:,ny); rtop(:,2) = y2d_sym(:,ny)
!                rbot(:,1) = rtop(:,1); rbot(:,2) = -rtop(:,2)
!                call tfi2d_struc(rleft,rright,rbot,rtop,x2d_tmp,y2d_tmp)
                ! step 3: integrate both original and TFI mesh
                allocate(x2d_tmp2(nx,ny_sym),y2d_tmp2(nx,ny_sym))
                call mesh_left_right(x2d_tmp,y2d_tmp,rleft,rright,x2d_tmp2,y2d_tmp2)
                deallocate(x2d_tmp,y2d_tmp)
                ! step 4: cut the top half
                ny_cut = ny+1
                x2d_sym = x2d_tmp2(1:nx,ny_cut:ny_sym)
                y2d_sym = y2d_tmp2(1:nx,ny_cut:ny_sym)
                deallocate(x2d_tmp2,y2d_tmp2)
            endif
        endif
        ! debug
!        print*,'Start writing Check_MeshSymTop.dat ...'
        open(unit=20,file="Check_MeshSymTop.dat",status='unknown')
        write(20,*) 'Variables = x,z'
        write(20,*) 'Zone I = ',nx,', J = ',ny,', F=POINT'
        do j=1,ny
        do i=1,nx
            write(20,*) x2d_sym(i,j),y2d_sym(i,j)
        enddo
        enddo
        close(20)
    end subroutine mesh_symmetric_top

    subroutine mesh_left_right(x2d,y2d,rleft,rright,x2d_out,y2d_out)
        real(8),dimension(:,:),intent(in) :: x2d,y2d
        real(8),dimension(:,:),intent(in) :: rleft,rright
        real(8),dimension(:,:),intent(out) :: x2d_out,y2d_out
        integer :: nx,ny,ny_new,inter,i,j,m,n
        real(8),dimension(2) :: p00,p01,p10,p11,point
        real(8),dimension(:,:),allocatable :: rbot,rtop
        real(8) :: time_begin

        nx = size(x2d,dim=1); ny = size(x2d,dim=2)
        ny_new = size(rleft,dim=1)
        if(size(x2d_out,dim=1).ne.nx .or. size(x2d_out,dim=2).ne.ny_new) then
            print*,'mesh_left_right: Inconsistent in-out mesh dimension.'
            stop
        endif
        ! initialize with TFI
        allocate(rbot(nx,2),rtop(nx,2))
        rbot(:,1) = x2d(:,1); rbot(:,2) = y2d(:,1)
        rtop(:,1) = x2d(:,ny); rtop(:,2) = y2d(:,ny)
        call tfi2d_struc(rleft,rright,rbot,rtop,x2d_out,y2d_out)

        time_begin = omp_get_wtime()
        !$OMP PARALLEL &
        !$OMP private(m,n,p00,p01,p10,p11,inter, point)
        !$OMP DO
        do j = 2,ny_new-1
        do i = 2,nx-1
            do n = 1,ny-1
            do m = 1,nx-1
                p00 = (/x2d(i,n),y2d(i,n)/)     ! original mesh
                p01 = (/x2d(i,n+1),y2d(i,n+1)/)     ! original mesh
                p10 = (/x2d_out(m,j),y2d_out(m,j)/)     ! TFI mesh
                p11 = (/x2d_out(m+1,j),y2d_out(m+1,j)/)     ! TFI mesh
                call get_intersect_point(p00, p01, p10, p11, inter, point)
                if( inter==1 ) then
                    x2d_out(i,j) = point(1)
                    y2d_out(i,j) = point(2)
                    exit
                endif
            enddo
            if( inter==1 ) exit
            enddo
        enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        print*,' '
        print*,'mesh_left_right OMP time [sec]: ',omp_get_wtime()-time_begin
        print*,' '
        deallocate(rbot,rtop)
    end subroutine mesh_left_right

    !----------------------------------------------
    ! Symmetric 2D mesh with respect to bottom boundary
    !----------------------------------------------
    subroutine mesh_symmetric_bot(x2d,y2d,nx,ny,x2d_sym,y2d_sym,itp)
        integer,intent(in) :: nx,ny
        real(8),dimension(nx,ny),intent(in) :: x2d,y2d
        real(8),dimension(nx,ny),intent(out) :: x2d_sym,y2d_sym
        integer,intent(in),optional :: itp
        integer :: itype,i,j
        real(8) :: y_sym

        if(any(abs(y2d(:,1)-y2d(1,1)).gt.real8_tol)) then
            print*,'mesh_symmetric_bot: Bottom boundary NOT a symmetric plane!'
            stop
        endif
        if(.not.present(itp)) then
            itype=0
        else
            itype = itp
        endif
        y_sym = y2d(1,1)
        x2d_sym = x2d
        do j=1,ny
        do i=1,nx
            y2d_sym(i,j) = 2.*y_sym-y2d(i,ny-j+1)
        enddo
        enddo
        ! debug
!        print*,'Start writing Check_MeshSymBot.dat ...'
        open(unit=20,file="Check_MeshSymBot.dat",status='unknown')
        write(20,*) 'Variables = x,z'
        write(20,*) 'Zone I = ',nx,', J = ',ny,', F=POINT'
        do j=1,ny
        do i=1,nx
            write(20,*) x2d_sym(i,j),y2d_sym(i,j)
        enddo
        enddo
        close(20)
    end subroutine mesh_symmetric_bot

    !----------------------------------------------
    ! Structured 2D mesh using TFI x2d(idim,jdim), &
    !                              y2d(idim,jdim)
    !----------------------------------------------
    subroutine tfi2d_struc(rleft,rright,rbot,rtop,x2d,y2d)
        real(8),dimension(:,:),intent(in) :: rleft,rright,rbot,rtop
        real(8),dimension(:,:),intent(out) :: x2d,y2d
        real(8) :: dx,dy,ksi,eta
        integer :: i,j,idim,jdim
        if(size(rleft,dim=1).ne.size(rright,dim=1) .or. size(rbot,dim=1).ne.size(rtop,dim=1)) then
            print*,'Dimension of opposite boundary NOT match. '
            stop
        endif
        idim = size(rbot,dim=1); jdim = size(rleft,dim=1)
        if(rbot(1,1).ne.rleft(1,1) .or. rbot(idim,1).ne.rright(1,1) .or. rright(jdim,1).ne.rtop(idim,1) .or. rleft(jdim,1).ne.rtop(1,1)) then
            print*,'x coordinate at four vertices NOT consistency. '
            print*,'bottom left: ',rbot(1,1),rleft(1,1)
            print*,'top left: ',rtop(1,1),rleft(jdim,1)
            print*,'top right: ',rtop(idim,1),rright(jdim,1)
            print*,'bottom right: ',rbot(idim,1),rright(1,1)
            stop
        endif
        if(abs(rbot(1,2)-rleft(1,2)).gt.real8_tol .or. abs(rbot(idim,2)-rright(1,2)).gt.real8_tol .or. &
           abs(rright(jdim,2)-rtop(idim,2)).gt.real8_tol .or. abs(rleft(jdim,2)-rtop(1,2)).gt.real8_tol) then
            print*,'y coordinate at four vertices NOT consistency. '
            print*,'bottom left: ',rbot(1,2),rleft(1,2)
            print*,'top left: ',rtop(1,2),rleft(jdim,2)
            print*,'top right: ',rtop(idim,2),rright(jdim,2)
            print*,'bottom right: ',rbot(idim,2),rright(1,2)
            stop
        endif
        x2d(1,:) = rleft(:,1); x2d(idim,:) = rright(:,1)
        x2d(:,1) = rbot(:,1); x2d(:,jdim) = rtop(:,1)
        y2d(1,:) = rleft(:,2); y2d(idim,:) = rright(:,2)
        y2d(:,1) = rbot(:,2); y2d(:,jdim) = rtop(:,2)
        dx = 1./dble(idim-1)
        dy = 1./dble(jdim-1)
        do j = 2,jdim-1
            eta = dble(j-1)*dy
            do i = 2,idim-1
                ksi = dble(i-1)*dx
                x2d(i,j) = (1.0-ksi)*rleft(j,1)+ksi*rright(j,1)+(1.0-eta)*rbot(i,1)+eta*rtop(i,1) &
                          -(1.0-ksi)*(1.0-eta)*rbot(1,1)-(1.0-ksi)*eta*rtop(1,1)-(1.0-eta)*ksi*rbot(idim,1)-ksi*eta*rtop(idim,1)
                y2d(i,j) = (1.0-ksi)*rleft(j,2)+ksi*rright(j,2)+(1.0-eta)*rbot(i,2)+eta*rtop(i,2) &
                          -(1.0-ksi)*(1.0-eta)*rbot(1,2)-(1.0-ksi)*eta*rtop(1,2)-(1.0-eta)*ksi*rbot(idim,2)-ksi*eta*rtop(idim,2)
            enddo
        enddo
    end subroutine tfi2d_struc

    !----------------------------------------------
    ! Structured 2D mesh using TFI x2d(jdim,idim), &
    !                              y2d(jdim,idim)
    !----------------------------------------------
    subroutine tfi2d_struc_shift(rleft,rright,rbot,rtop,x2d,y2d)
        real(8),dimension(:,:),intent(in) :: rleft,rright,rbot,rtop
        real(8),dimension(:,:),intent(out) :: x2d,y2d
        real(8) :: dx,dy,ksi,eta
        integer :: i,j,idim,jdim
        if(size(rleft,dim=1).ne.size(rright,dim=1) .or. size(rbot,dim=1).ne.size(rtop,dim=1)) then
            print*,'Dimension of opposite boundary NOT match. '
            stop
        endif
        idim = size(rbot,dim=1); jdim = size(rleft,dim=1)
        if(rbot(1,1).ne.rleft(1,1) .or. rbot(idim,1).ne.rright(1,1) .or. rright(jdim,1).ne.rtop(idim,1) .or. rleft(jdim,1).ne.rtop(1,1)) then
            print*,'x coordinate at four vertices NOT consistency. '
            stop
        endif
        if(rbot(1,2).ne.rleft(1,2) .or. rbot(idim,2).ne.rright(1,2) .or. rright(jdim,2).ne.rtop(idim,2) .or. rleft(jdim,2).ne.rtop(1,2)) then
            print*,'y coordinate at four vertices NOT consistency. '
            stop
        endif
        x2d(:,1) = rleft(:,1); x2d(:,idim) = rright(:,1)
        x2d(1,:) = rbot(:,1); x2d(jdim,:) = rtop(:,1)
        y2d(:,1) = rleft(:,2); y2d(:,idim) = rright(:,2)
        y2d(1,:) = rbot(:,2); y2d(jdim,:) = rtop(:,2)
        dx = 1./dble(idim-1)
        dy = 1./dble(jdim-1)
        do j = 2,jdim-1
            eta = dble(j-1)*dy
            do i = 2,idim-1
                ksi = dble(i-1)*dx
                x2d(j,i) = (1.0-ksi)*rleft(j,1)+ksi*rright(j,1)+(1.0-eta)*rbot(i,1)+eta*rtop(i,1) &
                          -(1.0-ksi)*(1.0-eta)*rbot(1,1)-(1.0-ksi)*eta*rtop(1,1)-(1.0-eta)*ksi*rbot(idim,1)-ksi*eta*rtop(idim,1)
                y2d(j,i) = (1.0-ksi)*rleft(j,2)+ksi*rright(j,2)+(1.0-eta)*rbot(i,2)+eta*rtop(i,2) &
                          -(1.0-ksi)*(1.0-eta)*rbot(1,2)-(1.0-ksi)*eta*rtop(1,2)-(1.0-eta)*ksi*rbot(idim,2)-ksi*eta*rtop(idim,2)
            enddo
        enddo
    end subroutine tfi2d_struc_shift

    !----------------------------------------------
    ! Create a 1D array with evenly spaced elements
    !----------------------------------------------
    subroutine apply_uniform(x,x_start,x_end,ilen)
        real(8),dimension(:),intent(out) :: x
        real(8),intent(in) :: x_start,x_end
        integer,intent(in) :: ilen
        real(8) :: dx
        integer :: i
        dx = (x_end-x_start)/dble(ilen-1)
        x(1:ilen) = [(x_start+dble(i-1)*dx,i=1,ilen)]
    end subroutine apply_uniform

    !----------------------------------------------
    ! Create a 1D array with a given distribution
    !----------------------------------------------
    subroutine apply_distribution(x,x_start,x_end,ilen,dist)
        real(8),dimension(:),intent(in) :: dist
        real(8),dimension(:),intent(out) :: x
        real(8),intent(in) :: x_start,x_end
        integer,intent(in) :: ilen
        real(8) :: dx
        integer :: i,ndist
        real(8),dimension(:),allocatable :: xindx_old,xindx_new,dist_new,dist_unit

        ndist = size(dist,dim=1)
        ! normalized 1D distribution
        allocate(dist_unit(ndist))
        call normalized1d(dist,dist_unit)
        ! interpolate
        allocate(dist_new(ilen))
        if(ndist.ne.ilen) then
            allocate(xindx_old(ndist),xindx_new(ilen))
            forall(i=1:ndist) xindx_old(i) = dble(i-1)/dble(ndist-1)
            forall(i=1:ilen) xindx_new(i) = dble(i-1)/dble(ilen-1)
            call SplineInterp(xindx_old,dist_unit,ndist, xindx_new,dist_new,ilen, 0)
            deallocate(xindx_old,xindx_new)
        else
            dist_new = dist
        endif
        ! apply distribution
        forall(i=1:ilen) x(i) = x_start+dist_new(i)*(x_end-x_start)
        x(ilen) = x_end
        deallocate(dist_new)
    end subroutine apply_distribution

    !----------------------------------------------
    ! Create a 1D array with symmetric part
    !----------------------------------------------
    subroutine copy_symmetric(x_in,nx,xunit_out)
        integer,intent(in) :: nx
        real(8),dimension(nx),intent(in) :: x_in
        real(8),dimension(2*nx-1),intent(out) :: xunit_out
        integer :: i
        real(8),dimension(:),allocatable :: xunit_in,x_out
        ! normalize 1D distribution
        allocate(xunit_in(nx))
        call normalized1d(x_in,xunit_in)
        ! copy and flipped
        allocate(x_out(2*nx-1))
        do i=1,nx
            x_out(i+nx-1) = xunit_in(i)
            x_out(i) = -xunit_in(nx-i+1)
        enddo
        ! normalized and return
        call normalized1d(x_out,xunit_out)
        deallocate(xunit_in,x_out)
    end subroutine copy_symmetric

    !----------------------------------------------
    ! normalized the 1D array to [0:1]
    !----------------------------------------------
    subroutine normalized1d(x_in,x_out)
        real(8),dimension(:),intent(in) :: x_in
        real(8),dimension(:),intent(out) :: x_out
        integer :: nx,i
        real(8) :: xscale
        nx = size(x_in)
        if(size(x_in).ne.size(x_out)) then
            print*,'Inconsistent array dimension. STOP!'
            stop
        endif
        xscale = x_in(nx)-x_in(1)
        forall(i=1:nx) x_out(i) = (x_in(i)-x_in(1))/xscale
    end subroutine normalized1d
end module MGRIDGEN2D
