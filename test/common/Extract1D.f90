!===============================================================================
! Mod_Extract1D — Estrazione grandezze 1D da output Tecplot ARES
!
! Variabili primitive: (p, h, u, v, w).  Lo stato termodinamico completo
! e' calcolato direttamente dalle tabelle 2D (p,h) del gas reale (FLINT).
!
! Uso:  extract1d  field.tec  wall.tec  1d.dat  INPUT_folder
!===============================================================================
module Mod_Extract1D
  implicit none
  private

  integer, parameter, public :: dp = selected_real_kind(15,307)
  integer, parameter :: MAXVAR = 80, LLEN = 1024

  ! ---------- Tipi Tecplot ----------
  type, public :: tec_zone_t
    character(len=64) :: name = ''
    integer :: ni=0, nj=0, nk=0, nx=0, ny=0, nz=0
    real(dp), allocatable :: xn(:,:,:), yn(:,:,:), zn(:,:,:)
    real(dp), allocatable :: v(:,:,:,:)
  end type

  type, public :: tec_file_t
    integer :: nzones=0, nvar=0
    character(len=64) :: varname(MAXVAR)
    type(tec_zone_t), allocatable :: zone(:)
  end type

  ! ---------- Tabelle EOS 2D gas reale ----------
  type, public :: eos_t
    ! griglia (p,h)
    real(dp) :: plo=0, dpt=0, hlo=0, dht=0
    integer  :: np=0, nh=0
    real(dp), allocatable :: rho(:,:), T(:,:), rh(:,:), cp(:,:)
    real(dp), allocatable :: s(:,:), snd(:,:)
    real(dp), allocatable :: mu(:,:), kth(:,:)
    ! inversa (p,T)->h costruita da eos%T (mimica FLINT/ph2pT):
    ! griglia T per colonna di p — Tmin2(i), deltaT(i), i=0..np
    real(dp), allocatable :: Tmin2(:), deltaT(:)
    real(dp), allocatable :: h_pT(:,:)        ! (0:np, 0:nh)
    logical :: have_pT = .false.
    logical :: ok = .false.
  end type

  public :: extract_1d, read_eos_tables

contains

  ! ===========================  TECPLOT READER  ============================

  subroutine read_next_float(iu, val, ierr)
    integer, intent(in) :: iu
    real(dp), intent(out) :: val
    integer, intent(out) :: ierr
    character(len=LLEN) :: line
    character :: c
    integer :: ios, p
    do
      read(iu,'(A)',iostat=ios) line
      if (ios /= 0) then; ierr=1; return; endif
      ! Find first non-blank character; skip lines that cannot be numeric data
      p=1; do while (p<=len_trim(line).and.line(p:p)==' '); p=p+1; end do
      if (p>len_trim(line)) cycle          ! blank line
      c=line(p:p)
      if (.not.(c=='+'.or.c=='-'.or.c=='.'.or.(c>='0'.and.c<='9'))) cycle
      read(line,*,iostat=ios) val
      if (ios == 0) then; ierr=0; return; endif
    end do
  end subroutine

  subroutine parse_varnames(line, names, nn)
    character(len=*), intent(in) :: line
    character(len=64), intent(out) :: names(MAXVAR)
    integer, intent(out) :: nn
    integer :: i, i1, ie, L
    nn = 0; L = len_trim(line); i = 1
    do while (i <= L)
      if (line(i:i) == '"') then
        i1 = i+1; ie = index(line(i1:),'"')
        if (ie > 0) then
          ie = i1+ie-2; nn = nn+1; names(nn) = line(i1:ie); i = ie+2
        else; exit; endif
      else; i = i+1; endif
    end do
  end subroutine

  subroutine extract_int(line, key, val)
    character(len=*), intent(in) :: line, key
    integer, intent(inout) :: val
    integer :: p, ios
    character(len=20) :: buf
    p = index(line, trim(key))
    if (p > 0) then
      p = p + len_trim(key); buf = ''
      do while (p <= len_trim(line))
        if (line(p:p)>='0' .and. line(p:p)<='9') then
          buf = trim(buf)//line(p:p)
        else if (len_trim(buf)>0) then; exit; endif
        p = p+1
      end do
      if (len_trim(buf)>0) read(buf,*,iostat=ios) val
    endif
  end subroutine

  subroutine to_upper(s)
    character(len=*), intent(inout) :: s
    integer :: i, c
    do i=1,len_trim(s); c=ichar(s(i:i))
      if (c>=ichar('a').and.c<=ichar('z')) s(i:i)=char(c-32)
    end do
  end subroutine

  subroutine parse_zone_header(line, zname, ni, nj, nk)
    character(len=*), intent(in) :: line
    character(len=64), intent(out) :: zname
    integer, intent(out) :: ni, nj, nk
    integer :: p1, p2
    character(len=LLEN) :: ul
    zname=''; ni=1; nj=1; nk=1; ul=line; call to_upper(ul)
    p1=index(line,'T='); if(p1==0) p1=index(line,'T =')
    if (p1>0) then
      p1=index(line(p1:),'=')+p1; p2=index(line(p1:),',')
      if (p2>0) then; zname=adjustl(line(p1:p1+p2-2))
      else; zname=adjustl(line(p1:)); endif
      zname=trim(adjustl(zname))
    endif
    call extract_int(ul,'I=',ni); call extract_int(ul,'J=',nj); call extract_int(ul,'K=',nk)
  end subroutine

  subroutine read_block_flat(iu, arr, nn, ierr)
    integer, intent(in) :: iu, nn
    real(dp), intent(inout) :: arr(*)
    integer, intent(out) :: ierr
    integer :: i
    do i=1,nn; call read_next_float(iu,arr(i),ierr); if(ierr/=0) return; end do
  end subroutine

  subroutine read_var_block(iu, varr, nv, sv, n1, n2, n3, ierr)
    integer, intent(in) :: iu, nv, sv, n1, n2, n3
    real(dp), intent(inout) :: varr(nv,n1,n2,n3)
    integer, intent(out) :: ierr
    integer :: i, j, k
    do k=1,n3; do j=1,n2; do i=1,n1
      call read_next_float(iu,varr(sv,i,j,k),ierr); if(ierr/=0) return
    end do; end do; end do
  end subroutine

  subroutine read_tec_ascii(filename, tec, ierr)
    character(len=*), intent(in) :: filename
    type(tec_file_t), intent(out) :: tec
    integer, intent(out) :: ierr
    integer :: iu,ios,iz,sv,n,ndir,ni,nj,nk,nx,ny,nz,nnode
    character(len=LLEN) :: line, varline
    character(len=64) :: allnames(MAXVAR)
    integer :: nall
    logical :: fv

    ierr=0
    open(newunit=iu,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios/=0) then; write(*,'(2A)') ' ERRORE: ',trim(filename); ierr=1; return; endif

    tec%nzones=0; fv=.false.; varline=''
    do
      read(iu,'(A)',iostat=ios) line; if(ios/=0) exit
      if (.not.fv) then
        if (index(line,'VARIABLES')>0 .or. index(line,'variables')>0 &
            .or. index(line,'Variables')>0) then
          n=index(line,'='); if(n>0) varline=trim(varline)//' '//trim(line(n+1:))
          fv=.true.; cycle
        endif
      endif
      if (fv .and. tec%nzones==0) then
        if (index(line,'ZONE')==0 .and. index(line,'Zone')==0 &
            .and. index(line,'zone')==0 .and. index(line,'"')>0) then
          varline=trim(varline)//' '//trim(line); cycle
        endif
      endif
      if (index(line,'ZONE')>0 .or. index(line,'Zone')>0 .or. index(line,'zone')>0) &
        tec%nzones = tec%nzones+1
    end do
    if (tec%nzones==0) then; ierr=2; close(iu); return; endif

    call parse_varnames(varline, allnames, nall)
    ndir=0
    do n=1,min(nall,3)
      if (allnames(n)=='x'.or.allnames(n)=='y'.or.allnames(n)=='z' &
          .or.allnames(n)=='X'.or.allnames(n)=='Y'.or.allnames(n)=='Z') then
        ndir=ndir+1
      else; exit; endif
    end do
    if (ndir==0) ndir=3
    tec%nvar = nall-ndir
    do n=1,tec%nvar; tec%varname(n)=allnames(ndir+n); end do

    allocate(tec%zone(tec%nzones)); rewind(iu); iz=0
    do
      read(iu,'(A)',iostat=ios) line; if(ios/=0) exit
      if (index(line,'ZONE')>0 .or. index(line,'Zone')>0 .or. index(line,'zone')>0) then
        iz=iz+1
        call parse_zone_header(line,tec%zone(iz)%name, &
          tec%zone(iz)%ni,tec%zone(iz)%nj,tec%zone(iz)%nk)
        tec%zone(iz)%nx=max(tec%zone(iz)%ni-1,1)
        tec%zone(iz)%ny=max(tec%zone(iz)%nj-1,1)
        tec%zone(iz)%nz=max(tec%zone(iz)%nk-1,1)
      endif
    end do

    rewind(iu)
    do iz=1,tec%nzones
      ni=tec%zone(iz)%ni; nj=tec%zone(iz)%nj; nk=tec%zone(iz)%nk
      nx=tec%zone(iz)%nx; ny=tec%zone(iz)%ny; nz=tec%zone(iz)%nz
      nnode=ni*nj*nk
      allocate(tec%zone(iz)%xn(ni,nj,nk))
      allocate(tec%zone(iz)%yn(ni,nj,nk))
      allocate(tec%zone(iz)%zn(ni,nj,nk))
      if (tec%nvar>0) allocate(tec%zone(iz)%v(tec%nvar,nx,ny,nz))
      call read_block_flat(iu,tec%zone(iz)%xn,nnode,ierr); if(ierr/=0) return
      if (ndir>=2) then
        call read_block_flat(iu,tec%zone(iz)%yn,nnode,ierr); if(ierr/=0) return
      else; tec%zone(iz)%yn=0; endif
      if (ndir>=3) then
        call read_block_flat(iu,tec%zone(iz)%zn,nnode,ierr); if(ierr/=0) return
      else; tec%zone(iz)%zn=0; endif
      do sv=1,tec%nvar
        call read_var_block(iu,tec%zone(iz)%v,tec%nvar,sv,nx,ny,nz,ierr)
        if(ierr/=0) return
      end do
    end do
    close(iu)
  end subroutine

  function find_var(tec, name) result(idx)
    type(tec_file_t), intent(in) :: tec
    character(len=*), intent(in) :: name
    integer :: idx, i
    character(len=64) :: a, b
    idx=0; b=name; call to_upper(b)
    do i=1,tec%nvar
      a=tec%varname(i); call to_upper(a)
      if(trim(a)==trim(b)) then; idx=i; return; endif
    end do
  end function

  ! ============================  GEOMETRIA  ================================
  ! I-face: diagonali incrociate con fattore 0.5 — coerente con Lib_Metrics.f90
  ! Hex volume: 5 tetraedri — stessa decomposizione di Lib_Metrics.f90

  pure subroutine quad_area_normal(p1,p2,p3,p4, area, normal)
    real(dp), intent(in) :: p1(3),p2(3),p3(3),p4(3)
    real(dp), intent(out) :: area, normal(3)
    real(dp) :: d1(3),d2(3),cr(3)
    d1 = p3-p1; d2 = p4-p2
    cr(1) = 0.5_dp*(d1(2)*d2(3)-d1(3)*d2(2))
    cr(2) = 0.5_dp*(d1(3)*d2(1)-d1(1)*d2(3))
    cr(3) = 0.5_dp*(d1(1)*d2(2)-d1(2)*d2(1))
    area = sqrt(cr(1)**2+cr(2)**2+cr(3)**2)
    if (area>0) then; normal=cr/area; else; normal=0; endif
  end subroutine

  pure function tet_vol(a,b,c,d) result(vol)
    real(dp), intent(in) :: a(3),b(3),c(3),d(3)
    real(dp) :: vol, e1(3),e2(3),e3(3)
    e1=b-a; e2=c-a; e3=d-a
    vol = abs(e1(1)*(e2(2)*e3(3)-e2(3)*e3(2)) &
            - e1(2)*(e2(1)*e3(3)-e2(3)*e3(1)) &
            + e1(3)*(e2(1)*e3(2)-e2(2)*e3(1))) / 6
  end function

  pure function hex_vol(a,b,c,d,e,f,g,h) result(vol)
    real(dp), intent(in) :: a(3),b(3),c(3),d(3),e(3),f(3),g(3),h(3)
    real(dp) :: vol
    ! Decomposizione ARES (Lib_Metrics.f90): 5 tetraedri
    !   tvol(1,2,3,5) + tvol(2,4,3,8) + tvol(5,8,6,2) + tvol(5,7,8,3) + tvol(5,8,2,3)
    ! Mapping Tecplot -> ARES:
    !   a=n0=N1, b=n1=N5, c=n2=N7, d=n3=N3, e=n4=N2, f=n5=N6, g=n6=N8, h=n7=N4
    vol = tet_vol(a,b,d,e) + tet_vol(b,c,d,g) &
        + tet_vol(e,g,f,b) + tet_vol(e,h,g,d) + tet_vol(e,g,b,d)
  end function

  ! ==========================  PARETE  =====================================

  subroutine parse_wall_name(name, ib, iface, ierr)
    character(len=*), intent(in) :: name
    integer, intent(out) :: ib, iface, ierr
    integer :: pb, pf, ios
    character(len=64) :: uname
    ierr=1; ib=0; iface=0
    uname=name; call to_upper(uname)
    ! Format ARES: "BLOCK:<n> FACE   <m>"
    pb=index(uname,'BLOCK:'); pf=index(uname,'FACE')
    if (pb>0 .and. pf>pb) then
      read(uname(pb+6:pf-1),*,iostat=ios) ib;    if(ios/=0) return
      read(uname(pf+4:),    *,iostat=ios) iface;  if(ios/=0) return
      ierr=0
      return
    endif
    ! Fallback: "B<n>F<m>"
    pb=index(uname,'B'); pf=index(uname,'F')
    if (pb>0 .and. pf>pb) then
      read(uname(pb+1:pf-1),*,iostat=ios) ib;    if(ios/=0) return
      read(uname(pf+1:),    *,iostat=ios) iface;  if(ios/=0) return
      ierr=0
    endif
  end subroutine

  function find_wall_zone(wall, ib, iface) result(idx)
    type(tec_file_t), intent(in) :: wall
    integer, intent(in) :: ib, iface
    integer :: idx, iz, wb, wf, ie
    idx=0
    do iz=1,wall%nzones
      call parse_wall_name(wall%zone(iz)%name, wb, wf, ie)
      if (ie==0 .and. wb==ib .and. wf==iface) then; idx=iz; return; endif
    end do
  end function

  ! Restituisce somme NON normalizzate: Aw, sum(|tau|*dA), sum(qw*dA), sum(Tw*dA)
  subroutine wall_station_sums(wz, iface, ic, it1,it2,it3,it_sc,iTw,iqw,ihs, &
                               Aw, tau_s, qw_s, Tw_s, hs_s)
    type(tec_zone_t), intent(in) :: wz
    integer, intent(in) :: iface, ic, it1, it2, it3, it_sc, iTw, iqw, ihs
    real(dp), intent(out) :: Aw, tau_s, qw_s, Tw_s, hs_s
    integer :: jc, kc
    real(dp) :: p1(3),p2(3),p3(3),p4(3), Aq, nrm(3), tw
    logical :: have_vec
    Aw=0; tau_s=0; qw_s=0; Tw_s=0; hs_s=0
    have_vec = (it1>0 .and. it2>0 .and. it3>0)

    if (iface==3 .or. iface==4) then
      do kc=1,wz%nz
        p1=(/wz%xn(ic,  1,kc  ),wz%yn(ic,  1,kc  ),wz%zn(ic,  1,kc  )/)
        p2=(/wz%xn(ic+1,1,kc  ),wz%yn(ic+1,1,kc  ),wz%zn(ic+1,1,kc  )/)
        p3=(/wz%xn(ic+1,1,kc+1),wz%yn(ic+1,1,kc+1),wz%zn(ic+1,1,kc+1)/)
        p4=(/wz%xn(ic,  1,kc+1),wz%yn(ic,  1,kc+1),wz%zn(ic,  1,kc+1)/)
        call quad_area_normal(p1,p2,p3,p4, Aq, nrm)
        tw=0
        if (have_vec) then
          tw=sqrt(wz%v(it1,ic,1,kc)**2+wz%v(it2,ic,1,kc)**2+wz%v(it3,ic,1,kc)**2)
        elseif (it_sc>0) then
          tw=wz%v(it_sc,ic,1,kc)
        endif
        Aw   = Aw   + Aq
        tau_s= tau_s + tw*Aq
        if(iqw>0) qw_s = qw_s + wz%v(iqw,ic,1,kc)*Aq
        if(iTw>0) Tw_s = Tw_s + wz%v(iTw,ic,1,kc)*Aq
        if(ihs>0) hs_s = hs_s + wz%v(ihs,ic,1,kc)*Aq
      end do
    else if (iface==5 .or. iface==6) then
      do jc=1,wz%ny
        p1=(/wz%xn(ic,  jc,  1),wz%yn(ic,  jc,  1),wz%zn(ic,  jc,  1)/)
        p2=(/wz%xn(ic+1,jc,  1),wz%yn(ic+1,jc,  1),wz%zn(ic+1,jc,  1)/)
        p3=(/wz%xn(ic+1,jc+1,1),wz%yn(ic+1,jc+1,1),wz%zn(ic+1,jc+1,1)/)
        p4=(/wz%xn(ic,  jc+1,1),wz%yn(ic,  jc+1,1),wz%zn(ic,  jc+1,1)/)
        call quad_area_normal(p1,p2,p3,p4, Aq, nrm)
        tw=0
        if (have_vec) then
          tw=sqrt(wz%v(it1,ic,jc,1)**2+wz%v(it2,ic,jc,1)**2+wz%v(it3,ic,jc,1)**2)
        elseif (it_sc>0) then
          tw=wz%v(it_sc,ic,jc,1)
        endif
        Aw   = Aw   + Aq
        tau_s= tau_s + tw*Aq
        if(iqw>0) qw_s = qw_s + wz%v(iqw,ic,jc,1)*Aq
        if(iTw>0) Tw_s = Tw_s + wz%v(iTw,ic,jc,1)*Aq
        if(ihs>0) hs_s = hs_s + wz%v(ihs,ic,jc,1)*Aq
      end do
    endif
  end subroutine

  ! ===========================  TABELLE EOS  ===============================

  subroutine find_file(folder, base, path, found)
    character(len=*), intent(in) :: folder, base
    character(len=512), intent(out) :: path
    logical, intent(out) :: found
    path=trim(folder)//'/'//trim(base); inquire(file=trim(path),exist=found)
    if (.not.found) then
      path=trim(folder)//'\'//trim(base); inquire(file=trim(path),exist=found)
    endif
  end subroutine

  ! Legge la rugosita' equivalente sand-grain ks [m] da input.ini.
  ! Ogni blocco wall ([qw], ...) puo' contenere una riga 'ks = ...'; viene
  ! restituito il valore massimo trovato (la parete piu' rugosa, rilevante per
  ! attrito e scambio termico).  Se il file manca o non c'e' alcun 'ks', hs=0.
  ! Gli esponenti Fortran 'd'/'D' (es. 1.60d-4) sono gestiti nativamente dal
  ! read list-directed.
  subroutine read_ks_from_ini(path, ks)
    character(len=*), intent(in) :: path
    real(dp), intent(out) :: ks
    integer :: iu, ios, p
    real(dp) :: val
    character(len=LLEN) :: line
    character(len=64)   :: key
    logical :: ex
    ks = 0.0_dp
    inquire(file=trim(path), exist=ex)
    if (.not.ex) then
      write(*,'(2A)') '   AVVISO: input.ini non trovato per hs: ', trim(path)
      return
    endif
    open(newunit=iu, file=trim(path), status='old', action='read', iostat=ios)
    if (ios/=0) return
    do
      read(iu,'(A)',iostat=ios) line; if(ios/=0) exit
      p = index(line,'=')
      if (p <= 1) cycle
      key = adjustl(line(1:p-1)); call to_upper(key)
      if (trim(key) /= 'KS') cycle
      read(line(p+1:),*,iostat=ios) val
      if (ios==0 .and. val > ks) ks = val
    end do
    close(iu)
  end subroutine

  ! Skip Tecplot ASCII header (TITLE / VARIABLES / ZONE ...) and extract
  ! I=, J= from the ZONE line.  Leaves file positioned at the first data row.
  ! Layout contract (FLINT, see Load_ThermoTransport.f90): I-axis = p, J-axis = h.
  ! POINT format here is written with j (=h) varying fastest, i.e. for each p
  ! value Nh+1 consecutive rows are written stepping through h.
  subroutine read_tec_table_header(iu, ni_zone, nj_zone, ierr)
    integer, intent(in)  :: iu
    integer, intent(out) :: ni_zone, nj_zone, ierr
    character(len=LLEN)  :: line, ul
    integer              :: ios
    ni_zone = 0; nj_zone = 0; ierr = 1
    do
      read(iu,'(A)',iostat=ios) line
      if (ios /= 0) return
      ul = line; call to_upper(ul)
      if (index(ul,'ZONE') > 0) then
        call extract_int(ul,'I=',ni_zone)
        call extract_int(ul,'J=',nj_zone)
        if (ni_zone > 0 .and. nj_zone > 0) ierr = 0
        return
      endif
    end do
  end subroutine

  subroutine read_eos_tables(folder, eos, ierr)
    character(len=*), intent(in) :: folder
    type(eos_t), intent(out) :: eos
    integer, intent(out) :: ierr
    integer :: iu, ios, i, j
    character(len=512) :: fn
    integer :: ni_zone, nj_zone, t_ni, t_nj
    real(dp), allocatable :: tmpblk(:,:)
    logical :: ex

    ierr=0

    ! --- thermo.dat  Tecplot BLOCK, I=p (fastest), J=h (slowest)
    !     Vars: Pressure, Enthalpy, Density, Temperature, dRho/dT, dRho/dh,
    !           Cp, Entropy, dRho/dp, SoundSpeed
    !     Storage: eos%rho(i,j), i=p index (0..Np), j=h index (0..Nh).
    call find_file(folder,'thermo.dat',fn,ex)
    if (.not.ex) then; write(*,'(A)') ' ERRORE: thermo.dat non trovato'; ierr=2; return; endif
    open(newunit=iu,file=trim(fn),status='old',action='read',iostat=ios)
    if(ios/=0) then; ierr=2; return; endif
    call read_tec_table_header(iu, ni_zone, nj_zone, ios)
    if (ios/=0) then
      write(*,'(A)') ' ERRORE: header Tecplot non valido in thermo.dat'
      ierr=2; close(iu); return
    endif
    eos%np = ni_zone - 1
    eos%nh = nj_zone - 1
    allocate(eos%rho(0:eos%np,0:eos%nh), eos%T(0:eos%np,0:eos%nh))
    allocate(eos%rh(0:eos%np,0:eos%nh),  eos%cp(0:eos%np,0:eos%nh))
    allocate(eos%s(0:eos%np,0:eos%nh),   eos%snd(0:eos%np,0:eos%nh))
    allocate(tmpblk(0:eos%np,0:eos%nh))
    ! Var 1: Pressure — i varia piu' veloce in BLOCK (I-axis=p)
    read(iu,*,iostat=ios) ((tmpblk(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    eos%plo = tmpblk(0,0);  eos%dpt = tmpblk(1,0) - tmpblk(0,0)
    ! Var 2: Enthalpy
    read(iu,*,iostat=ios) ((tmpblk(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    eos%hlo = tmpblk(0,0);  eos%dht = tmpblk(0,1) - tmpblk(0,0)
    ! Var 3: Density
    read(iu,*,iostat=ios) ((eos%rho(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 4: Temperature
    read(iu,*,iostat=ios) ((eos%T(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 5: dRho/dT (skip)
    read(iu,*,iostat=ios) ((tmpblk(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 6: dRho/dh
    read(iu,*,iostat=ios) ((eos%rh(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 7: Cp
    read(iu,*,iostat=ios) ((eos%cp(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 8: Entropy
    read(iu,*,iostat=ios) ((eos%s(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 9: dRho/dp (skip)
    read(iu,*,iostat=ios) ((tmpblk(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 10: SoundSpeed
    read(iu,*,iostat=ios) ((eos%snd(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    deallocate(tmpblk)
    close(iu)
    write(*,'(A,I0,A,I0)') '   thermo (p,h): ',eos%np+1,' x ',eos%nh+1
    write(*,'(A,1PE12.5,A,E12.5,A,E12.5)') '     p range: ', &
      eos%plo,' .. ',eos%plo+eos%np*eos%dpt,' dp=',eos%dpt
    write(*,'(A,1PE12.5,A,E12.5,A,E12.5)') '     h range: ', &
      eos%hlo,' .. ',eos%hlo+eos%nh*eos%dht,' dh=',eos%dht

    ! --- transport.dat  Tecplot BLOCK, I=p, J=h — vars: Pressure, Enthalpy, Viscosity, Conductivity
    !     Per [[ares-realfluid-thermo-layout]] must share orientation with thermo.dat.
    call find_file(folder,'transport.dat',fn,ex)
    if (.not.ex) then; write(*,'(A)') ' ERRORE: transport.dat non trovato'; ierr=2; return; endif
    open(newunit=iu,file=trim(fn),status='old',action='read',iostat=ios)
    if(ios/=0) then; ierr=2; return; endif
    call read_tec_table_header(iu, t_ni, t_nj, ios)
    if (ios/=0) then
      write(*,'(A)') ' ERRORE: header Tecplot non valido in transport.dat'
      ierr=2; close(iu); return
    endif
    if (t_ni-1 /= eos%np .or. t_nj-1 /= eos%nh) then
      write(*,'(A)') ' ERRORE: transport.dat ha griglia diversa da thermo.dat!'
      write(*,'(A,I0,A,I0,A,I0,A,I0)') '   thermo: ',eos%np+1,'x',eos%nh+1, &
        '  transport: ',t_ni,'x',t_nj
      ierr=2; close(iu); return
    endif
    allocate(eos%mu(0:eos%np,0:eos%nh), eos%kth(0:eos%np,0:eos%nh))
    allocate(tmpblk(0:eos%np,0:eos%nh))
    ! Var 1: Pressure (skip), Var 2: Enthalpy (skip)
    read(iu,*,iostat=ios) ((tmpblk(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    read(iu,*,iostat=ios) ((tmpblk(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 3: Viscosity
    read(iu,*,iostat=ios) ((eos%mu(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    ! Var 4: Conductivity
    read(iu,*,iostat=ios) ((eos%kth(i,j), i=0,eos%np), j=0,eos%nh)
    if(ios/=0) then; ierr=2; close(iu); return; endif
    deallocate(tmpblk)
    close(iu)
    write(*,'(A,1PE12.5,A,E12.5)') '   transport.dat: mu(0,0)=', &
      eos%mu(0,0),' k(0,0)=',eos%kth(0,0)

    ! --- (p,T)->h costruita in memoria invertendo eos%T per colonna di p,
    !     stessa logica di FLINT/ph2pT in Load_ThermoTransport.f90.
    call build_pT2h(eos)
    write(*,'(A,I0,A,I0)') '   pT2h (p,T) built: ',eos%np+1,' x ',eos%nh+1

    eos%ok = .true.
  end subroutine

  ! Costruisce la tabella inversa (p,T)->h a partire da eos%T(p,h).
  ! Per ogni colonna di p (indice i) la griglia T e' uniforme tra T_tab(i,0)
  ! e T_tab(i,Nh) con Nh+1 nodi.  h_pT(i,k) e' l'entalpia che produce
  ! T = Tmin2(i) + k*deltaT(i), trovata per inversione lineare di T_tab(i,:).
  subroutine build_pT2h(eos)
    type(eos_t), intent(inout) :: eos
    integer  :: i, j, k, Nh
    real(dp) :: T_target, dT_row, frac
    logical  :: found

    Nh = eos%nh
    allocate(eos%Tmin2(0:eos%np), eos%deltaT(0:eos%np))
    allocate(eos%h_pT(0:eos%np, 0:Nh))

    do i = 0, eos%np
      eos%Tmin2(i)  = eos%T(i, 0)
      eos%deltaT(i) = (eos%T(i, Nh) - eos%T(i, 0)) / real(Nh, dp)
      do k = 0, Nh
        T_target = eos%Tmin2(i) + k*eos%deltaT(i)
        found = .false.
        do j = 0, Nh - 1
          if (T_target >= eos%T(i,j) .and. T_target <= eos%T(i,j+1)) then
            dT_row = eos%T(i,j+1) - eos%T(i,j)
            if (dT_row > 0) then
              frac = (T_target - eos%T(i,j)) / dT_row
            else
              frac = 0.0_dp
            endif
            eos%h_pT(i,k) = (eos%hlo + j*eos%dht) + frac*eos%dht
            found = .true.
            exit
          endif
        end do
        if (.not.found) then
          if (k == 0)  eos%h_pT(i,k) = eos%hlo
          if (k == Nh) eos%h_pT(i,k) = eos%hlo + Nh*eos%dht
        endif
      end do
    end do
    eos%have_pT = .true.
  end subroutine

  ! =======================  INTERPOLAZIONE  ================================

  pure function bilerp(x,x0,dx,nx, y,y0,dy,ny, tab) result(val)
    real(dp), intent(in) :: x,x0,dx, y,y0,dy
    integer, intent(in)  :: nx, ny
    real(dp), intent(in) :: tab(0:nx,0:ny)
    real(dp) :: val, A,B,C,E, u,v
    integer :: i, j
    i = int((x-x0)/dx); j = int((y-y0)/dy)
    if (i<0) i=0; if (i>nx-1) i=nx-1
    if (j<0) j=0; if (j>ny-1) j=ny-1
    A = tab(i,j)
    B = (tab(i+1,j  )-tab(i,j))/dx
    C = (tab(i,  j+1)-tab(i,j))/dy
    E = (tab(i,j)-tab(i+1,j)-tab(i,j+1)+tab(i+1,j+1))/(dx*dy)
    u = x-x0-i*dx; v = y-y0-j*dy
    val = A + B*u + C*v + E*u*v
  end function

  ! lookup generico (p,h) -> proprieta'
  pure function ph(eos, tab, p, h) result(val)
    type(eos_t), intent(in) :: eos
    real(dp), intent(in) :: tab(0:eos%np,0:eos%nh), p, h
    real(dp) :: val
    val = bilerp(p,eos%plo,eos%dpt,eos%np, h,eos%hlo,eos%dht,eos%nh, tab)
  end function

  ! (p,T) -> h  (solo per Prw a parete).
  ! Griglia T per colonna di p: interpolo h sulle due colonne i, i+1 e poi in p.
  pure function pT2h(eos, p, Tv) result(h)
    type(eos_t), intent(in) :: eos
    real(dp), intent(in) :: p, Tv
    real(dp) :: h, h_i, h_ip1, up, Tk, fk
    integer :: i, k_i, k_ip1
    i = int((p - eos%plo)/eos%dpt)
    if (i < 0) i = 0
    if (i > eos%np-1) i = eos%np-1
    up = (p - (eos%plo + i*eos%dpt)) / eos%dpt
    ! colonna i
    if (eos%deltaT(i) > 0) then
      Tk = (Tv - eos%Tmin2(i)) / eos%deltaT(i)
    else
      Tk = 0.0_dp
    endif
    k_i = int(Tk)
    if (k_i < 0)         k_i = 0
    if (k_i > eos%nh-1)  k_i = eos%nh-1
    fk = Tk - real(k_i, dp)
    h_i = eos%h_pT(i,k_i) + fk*(eos%h_pT(i,k_i+1) - eos%h_pT(i,k_i))
    ! colonna i+1
    if (eos%deltaT(i+1) > 0) then
      Tk = (Tv - eos%Tmin2(i+1)) / eos%deltaT(i+1)
    else
      Tk = 0.0_dp
    endif
    k_ip1 = int(Tk)
    if (k_ip1 < 0)         k_ip1 = 0
    if (k_ip1 > eos%nh-1)  k_ip1 = eos%nh-1
    fk = Tk - real(k_ip1, dp)
    h_ip1 = eos%h_pT(i+1,k_ip1) + fk*(eos%h_pT(i+1,k_ip1+1) - eos%h_pT(i+1,k_ip1))
    h = h_i + up*(h_ip1 - h_i)
  end function

  ! ========================  EXTRACT_1D  ===================================

  subroutine extract_1d(field_file, wall_file, output_file, input_folder, ini_file)
    character(len=*), intent(in) :: field_file, wall_file, output_file, input_folder
    character(len=*), intent(in) :: ini_file

    type(tec_file_t) :: fld, wall
    type(eos_t) :: eos
    real(dp) :: hs_ks
    integer :: ierr, iu, ios, ib, ic, jc, kc, iter
    integer :: ni,nj,nk, nx,ny,nz
    logical :: have_wall

    ! indici variabili campo
    integer :: ip, ih, iu_f, iv_f, iw_f
    integer :: imil, ikl, irho, iT_f, isnd  ! variabili opzionali da campo
    ! indici variabili parete (it_sc = magnitudine scalare tau come fallback)
    integer :: it1, it2, it3, it_sc, iTw, iqw, ihs
    ! zone wall
    integer :: iw3,iw4,iw5,iw6

    real(dp) :: Area, Vol, portata, ds, slung, sign_i
    real(dp) :: p_sum, h_sum, h0_sum, T_sum
    real(dp) :: u_mass, v_mass, w_mass
    real(dp) :: p_b, h0_b, h_b, G
    real(dp) :: d_b, T_b, cp_b, s_b, a_b, mu_b, k_b
    real(dp) :: u_b, Tmean
    real(dp) :: p0_b, T0_b, Prw
    real(dp) :: Per, Diam, Re_b, Pr_b
    ! parete
    real(dp) :: Aw3,taus3,qs3,Tsm3, Aw4,taus4,qs4,Tsm4
    real(dp) :: Aw5,taus5,qs5,Tsm5, Aw6,taus6,qs6,Tsm6
    real(dp) :: Awt, tauw3,tauw4,tauw5,tauw6,tauw_a, qw3,qw4,qw5,qw6,qw_a, Twm3,Twm4,Twm5,Twm6,Twm_a
    ! rugosita' sand-grain per faccia (somma pesata area, media, aggregato), da wall.tec
    real(dp) :: hss3,hss4,hss5,hss6, hsf3,hsf4,hsf5,hsf6, hs_w
    ! geometria
    real(dp) :: pf1(3),pf2(3),pf3(3),pf4(3), Af, nf(3), un
    real(dp) :: n0(3),n1(3),n2(3),n3(3),n4(3),n5(3),n6(3),n7(3)
    real(dp) :: xm, pc, hc, uc, vc, wc, rc
    real(dp) :: cr(3), dir(3), ctr1(3), ctr2(3)
    ! Newton
    real(dp) :: F, dFdh, rho_h, rh_h, err_n, sp
    real(dp) :: hw, cpw, muw, kw
    ! periodicità settoriale
    real(dp) :: p_fac, nf_ref(3), C_cent(3), v1_pf(3), v2_pf(3), cp_pf(3), cn_pf, dn_pf, delta_phi

    real(dp), parameter :: tol=1e-5_dp, damp=0.5_dp

    ! ---------- 1) Lettura tabelle EOS ----------
    write(*,'(2A)') ' Tabelle EOS da: ', trim(input_folder)
    call read_eos_tables(input_folder, eos, ierr)
    if (ierr/=0 .or. .not.eos%ok) then
      write(*,'(A)') ' ERRORE: tabelle EOS non caricate.'; return
    endif

    ! ---------- 1b) Rugosita' equivalente sand-grain da input.ini ----------
    call read_ks_from_ini(ini_file, hs_ks)
    write(*,'(A,1PE12.5,2A)') '   hs (ks) = ', hs_ks, ' m  da ', trim(ini_file)

    ! ---------- 2) Lettura field.tec ----------
    write(*,'(2A)') ' Lettura: ', trim(field_file)
    call read_tec_ascii(field_file, fld, ierr)
    if (ierr/=0) then; write(*,'(A)') ' ERRORE lettura field.'; return; endif
    write(*,'(A,I0,A,I0,A)') '   ',fld%nzones,' blocchi, ',fld%nvar,' variabili'

    ! Diagnostica nomi variabili campo
    write(*,'(A)') '   Variabili campo:'
    do iter=1,fld%nvar
      write(*,'(A,I2,A,A,A)') '     ',iter,': "',trim(fld%varname(iter)),'"'
    end do

    ip   = find_var(fld,'p')
    ih   = find_var(fld,'h')
    iu_f = find_var(fld,'u')
    iv_f = find_var(fld,'v')
    iw_f = find_var(fld,'w')
    write(*,'(A,5I4)') '   Indici p,h,u,v,w: ',ip,ih,iu_f,iv_f,iw_f

    if (ip==0.or.ih==0.or.iu_f==0.or.iv_f==0.or.iw_f==0) then
      write(*,'(A)') ' ERRORE: variabili (p,h,u,v,w) non trovate.'; return
    endif

    ! Variabili opzionali gia' calcolate dal solutore (se presenti, usate al posto delle tabelle EOS)
    imil = find_var(fld,'mil');  if(imil==0) imil=find_var(fld,'visc'); if(imil==0) imil=find_var(fld,'mi')
    ikl  = find_var(fld,'kl');   if(ikl ==0) ikl =find_var(fld,'cond'); if(ikl ==0) ikl =find_var(fld,'k')
    irho = find_var(fld,'rho');  if(irho==0) irho=find_var(fld,'d')
    iT_f = find_var(fld,'T');    isnd = find_var(fld,'sound')
    write(*,'(A,5I4)') '   Indici mil,kl,rho,T,sound: ',imil,ikl,irho,iT_f,isnd

    ! ---------- 3) Lettura wall.tec (opzionale) ----------
    have_wall=.false.; inquire(file=trim(wall_file),exist=have_wall)
    if (have_wall) then
      write(*,'(2A)') ' Lettura: ', trim(wall_file)
      call read_tec_ascii(wall_file, wall, ierr)
      if (ierr/=0) then; have_wall=.false.
      else
        write(*,'(A,I0,A)') '   ',wall%nzones,' zone wall trovate:'
        do iter=1,wall%nzones
          write(*,'(A,I2,A,A,A,3(A,I0))') '     ',iter,': "', &
            trim(wall%zone(iter)%name),'", dims=', &
            'I=',wall%zone(iter)%ni,',J=',wall%zone(iter)%nj, &
            ',K=',wall%zone(iter)%nk
        end do
        write(*,'(A)') '   Variabili wall:'
        do iter=1,wall%nvar
          write(*,'(A,I2,A,A,A)') '     ',iter,': "',trim(wall%varname(iter)),'"'
        end do
        ! Cerca tauX con tutti i nomi possibili (diverse versioni di ARES)
        it1=find_var(wall,'TAUWx'); if(it1==0) it1=find_var(wall,'tauWx')
        if(it1==0) it1=find_var(wall,'taux');  if(it1==0) it1=find_var(wall,'tau_x')
        if(it1==0) it1=find_var(wall,'shx');   if(it1==0) it1=find_var(wall,'twx')
        it2=find_var(wall,'TAUWy'); if(it2==0) it2=find_var(wall,'tauWy')
        if(it2==0) it2=find_var(wall,'tauy');  if(it2==0) it2=find_var(wall,'tau_y')
        if(it2==0) it2=find_var(wall,'shy');   if(it2==0) it2=find_var(wall,'twy')
        it3=find_var(wall,'TAUWz'); if(it3==0) it3=find_var(wall,'tauWz')
        if(it3==0) it3=find_var(wall,'tauz');  if(it3==0) it3=find_var(wall,'tau_z')
        if(it3==0) it3=find_var(wall,'shz');   if(it3==0) it3=find_var(wall,'twz')
        ! Fallback: grandezza scalare del modulo di tau (se non trovati i vettori)
        it_sc=0
        if(it1==0.or.it2==0.or.it3==0) then
          it_sc=find_var(wall,'tauw');  if(it_sc==0) it_sc=find_var(wall,'tau')
          if(it_sc==0) it_sc=find_var(wall,'tauW'); if(it_sc==0) it_sc=find_var(wall,'tau_w')
          if(it_sc==0) it_sc=find_var(wall,'tw');   if(it_sc==0) it_sc=find_var(wall,'|tau|')
        endif
        iTw=find_var(wall,'Tw');   if(iTw==0) iTw=find_var(wall,'T_wall')
        iqw=find_var(wall,'qw');   if(iqw==0) iqw=find_var(wall,'q_wall')
        ! Rugosita' sand-grain di parete: ultima variabile in wall.tec ("hs").
        ! Se assente (ihs==0) la rugosita' resta 0 (per ogni parete).
        ihs=find_var(wall,'hs');   if(ihs==0) ihs=find_var(wall,'ks')
        write(*,'(A,7I4)') '   Indici tauX,tauY,tauZ,tau_sc,Tw,qw,hs: ',it1,it2,it3,it_sc,iTw,iqw,ihs
      endif
    endif

    ! ---------- 4) Apertura output ----------
    open(newunit=iu,file=trim(output_file),status='replace',action='write',iostat=ios)
    if (ios/=0) then; write(*,'(2A)') ' ERRORE output: ',trim(output_file); return; endif
    write(iu,'(A)') 'title = "1d flow data"'
    write(iu,'(A)') 'variables = "x[m]","s[m]","A[m^2]","P[m]","Dh[m]",' // &
      '"G[kg/s/m^2]","U[m/s]",' // &
      '"p[Pa]","p0[Pa]","T[K]","T0[K]","rho[kg/m^3]","h[J/kg]",' // &
      '"s[J/kg K]","Prw","cp[J/kg K]","mi[Pa/s]","k[W/m K]","Tm","Re",' // &
      '"Pr","tau3[Pa]","tau4[Pa]","tau5[Pa]","tau6[Pa]","tau[Pa]","qw3[W/m^2]"' // &
      ',"qw4[W/m^2]","qw5[W/m^2]","qw6[W/m^2]","qw[W/m^2]","Tw3[K]","Tw4[' // &
      'K]","Tw5[K]","Tw6[K]","Tw[K]","hs[m]"'

    ! ---------- 5) Loop blocchi ----------
    do ib=1,fld%nzones
      ni=fld%zone(ib)%ni; nj=fld%zone(ib)%nj; nk=fld%zone(ib)%nk
      nx=fld%zone(ib)%nx; ny=fld%zone(ib)%ny; nz=fld%zone(ib)%nz
      write(*,'(A,I0,A,I0,A,I0,A,I0)') '   Blocco ',ib,': ',nx,' x',ny,' x',nz
      write(iu,'(A,I0,A,I0,A)') 'ZONE T="B',ib,'", I=',nx,', F=POINT'

      ! Zone wall
      iw3=0; iw4=0; iw5=0; iw6=0
      if (have_wall) then
        iw3=find_wall_zone(wall,ib,3); iw4=find_wall_zone(wall,ib,4)
        iw5=find_wall_zone(wall,ib,5); iw6=find_wall_zone(wall,ib,6)
        write(*,'(A,4I4)') '   Wall zones (3,4,5,6): ',iw3,iw4,iw5,iw6
      endif

      ! Segno normale I-face (convenzione ARES: verso i crescenti)
      pf1=(/fld%zone(ib)%xn(1,1,1),fld%zone(ib)%yn(1,1,1),fld%zone(ib)%zn(1,1,1)/)
      pf2=(/fld%zone(ib)%xn(1,2,1),fld%zone(ib)%yn(1,2,1),fld%zone(ib)%zn(1,2,1)/)
      pf3=(/fld%zone(ib)%xn(1,2,2),fld%zone(ib)%yn(1,2,2),fld%zone(ib)%zn(1,2,2)/)
      pf4=(/fld%zone(ib)%xn(1,1,2),fld%zone(ib)%yn(1,1,2),fld%zone(ib)%zn(1,1,2)/)
      cr(1)=0.5_dp*((pf3(2)-pf1(2))*(pf4(3)-pf2(3))-(pf3(3)-pf1(3))*(pf4(2)-pf2(2)))
      cr(2)=0.5_dp*((pf3(3)-pf1(3))*(pf4(1)-pf2(1))-(pf3(1)-pf1(1))*(pf4(3)-pf2(3)))
      cr(3)=0.5_dp*((pf3(1)-pf1(1))*(pf4(2)-pf2(2))-(pf3(2)-pf1(2))*(pf4(1)-pf2(1)))
      ctr1=0.25_dp*(pf1+pf2+pf3+pf4)
      ctr2=0.25_dp*( &
        (/fld%zone(ib)%xn(2,1,1),fld%zone(ib)%yn(2,1,1),fld%zone(ib)%zn(2,1,1)/) &
       +(/fld%zone(ib)%xn(2,2,1),fld%zone(ib)%yn(2,2,1),fld%zone(ib)%zn(2,2,1)/) &
       +(/fld%zone(ib)%xn(2,2,2),fld%zone(ib)%yn(2,2,2),fld%zone(ib)%zn(2,2,2)/) &
       +(/fld%zone(ib)%xn(2,1,2),fld%zone(ib)%yn(2,1,2),fld%zone(ib)%zn(2,1,2)/))
      dir = ctr2-ctr1
      sign_i = sign(1.0_dp, cr(1)*dir(1)+cr(2)*dir(2)+cr(3)*dir(3))

      ! == Fattore di periodicita' settoriale (come TWODAX in wrt_TECfile.F) ==
      ! delthe = angolo del settore K = 2*pi/p_fac.
      ! Centro = media dei nodi sull'ASSE INTERNO (jc=1) all'I-face ic=1.
      ! (usare il centroide di tutti i nodi darebbe un centro a meta' settore,
      !  raddoppiando l'angolo apparente e dimezzando p_fac.)
      cn_pf = sqrt(cr(1)**2+cr(2)**2+cr(3)**2)
      if (cn_pf > 0) then; nf_ref = cr/cn_pf; else; nf_ref = (/1.0_dp,0.0_dp,0.0_dp/); endif
      ! Asse di simmetria = media dei nodi jc=1 (bordo interno / asse) sull'I-face ic=1
      C_cent = 0.0_dp
      do kc=1,nk
        C_cent = C_cent + (/fld%zone(ib)%xn(1,1,kc),fld%zone(ib)%yn(1,1,kc),fld%zone(ib)%zn(1,1,kc)/)
      end do
      C_cent = C_cent / real(nk, dp)
      ! Vettori dall'asse ai due bordi K della parete esterna (jc=nj)
      v1_pf = (/fld%zone(ib)%xn(1,nj,1), fld%zone(ib)%yn(1,nj,1), fld%zone(ib)%zn(1,nj,1)/) - C_cent
      v2_pf = (/fld%zone(ib)%xn(1,nj,nk),fld%zone(ib)%yn(1,nj,nk),fld%zone(ib)%zn(1,nj,nk)/) - C_cent
      ! Proiezione sul piano della sezione (rimozione componente assiale)
      v1_pf = v1_pf - dot_product(v1_pf, nf_ref)*nf_ref
      v2_pf = v2_pf - dot_product(v2_pf, nf_ref)*nf_ref
      ! Angolo tra i due vettori = angolo del settore
      cp_pf(1) = v1_pf(2)*v2_pf(3)-v1_pf(3)*v2_pf(2)
      cp_pf(2) = v1_pf(3)*v2_pf(1)-v1_pf(1)*v2_pf(3)
      cp_pf(3) = v1_pf(1)*v2_pf(2)-v1_pf(2)*v2_pf(1)
      cn_pf = sqrt(dot_product(cp_pf,cp_pf))
      dn_pf = dot_product(v1_pf, v2_pf)
      if (cn_pf > 1e-15_dp .or. abs(dn_pf) > 1e-15_dp) then
        delta_phi = atan2(cn_pf, dn_pf)   ! angolo settore [rad], in (0, pi]
      else
        delta_phi = 0.0_dp
      endif
      if (delta_phi > 1e-8_dp) then
        p_fac = 2.0_dp*acos(-1.0_dp) / delta_phi
      else
        p_fac = 1.0_dp
      endif
      p_fac = max(p_fac, 1.0_dp)
      write(*,'(A,1PE10.3,A,E10.3)') '   Settore: angolo=',delta_phi,' rad  p_fac=',p_fac

      slung = 0; h_b = 0

      ! ---------- Loop stazioni I ----------
      do ic=1,nx

        ! == A) Accumulo sulla sezione ==
        Area=0; Vol=0; portata=0
        p_sum=0; h_sum=0; h0_sum=0; T_sum=0
        u_mass=0; v_mass=0; w_mass=0

        do kc=1,nz; do jc=1,ny
          ! I-face (sinistra cella ic)
          pf1=(/fld%zone(ib)%xn(ic,jc,  kc  ),fld%zone(ib)%yn(ic,jc,  kc  ),fld%zone(ib)%zn(ic,jc,  kc  )/)
          pf2=(/fld%zone(ib)%xn(ic,jc+1,kc  ),fld%zone(ib)%yn(ic,jc+1,kc  ),fld%zone(ib)%zn(ic,jc+1,kc  )/)
          pf3=(/fld%zone(ib)%xn(ic,jc+1,kc+1),fld%zone(ib)%yn(ic,jc+1,kc+1),fld%zone(ib)%zn(ic,jc+1,kc+1)/)
          pf4=(/fld%zone(ib)%xn(ic,jc,  kc+1),fld%zone(ib)%yn(ic,jc,  kc+1),fld%zone(ib)%zn(ic,jc,  kc+1)/)
          call quad_area_normal(pf1,pf2,pf3,pf4, Af, nf)
          nf = nf*sign_i   ! convenzione ARES
          Area = Area + Af

          ! Volume esaedro
          n0=(/fld%zone(ib)%xn(ic,  jc,  kc  ),fld%zone(ib)%yn(ic,  jc,  kc  ),fld%zone(ib)%zn(ic,  jc,  kc  )/)
          n1=(/fld%zone(ib)%xn(ic+1,jc,  kc  ),fld%zone(ib)%yn(ic+1,jc,  kc  ),fld%zone(ib)%zn(ic+1,jc,  kc  )/)
          n2=(/fld%zone(ib)%xn(ic+1,jc+1,kc  ),fld%zone(ib)%yn(ic+1,jc+1,kc  ),fld%zone(ib)%zn(ic+1,jc+1,kc  )/)
          n3=(/fld%zone(ib)%xn(ic,  jc+1,kc  ),fld%zone(ib)%yn(ic,  jc+1,kc  ),fld%zone(ib)%zn(ic,  jc+1,kc  )/)
          n4=(/fld%zone(ib)%xn(ic,  jc,  kc+1),fld%zone(ib)%yn(ic,  jc,  kc+1),fld%zone(ib)%zn(ic,  jc,  kc+1)/)
          n5=(/fld%zone(ib)%xn(ic+1,jc,  kc+1),fld%zone(ib)%yn(ic+1,jc,  kc+1),fld%zone(ib)%zn(ic+1,jc,  kc+1)/)
          n6=(/fld%zone(ib)%xn(ic+1,jc+1,kc+1),fld%zone(ib)%yn(ic+1,jc+1,kc+1),fld%zone(ib)%zn(ic+1,jc+1,kc+1)/)
          n7=(/fld%zone(ib)%xn(ic,  jc+1,kc+1),fld%zone(ib)%yn(ic,  jc+1,kc+1),fld%zone(ib)%zn(ic,  jc+1,kc+1)/)
          Vol = Vol + hex_vol(n0,n1,n2,n3,n4,n5,n6,n7)

          ! Variabili primitive della cella
          pc = fld%zone(ib)%v(ip,  ic,jc,kc)
          hc = fld%zone(ib)%v(ih,  ic,jc,kc)
          uc = fld%zone(ib)%v(iu_f,ic,jc,kc)
          vc = fld%zone(ib)%v(iv_f,ic,jc,kc)
          wc = fld%zone(ib)%v(iw_f,ic,jc,kc)

          ! rho: usa direttamente dal campo se disponibile, altrimenti EOS
          if (irho>0) then
            rc = fld%zone(ib)%v(irho,ic,jc,kc)
          else
            rc = ph(eos, eos%rho, pc, hc)
          endif

          un = uc*nf(1)+vc*nf(2)+wc*nf(3)

          ! Accumulo
          portata = portata + rc*un*Af
          u_mass  = u_mass  + rc*un*uc*Af
          v_mass  = v_mass  + rc*un*vc*Af
          w_mass  = w_mass  + rc*un*wc*Af
          p_sum   = p_sum   + pc*Af
          h_sum   = h_sum   + hc*Af
          h0_sum  = h0_sum  + (hc+0.5_dp*(uc**2+vc**2+wc**2))*rc*un*Af
          ! T: usa campo se disponibile, altrimenti EOS
          if (iT_f>0) then
            T_sum = T_sum + fld%zone(ib)%v(iT_f,ic,jc,kc)*Af
          else
            T_sum = T_sum + ph(eos,eos%T,pc,hc)*Af
          endif
        end do; end do

        ! Diagnostica prima stazione
        if (ic==1 .and. ib==1) then
          write(*,'(A)') '   --- STAZIONE 1 DIAGNOSTICA ---'
          write(*,'(A,1PE12.5,A,E12.5,A,E12.5)') '   Area/Vol/portata: ',Area,', ',Vol,', ',portata
          write(*,'(A,3(1PE12.5))') '   n(1,1,1): ', &
            fld%zone(ib)%xn(1,1,1),fld%zone(ib)%yn(1,1,1),fld%zone(ib)%zn(1,1,1)
          write(*,'(A,3(1PE12.5))') '   n(1,2,1): ', &
            fld%zone(ib)%xn(1,2,1),fld%zone(ib)%yn(1,2,1),fld%zone(ib)%zn(1,2,1)
          write(*,'(A,5(1PE12.5))') '   cella(1,1,1) p,h,u,v,w: ', &
            fld%zone(ib)%v(ip,1,1,1),fld%zone(ib)%v(ih,1,1,1), &
            fld%zone(ib)%v(iu_f,1,1,1),fld%zone(ib)%v(iv_f,1,1,1),fld%zone(ib)%v(iw_f,1,1,1)
        endif

        ! Normalizzazione
        ds = Vol/Area
        if (ic==1) then; slung=0.5_dp*ds; else; slung=slung+ds; endif
        p_b   = p_sum/Area
        h0_b  = h0_sum/portata
        Tmean = T_sum/Area
        G     = portata/Area

        ! == B) Newton per h_bulk ==
        ! h_b + 0.5*(G/rho(p_b,h_b))^2 = h0_b
        if (ic==1) h_b = h_sum/Area   ! guess iniziale: media area di h
        do iter=1,200
          rho_h = ph(eos, eos%rho, p_b, h_b)
          rh_h  = ph(eos, eos%rh,  p_b, h_b)
          F     = h_b + 0.5_dp*(G/rho_h)**2 - h0_b
          dFdh  = 1.0_dp - G**2*rh_h/rho_h**3
          err_n = abs(F/dFdh)/max(abs(h_b),1.0_dp)
          h_b   = h_b - damp*F/dFdh
          if (err_n < tol) exit
        end do

        ! == C) Grandezze bulk da tabelle (p_b, h_b) ==
        d_b  = ph(eos, eos%rho, p_b, h_b)
        T_b  = ph(eos, eos%T,   p_b, h_b)
        cp_b = ph(eos, eos%cp,  p_b, h_b)
        s_b  = ph(eos, eos%s,   p_b, h_b)
        a_b  = ph(eos, eos%snd, p_b, h_b)
        u_b  = G/d_b

        ! mu_b e k_b dalle tabelle EOS valutate alle grandezze bulk (p_b, h_b)
        mu_b = ph(eos, eos%mu,  p_b, h_b)
        k_b  = ph(eos, eos%kth, p_b, h_b)

        ! Diagnostica bulk prima stazione
        if (ic==1 .and. ib==1) then
          write(*,'(A,5(1PE12.5))') '   bulk p,h,G,d,T: ',p_b,h_b,G,d_b,T_b
          write(*,'(A,4(1PE12.5),I4)')  '   mu,k,cp,a,Newton: ',mu_b,k_b,cp_b,a_b,iter
          write(*,'(A,4(1PE12.5))') '   Per,Diam: (dopo D)'
        endif

        ! == D) Parete ==
        Aw3=0;taus3=0;qs3=0;Tsm3=0; Aw4=0;taus4=0;qs4=0;Tsm4=0
        Aw5=0;taus5=0;qs5=0;Tsm5=0; Aw6=0;taus6=0;qs6=0;Tsm6=0
        hss3=0;hss4=0;hss5=0;hss6=0
        if (have_wall) then
          if(iw3>0) call wall_station_sums(wall%zone(iw3),3,ic,it1,it2,it3,it_sc,iTw,iqw,ihs, Aw3,taus3,qs3,Tsm3,hss3)
          if(iw4>0) call wall_station_sums(wall%zone(iw4),4,ic,it1,it2,it3,it_sc,iTw,iqw,ihs, Aw4,taus4,qs4,Tsm4,hss4)
          if(iw5>0) call wall_station_sums(wall%zone(iw5),5,ic,it1,it2,it3,it_sc,iTw,iqw,ihs, Aw5,taus5,qs5,Tsm5,hss5)
          if(iw6>0) call wall_station_sums(wall%zone(iw6),6,ic,it1,it2,it3,it_sc,iTw,iqw,ihs, Aw6,taus6,qs6,Tsm6,hss6)
        endif

        Awt = Aw3+Aw4+Aw5+Aw6
        tauw3=0;tauw4=0;tauw5=0;tauw6=0;tauw_a=0
        qw3=0;qw4=0;qw5=0;qw6=0;qw_a=0
        Twm3=0;Twm4=0;Twm5=0;Twm6=0;Twm_a=0
        hsf3=0;hsf4=0;hsf5=0;hsf6=0; hs_w=0
        if(Aw3>0) then; tauw3=taus3/Aw3; qw3=qs3/Aw3; Twm3=Tsm3/Aw3; hsf3=hss3/Aw3; endif
        if(Aw4>0) then; tauw4=taus4/Aw4; qw4=qs4/Aw4; Twm4=Tsm4/Aw4; hsf4=hss4/Aw4; endif
        if(Aw5>0) then; tauw5=taus5/Aw5; qw5=qs5/Aw5; Twm5=Tsm5/Aw5; hsf5=hss5/Aw5; endif
        if(Aw6>0) then; tauw6=taus6/Aw6; qw6=qs6/Aw6; Twm6=Tsm6/Aw6; hsf6=hss6/Aw6; endif
        if(Awt>0) then
          tauw_a=(taus3+taus4+taus5+taus6)/Awt
          qw_a=(qs3+qs4+qs5+qs6)/Awt
          Twm_a=(Tsm3+Tsm4+Tsm5+Tsm6)/Awt
        endif
        ! Rugosita' di parete rappresentativa della stazione: la parete piu' rugosa
        ! tra le facce presenti (le pareti senza dato 'hs' contribuiscono 0).
        hs_w = max(hsf3,hsf4,hsf5,hsf6)

        Per=0; Diam=0
        if (ds>0) Per=Awt/ds
        if (Per>0) Diam=4*Area/Per

        ! == E) Re, Pr, Prw ==
        Re_b=0; Pr_b=0; Prw=0
        if (Diam>0 .and. mu_b>0) Re_b = G*Diam/mu_b
        if (mu_b>0 .and. k_b>0)  Pr_b = mu_b*cp_b/k_b

        ! Diagnostica parete/Re/Pr prima stazione
        if (ic==1 .and. ib==1) then
          write(*,'(A,4(1PE12.5))') '   Aw3456: ',Aw3,Aw4,Aw5,Aw6
          write(*,'(A,4(1PE12.5))') '   Awt,ds,Per,Diam: ',Awt,ds,Per,Diam
          write(*,'(A,2(1PE12.5))') '   Re,Pr: ',Re_b,Pr_b
          write(*,'(A)') '   ---------------------------------'
        endif

        ! Prw: proprieta' a (p_b, h_wall) dove h_wall=pT2h(p_b, Twm_a).
        if (Twm_a>0 .and. eos%have_pT) then
          hw  = pT2h(eos, p_b, Twm_a)
          cpw = ph(eos, eos%cp,  p_b, hw)
          muw = ph(eos, eos%mu,  p_b, hw)
          kw  = ph(eos, eos%kth, p_b, hw)
          if (kw>0) Prw = muw*cpw/kw
        endif

        ! == F) Condizioni totali ==
        ! h0 e' noto (h0_b). Cerco p0 tale che s(p0, h0_b) = s_b  [Newton 1D]
        p0_b = p_b + 0.5_dp*d_b*u_b**2   ! guess
        do iter=1,500
          F    = ph(eos, eos%s, p0_b, h0_b) - s_b
          sp   = (ph(eos, eos%s, p0_b*(1+tol), h0_b) &
               -  ph(eos, eos%s, p0_b,         h0_b)) / (p0_b*tol)
          if (abs(sp)<1e-30_dp) exit
          err_n = abs(F/sp)/max(p0_b,1.0_dp)
          p0_b  = p0_b - damp*F/sp
          if (err_n < tol) exit
        end do
        T0_b = ph(eos, eos%T, p0_b, h0_b)

        ! == G) Coordinata assiale ==
        xm = 0.5_dp*(fld%zone(ib)%xn(ic,1,1)+fld%zone(ib)%xn(ic+1,1,1))

        ! == H) Scrittura (37 colonne, formato wrt_1d) ==
        write(iu,'(1P,37E14.6)') &
          xm, slung, Area*p_fac, Per*p_fac, Diam, &  !  1-5
          G, u_b, &                            !  6-7
          p_b, p0_b, T_b, T0_b, d_b, h_b, &  !  8-13
          s_b, Prw, cp_b, mu_b, k_b, Tmean, & ! 14-19
          Re_b, Pr_b, &                        ! 20-21
          tauw3,tauw4,tauw5,tauw6,tauw_a, &    ! 22-26
          qw3,qw4,qw5,qw6,qw_a, &            ! 27-31
          Twm3,Twm4,Twm5,Twm6,Twm_a, &       ! 32-36
          hs_w                                 ! 37: hs (rugosita' da wall.tec; 0 se assente)

      end do ! ic
    end do ! ib

    close(iu)
    write(*,'(2A)') ' Scritto: ', trim(output_file)
  end subroutine extract_1d

end module Mod_Extract1D


! ===========================  MAIN  =======================================
program Extract1D_main
  use Mod_Extract1D
  implicit none
  character(len=256) :: f_field, f_wall, f_out, f_input, f_ini

  f_field='OUTPUT/field.tec'; f_wall='OUTPUT/wall.tec'
  f_out='OUTPUT/1d.dat';      f_input='INPUT'
  f_ini='input.ini'

  if (command_argument_count()>=1) call get_command_argument(1,f_field)
  if (command_argument_count()>=2) call get_command_argument(2,f_wall)
  if (command_argument_count()>=3) call get_command_argument(3,f_out)
  if (command_argument_count()>=4) call get_command_argument(4,f_input)
  if (command_argument_count()>=5) call get_command_argument(5,f_ini)

  call extract_1d(trim(f_field), trim(f_wall), trim(f_out), trim(f_input), trim(f_ini))
end program Extract1D_main
