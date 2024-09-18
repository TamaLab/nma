      module diagarpack
c
c Diagonalization of hessian matrix using ARPACK
c
c Copyright (c) 2024 Florence Tama, Osamu Miyashita
c Developed at Nagoya University, Nagoya, Japan,
c              RIKEN, Kobe, Japan
c
c References:
c
c https://www.arpack.org/home
c  
c Danny Sorensen at sorensen@caam.rice.edu
c Richard Lehoucq at rblehou@sandia.gov
c Chao Yang at cyang@lbl.gov
c Kristi Maschhoff at kristyn@tera.com
c           
      implicit none
      
      contains

      subroutine do_diag(element_i, element_j, element_s,
     .     n_nz, n, nev, ncv, tol, v, d, nconv)

      implicit none

      integer, intent(in) :: element_i(n_nz)   ! index i of sparse matrix element
      integer, intent(in) :: element_j(n_nz)   ! index j of sparse matrix element
      real*8, intent(in) :: element_s(n_nz)   ! value of sparse matrix element i, j

      integer          iparam(11), ipntr(11)

      integer, intent(in) :: n_nz   ! number of sparse matrix elements
      integer, intent(in) :: n   ! dimension of sparse matrix, i.e., maximum index
      integer, intent(in) :: ncv   ! number of requested eigen vectors
      integer, intent(in) :: nev
      real*8, intent(in) :: tol

      real*8, intent(out) :: v(n, ncv)
      real*8, intent(out) :: d(ncv, 2)

c     
c     %-----------------%
c     | Local Variables |
c     %-----------------%
c
      character        bmat*1, which*2
      integer          ido, lworkl, info, ierr,
     &                 j, ishfts, maxitr, mode1, nconv
      logical          rvec
      Real*8           sigma

      real*8 :: workl(1:ncv*(ncv+8)) ! Declare the array as allocatable
      real*8 :: workd(1:3*n)  ! Declare the array as allocatable
      real*8 :: ax(n)   ! storage for vector A*x
      real*8 :: resid(n)   ! storage for residual vector
      logical :: select(ncv)   ! used as a workspace

c
c     %------------%
c     | Parameters |
c     %------------%
c
      Real*8           zero
      parameter        (zero = 0.0D+0)
c  
c     %-----------------------------%
c     | BLAS & LAPACK routines used |
c     %-----------------------------%
c
      Real*8           dnrm2
      external         dnrm2, daxpy
c
c     %--------------------%
c     | Intrinsic function |
c     %--------------------%
c
      intrinsic        abs
c timing function not standard fortan77 and really serves no function, removed.
c clb3
c//      intrinsic        cpu_time //c
c//      real*8           sec  //c

c     %-------------------------------------------------%
c     | The following include statement and assignments |
c     | initiate trace output from the internal         |
c     | actions of ARPACK.  See debug.doc in the        |
c     | DOCUMENTS directory for usage.  Initially, the  |
c     | most useful information will be a breakdown of  |
c     | time spent in the various stages of computation |
c     | given by setting msaupd = 1.                    |
c     %-------------------------------------------------%
c
      include 'debug.h'
      ndigit = -3
      logfil = 6
      msgets = 0
      msaitr = 0 
      msapps = 0
      msaupd = 1
      msaup2 = 0
      mseigt = 0
      mseupd = 0
      bmat  = 'I'
      which = 'SM'

c
c     %-----------------------------------------------------%
c     |                                                     |
c     | Specification of stopping rules and initial         |
c     | conditions before calling DSAUPD                    |
c     |                                                     |
c     | TOL  determines the stopping criterion.             |
c     |                                                     |
c     |      Expect                                         |
c     |           abs(lambdaC - lambdaT) < TOL*abs(lambdaC) |
c     |               computed   true                       |
c     |                                                     |
c     |      If TOL .le. 0,  then TOL <- macheps            |
c     |           (machine precision) is used.              |
c     |                                                     |
c     | IDO  is the REVERSE COMMUNICATION parameter         |
c     |      used to specify actions to be taken on return  |
c     |      from DSAUPD. (See usage below.)                |
c     |                                                     |
c     |      It MUST initially be set to 0 before the first |
c     |      call to DSAUPD.                                | 
c     |                                                     |
c     | INFO on entry specifies starting vector information |
c     |      and on return indicates error codes            |
c     |                                                     |
c     |      Initially, setting INFO=0 indicates that a     | 
c     |      random starting vector is requested to         |
c     |      start the ARNOLDI iteration.  Setting INFO to  |
c     |      a nonzero value on the initial call is used    |
c     |      if you want to specify your own starting       |
c     |      vector (This vector must be placed in RESID.)  | 
c     |                                                     |
c     | The work array WORKL is used in DSAUPD as           | 
c     | workspace.  Its dimension LWORKL is set as          |
c     | illustrated below.                                  |
c     |                                                     |
c     %-----------------------------------------------------%
c
      lworkl = ncv*(ncv+8)

      info = 0
      ido = 0
c
c     %---------------------------------------------------%
c     | Specification of Algorithm Mode:                  |
c     |                                                   |
c     | This program uses the exact shift strategy        |
c     | (indicated by setting PARAM(1) = 1).              |
c     | IPARAM(3) specifies the maximum number of Arnoldi |
c     | iterations allowed.  Mode 1 of DSAUPD is used     |
c     | (IPARAM(7) = 1). All these options can be changed |
c     | by the user. For details see the documentation in |
c     | DSAUPD.                                           |
c     %---------------------------------------------------%
c
      ishfts = 1
      maxitr = 300 
      mode1  = 1
c
      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode1
c
c     %------------------------------------------------%
c     | M A I N   L O O P (Reverse communication loop) |
c     %------------------------------------------------%
c
 10   continue
c
c        %---------------------------------------------%
c        | Repeatedly call the routine DSAUPD and take | 
c        | actions indicated by parameter IDO until    |
c        | either convergence is indicated or maxitr   |
c        | has been exceeded.                          |
c        %---------------------------------------------%
c
      call dsaupd ( ido, bmat, n, which, nev, tol, resid, 
     &     ncv, v, n, iparam, ipntr, workd, workl,
     &     lworkl, info )

c
         if (ido .eq. -1 .or. ido .eq. 1) then
c
c           %--------------------------------------%
c           | Perform matrix vector multiplication |
c           |              y <--- OP*x             |
c           | The user should supply his/her own   |
c           | matrix vector multiplication routine |
c           | here that takes workd(ipntr(1)) as   |
c           | the input, and return the result to  |
c           | workd(ipntr(2)).                     |
c           %--------------------------------------%
c

            call av(n, workd(ipntr(1)), workd(ipntr(2)),
     $           n_nz, element_i, element_j, element_s)

c
c           %-----------------------------------------%
c           | L O O P   B A C K to call DSAUPD again. |
c           %-----------------------------------------%
c
            go to 10
c
         end if 
c
c     %----------------------------------------%
c     | Either we have convergence or there is |
c     | an error.                              |
c     %----------------------------------------%
c
      if ( info .lt. 0 ) then
c
c        %--------------------------%
c        | Error message. Check the |
c        | documentation in DSAUPD. |
c        %--------------------------%
c
         print *, ' '
         print *, ' Error with _saupd, info = ', info
         print *, ' Check documentation in _saupd '
         print *, ' '
c
      else 
c
c        %-------------------------------------------%
c        | No fatal errors occurred.                 |
c        | Post-Process using DSEUPD.                |
c        |                                           |
c        | Computed eigenvalues may be extracted.    |  
c        |                                           |
c        | Eigenvectors may be also computed now if  |
c        | desired.  (indicated by rvec = .true.)    | 
c        |                                           |
c        | The routine DSEUPD now called to do this  |
c        | post processing (Other modes may require  |
c        | more complicated post processing than     |
c        | mode1.)                                   |
c        |                                           |
c        %-------------------------------------------%
c           
          rvec = .true.
c
          call dseupd ( rvec, 'All', select, d, v, n, sigma, 
     &         bmat, n, which, nev, tol, resid, ncv, v, n, 
     &         iparam, ipntr, workd, workl, lworkl, ierr )
c
c         %----------------------------------------------%
c         | Eigenvalues are returned in the first column |
c         | of the two dimensional array D and the       |
c         | corresponding eigenvectors are returned in   |
c         | the first NCONV (=IPARAM(5)) columns of the  |
c         | two dimensional array V if requested.        |
c         | Otherwise, an orthogonal basis for the       |
c         | invariant subspace corresponding to the      |
c         | eigenvalues in D is returned in V.           |
c         %----------------------------------------------%
c
          if ( ierr .ne. 0) then
c
c            %------------------------------------%
c            | Error condition:                   |
c            | Check the documentation of DSEUPD. |
c            %------------------------------------%
c
             print *, ' '
             print *, ' Error with _seupd, info = ', ierr
             print *, ' Check the documentation of _seupd. '
             print *, ' '
c
          else
c
             nconv =  iparam(5)
             do 20 j=1, nconv
c
c               %---------------------------%
c               | Compute the residual norm |
c               |                           |
c               |   ||  A*x - lambda*x ||   |
c               |                           |
c               | for the NCONV accurately  |
c               | computed eigenvalues and  |
c               | eigenvectors.  (iparam(5) |
c               | indicates how many are    |
c               | accurate to the requested |
c               | tolerance)                |
c               %---------------------------%
c
                call av(n, v(1,j), ax,
     $               n_nz, element_i, element_j, element_s)
                call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))
c
 20          continue
c
c            %-----------------------------%
c            | Display computed residuals. |
c            %-----------------------------%
c
             call dmout(6, nconv, 2, d, ncv, -6,
     &            'Ritz values and relative residuals')

          end if
c
c         %-------------------------------------------%
c         | Print additional convergence information. |
c         %-------------------------------------------%
c
          if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
          else if ( info .eq. 3) then
             print *, ' ' 
             print *, ' No shifts could be applied during implicit',
     &                ' Arnoldi update, try increasing NCV.'
             print *, ' '
          end if      
c
          print *, ' '
          print *, ' END DIAGARPACK'
          print *, ' ====== '
          print *, ' '
          print *, ' Size of the matrix is                  ', n
          print *, ' The number of Ritz values requested is ', nev
          print *, ' The number of Arnoldi vectors generated',
     &             ' (NCV) is                               ', ncv
          print *, ' What portion of the spectrum:          ', which
          print *, ' The number of converged Ritz values is ', 
     &               nconv 
          print *, ' The number of Implicit Arnoldi update',
     &             ' iterations taken is                    ', iparam(3)
          print *, ' The number of OP*x is                  ', iparam(9)
          print *, ' The convergence criterion is           ', tol
          print *, ' '
c
      end if
c
      return
      end

c 
c ------------------------------------------------------------------
c     matrix vector subroutine
c
c     Computes w <--- OP*v, where OP is the n*n by n*n block 
c     tridiagonal matrix
c
      subroutine av(n, v, w, n_nz, ei, ej, es)
      integer          n, n_nz
      integer          ei(n_nz), ej(n_nz)
      real*8 es(n_nz)
      real*8 v(n), w(n)
c     
      integer k, i, j
      
      do k = 1, n
         w(k)=0
      enddo
c$$$      print *,'w', w(1:n)
c$$$      print *,'v', v(1:n)
      do k = 1, n_nz
         i = ei(k)
         j = ej(k)
         if( i .eq. j ) then
            w(i) = w(i)+es(k)*v(i)
         else
            w(i) = w(i)+es(k)*v(j)
            w(j) = w(j)+es(k)*v(i)
         endif
      enddo

      return
      end

c-------------------------------------------------------------------
      subroutine read_matrix(n_nz, ei, ej, es, fill)
      
            implicit none
            integer         n_nz, ei(*), ej(*)
            real*8            es(*)
            
            integer         m, a, b
            real*8          x
            logical fill
            integer iostat
            logical file_exists
            logical formatted
      
            inquire(FILE='matrice.sdij', EXIST=file_exists)
            if(file_exists) then
                formatted = .false.
                print *, 'reading binary matrice.sdij file'
            else
                inquire(FILE='matrice.sdijf', EXIST=file_exists)
                if(file_exists) then
                    formatted = .true.
                    print *, 'reading formatted matrice.sdijf file'
                else
                    print *, 'no matrice.sdij or matrice.sdijf file'
                    stop
                end if
            end if
      
            if(formatted) then
                open(UNIT=10,FILE='matrice.sdijf',STATUS='OLD',
     &              FORM='FORMATTED')
            else
                open(UNIT=10,FILE='matrice.sdij',STATUS='OLD',
     &              FORM='UNFORMATTED')
            end if
      
            m=0
            do while (.true.)
                if(formatted) then
                    read(10,*,iostat=iostat) a, b, x
                else
                    read(10,iostat=iostat) a, b, x
                endif
                if( iostat < 0 ) then
                    print *, m, 'lines read from matrice.sdij'
                    exit
                endif
                m = m + 1
                if(fill) then
                    ei(m) = a
                    ej(m) = b
                    es(m) = x
                end if
            end do
      
            n_nz = m
      
            close(10)
            return
            end

      function num_nonzero(si)
            ! Return the number of elements. Assumes the first si == 0 is the end
            implicit none
            integer ::     num_nonzero

            integer, intent(in) ::      si(:)

            integer          m
      
            do m=1, size(si)
                  if(si(m) .eq. 0 ) then
                        exit
                  endif
            enddo
            num_nonzero = m - 1
            
      end
      
      function matrix_dim( n_nz, ei, ej)
            implicit none
            integer          matrix_dim
            integer          n_nz, ei(*), ej(*)
            
            integer          m, max
      
            max=0
            do m=1, n_nz
               if( max .lt. ei(m) ) then
                  max = ei(m)
               endif
               if( max .lt. ej(m) ) then
                  max = ej(m)
               endif
            enddo
            matrix_dim = max
      
      end
      
      function matrix_tr( n_nz, ei, ej, es )
            implicit none
            real*8 matrix_tr
            integer          n_nz, ei(*), ej(*)
            real*8 es(*)
            
            integer          m
            real*8 tr
      
            tr=0
            do m=1, n_nz
               if( ei(m) .eq. ej(m) ) then
                  tr = tr + es(m)
               endif
            enddo
            matrix_tr = tr
      end
      
      end module diagarpack
