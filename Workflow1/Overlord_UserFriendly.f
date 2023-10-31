      module holdingCell
         double precision Ri_in, Ri_out, Rtol, tolz, slope, rise, rlinetest, diff, rmindiff, delP, elemLen,
     1   rmu, Q_max, ubar, Tau, R_cur, rhalfR, delR, timeinc, thetaval, elementnum, temptol, Q, Qslope, Qlast
        integer icounter, icounter2, ifinal, ind2, ind4, icounter3, ind_umat, ind_umat2(2), ind_surf, numQ,
     1   mastercount, ind_press, icheck, icheck_dload, ind_dload2(2), ind_pu, IFLAG, ind_ord,numelem,numel,linenum
         DIMENSION x_coord(100000), y_coord(100000), z_coord(100000), R_val(100000), P_val(100000,2), elnum(100000),
     1   elemarray(100000,4),tempelemarray(100000,4), profile(100000), props(3), Q_prof(10000,2), step_time(1)
         REAL rL, rkval
         CHARACTER(500) directory
         CHARACTER(500) tempchar
         CHARACTER(500) checkfile
         CHARACTER(500) coordfile
         CHARACTER(500) getinputfile
         CHARACTER(5) searchinp
         TYPE :: inp_read
            double precision :: startinc
            double precision :: steptime
            double precision :: smallinc
            double precision :: finalinc
         END TYPE inp_read
         TYPE(inp_read), DIMENSION(1) :: report

      end module holdingCell
C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     ! References:
C     ! https://www.researchgate.net/post/Is-anyone-experienced-with-UExternalDB-subroutines-in-Abaqus
C     ! https://abaqus-docs.mit.edu/2017/English/SIMAINPRefResources/uwavexx3.f
C     ! https://www.reddit.com/r/fea/comments/idc2x/abaqus_user_subroutines_and_opening_files/
C     ! Skeleton script by Collette Gillaspie
C     ! Code by Hannah Stroud

C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C     !***************************************************
C     !                      READ ME                      
C     !***************************************************

C     ! Enough of this code is hardcoded to be a pain in the ass. You need the number of elements in the model
C     ! (or possibly just some number greater than that number) for UTRACLOAD and the total STEP time (if different 
C     ! from 1) for UTRACLOAD. It is possible to have fortran open the input file to get this info. See locate subroutine.
C     ! You also need the filepath for coord and check files.

C     ! The way it is right now, this code (UTRACLOAD implementing H-P pipeflow) CANNOT
C     ! be run on multiple processors.

C     ! It also uses slightly dumb tolerances (ztol and one applied to thetaval)
C     ! to get COORD information and apply the desired pressure load. If things are not working, check if atan
C     ! is giving a singularity in theta (i.e. 2 theta values are counted the same even though they are 180deg
C     ! apart-- this is easily fixable if you care)... you can get rid of this issue with math or by making the number
C     ! of elements around your pipe odd.
C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

C     !***************************************************
C     ! User Subroutine to Manipulate User External Files
C     !***************************************************
      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C     ! LOP=0 indicates that the user subroutine is being called at the start of the analysis.
C     ! LOP=1 indicates that the user subroutine is being called at the start of the current analysis increment...
C     ! The user subroutine can be called rmultiple times at the beginning of an analysis increment if the...
C     ! increment fails to converge and a smaller time increment is required.
C     ! LOP=2 indicates that the user subroutine is being called at the end of the current analysis increment...
C     ! When LOP=2, all information that you need to restart the analysis should be written to external files.
C     ! LOP=3 indicates that the user subroutine is being called at the end of the analysis.
C     ! LOP=4 indicates that the user subroutine is being called at the beginning of a restart analysis. When LOP=4,...
C     ! all necessary external files should be opened and properly positioned and all information required for the...
C     ! restart should be read from the external files.
C     ! LOP=5 indicates that the user subroutine is being called at the start of a step. The KSTEP argument contains...
C     ! the current step number.
C     ! LOP=6 indicates that the user subroutine is being called at the end of a step. The KSTEP argument contains...
C     ! the current step number.
C     ! LRESTART=0 indicates that an analysis restart file is not being written for this increment.
C     ! LRESTART=1 indicates that an analysis restart file is being written for this increment.
C     ! LRESTART=2 indicates that an analysis restart file is being written for this increment and that only one...
C     ! increment is being retained per step so that the current increment overwrites the previous increment in the...
C     ! restart file (see Restarting an analysis).
C     ! TIME(1) = Value of current step time.
C     ! TIME(2) = Value of current total time.
C     ! DTIME = Time increment.
C     ! KSTEP = Current step number. When LOP=4, KSTEP gives the restart step number.
C     ! KINC = Current increment number. When LOP=4, KINC gives the restart increment number.
C     ! NOTE on LOP: Dummy variables are passed in & out of subroutines/functions by default; e.g. If an argument 
C     ! passed to a subroutine is modified in the subroutine, that modified value is returned back to the main 
C     ! program. LOP appears to be intended to inform UEXTERNALDB where in the analysis solution the program is when 
C     ! the subroutine is called. Use this information to decide what you wish to do in the subroutine, e.g. a case 
C     ! structure or if statement block.

      use holdingCell

      INCLUDE 'ABA_PARAM.INC'
      logical found


      DIMENSION TIME(2)

      write(*,*) 'We are in UEXTERNALDB before Case Selector, this is LOP number', LOP
        rL = 1.0
        rkval = 0.005

        ubar = 6.0
        rmu = 1.793e-3
        Q_max = 0.471238898038469*1e-3
        props=(/rmu,Q,ubar/)
        timeinc = 0.01
        numelem=25760
        mastercount=0
        directory ='C:\ALETube\PressureUpdate\LargeStep\DeformErode'
        directory=trim(directory)
        tempchar='\check.txt'
        checkfile='C:\ALETube\PressureUpdate\LargeStep\DeformErode\check.txt'
        tempchar='\coords.txt'
        coordfile='C:\ALETube\PressureUpdate\LargeStep\DeformErode\coords.txt'
        getinputfile='C:\ALETube\PressureUpdate\LargeStep\DeformErode\JobErode.inp'
        searchinp='*Step'
        searchinp=trim(searchinp)

        write(*,*)searchinp
        
      SELECTCASE(LOP)
C        ! LOP=0 indicates that the subroutine is being called at the start of the analysis.
         CASE(0)
            write(*,*) 'We are in LOP=0'
            write(*,*) 'The current UEXTERNALDB time is:', time
            
            open(13, file=checkfile, status='old')
            close(13, status='delete')
            open(13, file=checkfile, status='new')
            close(13)
            
            open(12, file=coordfile, status='old')
            close(12, status='delete')
            open(12, file=coordfile, status='new')
            close(12)
            
            !Get step time to develop Q profile
            open(11, file=getinputfile, action='read')
            call locate(11, searchinp, found, linenum)
            linenum=linenum+1
            if(found)then
              read(11,*)
              read(11,*)report
            endif
            close(11)
            step_time=report%steptime
            
            ! INPUT Q PROFILE for t=1,step_time-- Q_prof(:,1)=time, Q_prof(:,2)=value of Q
            Q_prof(1,1)=0.
            Q_prof(2,1)=0.5*step_time(1)
            Q_prof(3,1)=step_time(1)

            Q_prof(1,2)=0.
            Q_prof(2,2)=Q_max
            Q_prof(3,2)=Q_max
            numQ = 2 !number of intervals for which Q is defined
            Qlast=0.0
C        ! LOP=1 indicates that the subroutine is being called at the start of the current analysis 
C        ! increment. The subroutine can be called rmultiple times at the beginning of an analysis 
C        ! increment if the increment fails to converge and a smaller time increment is required.
         CASE(1)
            write(*,*) 'We are in LOP=1'
            write(*,*) 'The current UEXTERNALDB time is:', time, KINC
            
            ! Calculate current value of Q
            ! Search over Q_prof(:,1) to see where time(2) falls
            ! Q is linear between segments
            do i=1,numQ
              if(time(2).ge.Q_prof(i,1) .and. time(2).lt.Q_prof(i+1,1))then
                exit
              endif
            enddo
            Qslope=(Q_prof(i+1,2)-Q_prof(i,2))/(Q_prof(i+1,1)-Q_prof(i,1))
            Q=Qslope*time(2)+Q_prof(i,2)
            
            icounter=0
            icounter2=1
            icheck=0
            open(13, file=checkfile, status='old', action="write")
            write(13,*)'Increment',KINC
            
C        ! LOP=2 indicates that the subroutine is being called at the end of the current analysis 
C        ! increment. When LOP=2, all information that you need to restart the analysis should be 
C        ! written to external files.
         CASE(2)
            write(*,*) 'We are in LOP=2'
            write(*,*) 'The current UEXTERNALDB time is:', time, KINC
            
            !Sort the array of element coordinates (R, THETA, Z, ELNUM) in increasing order of Z
            call sort(elemarray,ifinal,4,3)
            
            !Initialize pressure profile
            do i=1,numelem
              P_val(i,:)=(/0.,0./)
            enddo
            write(*,*)'before tolz'
            tolz=1.
            temptol=1.
            do j=2,ifinal
              temptol=elemarray(j,3)-elemarray(j-1,3)
              if(temptol.LT.tolz)then
                tolz=temptol
              endif
            enddo
            write(*,*)'tolz',tolz
            
            !Calculate pressure profile
            call getPressure(elemarray,ifinal,rmu,Q,P_val)
            open(13, file=checkfile, status='old', action="write")
            do i=1,ifinal
              write(13,*)P_val(i,:)
            enddo
            
C            
            open(12, file=coordfile, status='old', action="write")
            write(12,*)'Increment',KINC
            do i=1,ifinal
              write(12,*)elemarray(i,:)
            enddo
            do j=1,numelem
              tempelemarray(j,:)=(/0.,0.,0.,0./)
              elemarray(j,:)=(/0.,0.,0.,0./)
            enddo
            write(*,*)'icounter',icounter !calls to UTRAC per increment
            write(*,*)'icounter2',icounter2 !calls to UTRAC per iteration
            Qlast=Q
C        ! LOP=3 indicates that the user subroutine is being called at the end of the analysis.
         CASE(3)
            write(*,*) 'We are in LOP=3'
            write(*,*) 'The current UEXTERNALDB time is:', time
C        ! LOP=4 indicates that the user subroutine is being called at the beginning of a restart analysis. When LOP=4,...
C        ! all necessary external files should be opened and properly positioned and all information required for the...
C        ! restart should be read from the external files.
         CASE(4)
            write(*,*) 'We are in LOP=4'
            write(*,*) 'The current UEXTERNALDB time is:', time
C        ! LOP=5 indicates that the user subroutine is being called at the start of a step. The KSTEP 
C        ! argument contains the current step number.
         CASE(5)
            write(*,*) 'We are in LOP=5'
            write(*,*) 'The current UEXTERNALDB time is:', time
C        ! LOP=6 indicates that the user subroutine is being called at the end of a step. The KSTEP argument contains...
C        ! the current step number.
         CASE(6)
            write(*,*) 'We are in LOP=6'
            write(*,*) 'The current UEXTERNALDB time is:', time
      END SELECT\
      
      write(*,*) 'We are in UEXTERNALDB after Case Selector'
C     !************************
C     ! FORMATS
C     !************************
  530 FORMAT(E12.4,E12.4,E12.4,E12.4)
      return
      end subroutine UEXTERNALDB
      

C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
C     !**********************
C     ! User Subroutine UMESHMOTION
C     !**********************
      subroutine umeshmotion(uref,ulocal,node,nndof,lnodetype,alocal,
     *     ndim,time,dtime,pnewdt,kstep,kinc,
     *     kmeshsweep,jmatyp,jgvblock,lsmooth)
      use holdingCell

      INCLUDE 'ABA_PARAM.INC'

C
      CHARACTER*80 PARTNAME
      DIMENSION ARRAY(1000)
      DIMENSION ULOCAL(*)
      DIMENSION JGVBLOCK(*),JMATYP(*)
      DIMENSION ALOCAL(NDIM,*)
      PARAMETER (NELEMMAX=1000000)
      DIMENSION JELEMLIST(NELEMMAX),JELEMTYPE(NELEMMAX)
      PARAMETER (ALamda1=0.02626D0,ALamda2=1.0e-12)
      PARAMETER (TRCON = 0.001)
      PARAMETER (ELINC = 0.1D0)
      DIMENSION TIME(2)
      DIMENSION CFLOCAL(3),CFGLOBAL(3)
      DIMENSION ASUB(3,3), AMAT(3,3), AMATT(3,3)
      REAL*8    CF1, CF2, CF3, CORR, CORR2, R_UMESH

      !write(*,*) 'We are in UMESHMOTION'

      LOCNUM = 0
      JRCD = 0
      LTRN = 0
      PARTNAME = ' '
      PEEQ = 0.0D0
      FPEEQ= 0.0D0
      FPOROSITY=0.00D0
      FRVF = 0.00D0
      CHARLENGTH = 0.005
      JTYP = 1
      DO I=1,3
        DO J=1,3
          ASUB(I,J)=0.0D0
          AMAT(I,J)=0.0D0
          AMATT(I,J)=0.0D0
        ENDDO
      ENDDO

      CALL GETPARTINFO(NODE,0,PARTNAME,LOCNUM,JRCD)
      NELEMS = NELEMMAX
      CALL GETNODETOELEMCONN(NODE,NELEMS,JELEMLIST,JELEMTYPE,JRCD,
     $     JGVBLOCK)
      CALL GETVRN(NODE,'COORD',ARRAY,JRCD,JGVBLOCK,LTRN)
      ! CF1 = ARRAY(1)
      ! CF2 = ARRAY(2)
      ! CF3 = ARRAY(3)
      R_UMESH=SQRT(ARRAY(1)**2+ARRAY(2)**2)
      ! CFGLOBAL(1)=CF1
      ! CFGLOBAL(2)=CF2
      ! CFGLOBAL(3)=CF3
      
      ! DO I=1,NDIM
        ! DO J=1,NDIM
          ! CFLOCAL(I)=CFLOCAL(I)+CFGLOBAL(J)*ALOCAL(J,I)
        ! END DO
      ! END DO

      !TAU=SQRT(CFLOCAL(3)**2+CFLOCAL(2)**2)
      Tau=4.*rmu*Qlast/(3.14159*R_UMESH**3)
      UREF=UREF*SQRT(2.)
      DO I=1,3
        DO J=1,3
          ASUB(I,J)=ALOCAL(I,J)
        ENDDO
      ENDDO
      
      ASUB(3,3)=0.
      CORR2=0.
      DO I=1,3
        DO J=1,3
          CORR2=CORR2+ASUB(I,J)**2.
        ENDDO
      ENDDO

          
C     CORR=0.5*SQRT(CORR2)
      CORR=sqrt(ALOCAL(3,1)**2.+ALOCAL(3,2)**2.)

      
      IF(LNODETYPE.EQ.3)THEN
        SURFV=UREF*1./CORR
      ELSEIF(LNODETYPE.EQ.4)THEN
        SURFV = UREF
      ELSEIF(LNODETYPE.EQ.5)THEN
        SURFV = UREF/SQRT(2.)
      ENDIF

        SURFV = SURFV*TAU!+ALamda2

          ULOCAL(NDIM) = ULOCAL(NDIM)-SURFV

C       SHARE OF ELEMENT CONSUMED IN NEXT INCREMENT
        PELEM = DTIME * SURFV / CHARLENGTH
        PNEWDT1 = 0.0D0
        IF (PELEM.GT.ELINC.AND.SURFV.GT.0.0D0) THEN
          PNEWDT1 = CHARLENGTH * ELINC / SURFV
          IF (PNEWDT1.LT.DTIME) THEN
            WRITE (7,*) 'CHANGING TIME INCREMENT FROM ',DTIME
            PNEWDT = PNEWDT1
            WRITE (7,*) 'TO ',PNEWDT
            WRITE (7,*) 'BASED ON NODE ',LOCNUM
            PNEWDT = PNEWDT / DTIME
          END IF
        END IF
        LSMOOTH = 1

      return
      end



         
C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     user subroutine for traction load
c
      subroutine utracload(alpha,t_user,kstep,kinc,times,noel,npt,
     1     coords,dircos,jltyp,surface_name)
      use holdingCell

      INCLUDE 'ABA_PARAM.INC'
c
      dimension t_user(3),times(*),coords(*), dircos(3,*)
      real Pmatch, Pmatch1, Zmatch
      character*80 surface_name
c
      !write(*,*) 'We are in UTRACLOAD',ifinal
      ! if we go through the elements again (i.e. since abaqus calls elements in numeric
      ! order, noel that is less than the previous noel indicates a new call) reset the 
      ! counter for elemarray such that elements are only stored once
      if(noel.lt.numel)then
        !if(kinc.le.2)then
          ifinal=icounter2-1
        !endif
        icounter2=1
      endif
      numel=noel
      elementnum=dble(noel)
      PI=2.*ASIN(1.D0)
      THETA=ATAN(COORDS(2)/COORDS(1))
      R=SQRT(COORDS(1)**2+COORDS(2)**2)
      Z=COORDS(3)

      tempelemarray(:,:)=elemarray(:,:)
C   If this is the first call to UTRAC in the increment, store theta
      if(icounter.eq.0)then
        thetaval=THETA
      endif

C   Track Z and R for one theta value to get R profile sorted by Z
      if(icheck.lt.1)then
      if(THETA.LE.thetaval+1e-4 .AND. THETA.GE.thetaval-1e-4)then
        
        elemarray(icounter2,:)=(/R,THETA,Z,elementnum/)
        
        icounter2=1+icounter2
      endif
      endif
      ! if(kstep .eq. 1) then
         ! alpha=1.0*times(1)
      ! else if(kstep .eq. 2) then
         ! alpha=1.0d0
      ! endif
      Pmatch=0.
      Zmatch=1
      !Find Pmatch based on Z coord
      do i=1,numelem
        if(P_val(i,2).le.Z+0.99*tolz .and. P_val(i,2).ge.Z-0.99*tolz)then
          Pmatch=P_val(i,1)
          Zmatch=P_val(i,2)
        endif
      enddo
      Pmatch1=35001.*0.5*Z!*times(1)
      !Apply pressure from P_val
      if(kinc.eq.1)then
        alpha=0.
      else
        alpha=Pmatch
      endif
c
      t_user(1)=1.0
      t_user(2)=0.0
      t_user(3)=0.0
      icounter=1+icounter
      if(icounter.gt.numelem*8)then
        icheck=1
      endif

         
      RETURN
      END
C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     !************************
C     ! Utility Routines
C     !************************
        SUBROUTINE sort(arr,len1,len2,indsort)

        DIMENSION arr(100000,len2), temp(1,len2)
        INTEGER iadjust, len1, len2, indsort, i, k


c        write(*,*)'in sort subroutine'
        iadjust = len1 - 3
c        write(*,*)iadjust, len1, len2, indsort
c        write(*,*)arr(1,:)
        temp(1,:)=(/0.,0.,0.,0./)
c        write(*,*)temp(1,:)


        DO 20, I = 1,len1
           DO 10, k = len1,len1 - iadjust -1, -1
               IF ( arr(k,indsort) .LT. arr(k - 1,indsort) ) THEN
C                 SWAP if A(k,indsort) < A(k-1,indsort)
                  do 30 l=1,len2
                  temp(1,l) = arr(k,l)
                  arr(k,l) = arr(k-1,l)
                  arr(k-1,l) = temp(1,l)

   30             continue
               END IF
   10     CONTINUE

   20  CONTINUE
      RETURN
      END

C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        SUBROUTINE getPressure(arr,len1,rmu,Q,press)

        DIMENSION arr(100000,4), props(3), press(100000,2)
        REAL*8 const, rmu, Q, ubar, delP, delZ

        write(*,*)'in pressure subroutine'

        
        const=8.*rmu*Q/3.14159
        
        do 10 i=2,len1
          delZ=arr(i,3)-arr(i-1,3)
          delP=const*delZ/arr(i,1)**4
          press(i,1)=press(i-1,1)+delP
          press(i,2)=arr(i,3)
   10   continue

      RETURN
      END
C     ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine locate(fileid, keyword, have_data, line_num)
          
          integer,intent(in)          :: fileid               ! File unit number
          character(5) :: keyword              ! Input keyword, 5 characters here
          logical,intent(out)         :: have_data            ! Flag: input found
          integer                     :: line_num             ! Line number of keyword

      
          character*(100)             :: linestring           ! First 100 chars
          integer                     :: keyword_length       ! Length of keyword
          integer                     :: io                   ! File status flag
          keyword=trim(keyword)
          keyword_length = len(keyword)
          rewind(fileid)
          line_num=1
          k=1
          ! Loop until end of file or keyword found
          do
              ! Read first 100 characters of line
              read (fileid,'(a)',iostat=io) linestring
                
              ! If end of file is reached, exit
              if (io.ne.0) then 
                  have_data = .false.
                  exit 
              end if
              
              ! If the first characters match keyword, exit
              if (linestring(1:keyword_length).eq.keyword) then
                  have_data = .true.
                  line_num = k
                  exit
              endif
              k=k+1
          end do
      
      RETURN
      END