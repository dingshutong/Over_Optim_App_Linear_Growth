cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  ----- Using C. Fernandez, E. Ley & M.J.F. Steel (2001) F77 code -------
c  ----- To check the flip-flopping behavior of variable selection -------
c  ----- Modified by Shutong Ding ----------------------------------------
c  ----- 2015-05-22 ------------------------------------------------------    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc    
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                            c
c   VERSION: 06.08.27                                                        c
c                                                                            c
c   Modified version of (January 22, 2001) f77 code used in                  c
c                                                                            c
c   C. Fernandez, E. Ley & M.J.F. Steel (2001)                               c
c     "Model Uncertainty in Cross-country Growth Regressions"                c
c     Journal of Applied Econometrics                                        c
c                                                                            c
c   Updates:                                                                 c
c                                                                            c
c   [1] This version works for # regressors =< 104.                          c
c   [2] This version computes jointness measures                             c
c       (see Ley and Steel (2006) "Jointness in Bayesian Model               c
c       Selection---with Application to Growth Regression")                  c
c                                                                            c
c   INPUT:                                                                   c
c         fls.par (first line contains namenamename)                         c
c         fdatname.dat  (grk41t72.dat or grk67t88.dat)                       c
c                                                                            c
c   OUTPUT:                                                                  c
c      namenamename.out                                                      c
c      namenamenamebj.dat  (contains bi-jointness measures)                  c
c      namenamenametj.dat  (contains tri-jointness measures)                 c
c                                                                            c
c   NOTE COMPILING OPTIONS:                                                  c
c     1. must link with unix libraries for wr_date() to work,                c
c        if not available simply remove it or change its body to do          c
c        nothing (ie, just 'return')                                         c
c     2. must compile using static storage                                   c
c     3. optimized code will run significantly faster                        c
c                                                                            c
c                                                                            c
c   fls.par file                                                             c
c   ============                                                             c
c                                                                            c
c   outfilename0            out namefile, change at wish, char12             c
c   grk41t72.dat            dat namefile, grk67t88.dat is the other  file    c
c   -1665432                random number seed, enter any integer            c
c   9                       integer (1-9) specifiying prior                  c
c   100000                  warmup draws---set to small number for testing   c
c   200000                  chain draws                                      c
c   F                       standardise xs ??                                c
c   F                       lps?                                             c
c   5                       lpsloop? (must be < 26)                          c
c   0.85d0                  real taking care of sample split                 c
c   T                       wrpost?                                          c
c   T                       dogm?---do g&m stuff? convergence                c
c   T                       dojoint?---do jointness stuff?                   c
c                                                                            c
c   Data Files                                                               c
c   ==========                                                               c
c                                                                            c
c   The structure of the *.dat file is 1st line with kreg, 2nd line          c
c   with nobs, then kreg lines w/ varnames, and then nobs lines each         c
c   with (kreg+1) data values: (y,Z):                                        c
c   ................................................................         c
c   kreg       Note: kreg must be =< 104                                     c
c   nobs       No limit but if nobs > 88 change maxn=88 throughout           c
c   <# kreg lines with the names of the kreg vars>                           c
c   <# nobs lines with kreg+1 format-free variables: (y,Z) >                 c
c                                                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
      program bma
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8 (a-h,o-z)
      
 !     include 'mkl_vsl.f77' ! to use the MKL library VSL

      integer fout,foutv,foutm,foutbj,fouttj
      parameter (fout=16,foutv = 9, foutm = 11, foutbj=22,fouttj=23)
      parameter (maxk=104,maxn=88,maxm=200000,maxnf=50) 
      logical fail,wrpost,dojoint,dogm,lpsloop
      logical standard,fix_start
      integer kj(maxm),mmodel(maxk),idx(maxm),models(maxk)
      character*12 regname(maxk),fdatname
  
  
      real*8 z(maxn,maxk),y(maxn),ztz(maxk,maxk),yf(maxnf),
     &   zf(maxk,maxnf),midx(maxm,2),
     &   bayesf(maxm),freq(maxm),gj(maxm),dstar(maxm),
     &   bstar(maxk,maxm),xdata(maxn,maxk),ydata(maxn)
     
      
c
c  if you want to do more than 25 sample splits change the 25
c  below and also in the setup() statement: ilps = min(ilps,25)   
c      

       real*8 avelps(25),bestlps(25),fulllps(25),nulllps(25)

c  some new definition in bootstrap
      integer bootrep, maxboot, bi,brng, errcode,method,prev_acc
      integer max_top, top_num,ir(maxk),kjj,acc_count,fail_count
      parameter(maxboot = 10000, max_top = 100)
      integer boot_inds(maxn), rands(maxboot)
      real*8  top_midx(max_top,2)
      real*8  xdata_org(maxn,maxk),ydata_org(maxn)
      real*8  beta(maxk),percvisits, top_mod_post(max_top)
  !    TYPE (VSL_STREAM_STATE) :: stream
 !     common method, stream
c
cccccccccccccccccccccccccc
ccc   setup            ccc
cccccccccccccccccccccccccc
c
      acc_count = 0
      fail_count = 0
       call setup(               
     &      iprior,idum,initrep,mnumrep,bootrep,top_num,
     &      fdatname,standard,
     &      lpsloop,ilps,split,
     &      wrpost,dogm,dojoint,
     &      fout,foutv,foutm,foutbj,fouttj,maxboot,max_top,fix_start)
           
       call readdata(fout,fdatname,
     &       lpsloop,ntot,split,nobs,nf,standard,
     &       kreg,regname,xdata_org,ydata_org)
       

c     Generating random indexes for all bootstrap  
      ! generate the random seed first
!     ***** Initializing *****
 !     method = VSL_RNG_METHOD_UNIFORM_STD
 !     brng=VSL_BRNG_MT19937
 !     errcode=vslnewstream( stream, brng,  idum )
      ! if (errcode .ne. 0) then 
      !   write(fout,*) '!!!!! ERROR: bootstrap seed genrates fail !!!!'
      !write(*,*)    '!!!!! ERROR: bootstrap seed genrates fail !!!!'  
       !  stop
      !end if    
!      do 5102 bi = 1,bootrep
 !        call gen_srs_with_replace (nobs, nobs, 
 !    &   boot_inds(1:nobs,bi))        
 !       call ran_srs(idum,nobs,nobs,boot_inds(1:nobs,bi))
!5102   end do

  !     ***** Deinitialize *****
  !    errcode=vsldeletestream( stream ) 
      
   !   if (errcode .ne. 0) then 
   !      write(fout,*) '!!!!! ERROR: bootstrap seed genrates fail !!!!'
   !   write(*,*)    '!!!!! ERROR: bootstrap seed genrates fail !!!!'  
   !      stop
   !   end if    
    
 !      output the random seed of bootstrap
          open(unit=3, file='bootstrap_seed_row.out',
     &    ACTION="write", STATUS="replace")
          open(unit=13, file='bootstrap_seed_failed.out',
     &    ACTION="write", STATUS="replace")
 !         do bi = 1,bootrep
        !     write(3,'(100i3)') (boot_inds(i,bi), i=1,nobs )        
!          end do
!          close(3)
 
    ! run the MC3 on the original data once to find the optimial models
        call wr_date(fout,.true.,.true.)
        call wr_date(foutm,.true.,.true.)
        write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        write(fout,*) 'run the MC3 on the original data once' 
        write(fout,*) 'in order to find the optimial models' 
        write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
               
      call splitdata(fout,regname,idum,ntot,nobs,nf,kreg,
     &      xdata_org,z,ztz,ydata_org,y,avey,ssqyn,lpsloop,yf,zf,fail)
      
c      write(*,*) 'Splitdata:',fail
      
      call ols(fout,regname,nobs,kreg,z,ztz,y,avey,fail)
      
c      write(*,*) 'OLS:',fail
      
    !  call wr_date(fout,.true.,.true.)
      
      if (fail) then
         write(fout,*) '!!! setup block fail !!!'
         write(*,*)    '!!! setup block fail !!!' 
         close(fout)
         stop  
      else
         write(*,*)    ' ... setup block done ... '         
      endif
c
cccccccccccccccccccccccccc
ccc   run chain        ccc
cccccccccccccccccccccccccc
c
      call runchain(fix_start,iprior,idum,initrep,mnumrep,nobs,kreg,
     &   ztz,z,y,ssqyn,
     &   freq,bayesf,gj,dstar,bstar,midx,kj,
     &   imax,nvout,fout,fjout,fjsum,fail)

      write(*,*) '... chain done! ',imax
      if (fail) then
         write(fout,*) ' !!! fail running chain!!'
         stop
      endif
      if (nvout.gt.0) then
         write(fout,*) ' !!! ',nvout,' visits out!!!'
         write(fout,*) ' !!! mass',fjout,' visits out!!!'
      endif
 !     call wr_date(fout,.true.,.true.)
c
ccccccccccccccccccccccccccc
ccc  write chain info   ccc
ccccccccccccccccccccccccccc
c
      write(*,*) '... calling wrchainfo() ...' 

      
      call wrchainfo(regname,kreg,bayesf,midx,bstar,gj,kj,
     &      imax,nvout,ibm,icut,idx,fout, .false.,beta)
      
      ! idx save the desending order of models
      write(foutm,*) '#################################################'
      write(foutm,*) '# Top models have been picked by original data  #'
      write(foutm,*) '#################################################'
      call barra('-',55,foutm)
      do i = 1,top_num    
                 ! print the top models information in foutm
         call g2model2(midx(idx(i),1),midx(idx(i),2),
     &         kreg,models,kjj,ir)      
            percvisits = bayesf(idx(i))*100.0d0
            write(foutm,1467) i, percvisits,(ir(k),k=1,kjj)
 1467       format(i3,t5,f6.2,'%',x,'|',100(:i3))         
         ! return the indexes of top models
           top_midx(i,1:2) = midx(idx(i), 1:2)
      end do    
      call barra('-',55,foutm)
          
        call wr_date(fout,.true.,.true.)
        call wr_date(foutv,.true.,.true.)
        call wr_date(foutm,.true.,.true.)
        
        write(fout,*) '==============================================='
        write(fout,*) 'Boostrap starts' 
        write(fout,*) '==============================================='
        write(foutv,*) '==============================================='
        write(foutv,*) 'Boostrap starts' 
        write(foutv,*) '==============================================='
        write(foutm,*) '==============================================='
        write(foutm,*) 'Boostrap starts' 
        write(foutm,*) '==============================================='
                    
      
       write(*,*) '==============================================='
       write(*,*) 'Boostrap starts' 
       write(*,*) '==============================================='     
       
      call barra('-',100,foutv) 
      write (foutv,'(100a20)') 'Bootstrap.rep', (regname(i), i=1,kreg)  
      call barra('-',100,foutv) 
      
      call barra('-',100,foutm) 
      write (foutm,*) 'Bootstrap.rep\Top models from 1 to ', top_num 
      call barra('-',100,foutm) 
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
   ! Boostrap starts from here, repeating the entire process      
      do 5189 bi = 1,maxboot      
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
        if (acc_count .eq. bootrep) then
            exit 
        end if 
        prev_acc = acc_count
        ! generate srs sampler
1881    call ran_srs(idum,nobs,nobs,boot_inds(1:nobs)) 
          
        write(fout,*) ' '
        write(fout,*) '------------------------------'
        write(fout,*) 'For the data replicate: ', bi
        write(fout,*) '------------------------------'
        write(*,*)    ' '
        write(*,*) '------------------------------'
        write(*,*) 'For the data replicate: ', bi
        write(*,*) '------------------------------'
        xdata = xdata_org
        xdata(1:nobs,1:maxk) = xdata_org(boot_inds(1:nobs), 1:maxk)
        ydata = ydata_org
        ydata(1:nobs) = ydata_org(boot_inds(1:nobs))  
        
          
      call splitdata(fout,regname,idum,ntot,nobs,nf,kreg,
     &      xdata,z,ztz,ydata,y,avey,ssqyn,lpsloop,yf,zf,fail)
      
c      write(*,*) 'Splitdata:',fail
      
      call ols(fout,regname,nobs,kreg,z,ztz,y,avey,fail)
      
c      write(*,*) 'OLS:',fail
      
      call wr_date(fout,.true.,.true.)
      
      if (fail) then
         write(fout,*) '!!! setup block fail !!!'
         write(*,*)    '!!! setup block fail !!!' 
         fail_count = fail_count + 1
         write(13,'(100i3)') (boot_inds(i), i=1,nobs )    
         go to 1881
         !close(fout)
         !stop
      else
         write(fout,*) '... setup block done ...'
         write(*,*)    ' ... setup block done ... '         
         acc_count = acc_count + 1   
      endif
c
cccccccccccccccccccccccccc
ccc   run chain        ccc
cccccccccccccccccccccccccc
c
      call runchain(fix_start,iprior,idum,initrep,mnumrep,nobs,kreg,
     &   ztz,z,y,ssqyn,
     &   freq,bayesf,gj,dstar,bstar,midx,kj,
     &   imax,nvout,fout,fjout,fjsum,fail)

      write(*,*) '... chain done! ',imax
      if (fail) then
         write(fout,*) ' !!! fail running chain!!'
         acc_count = acc_count - 1
         fail_count = fail_count + 1
         write(13,'(100i3)') (boot_inds(i), i=1,nobs )   
         go to 1881
        ! stop
      endif
      if (nvout.gt.0) then
         write(fout,*) ' !!! ',nvout,' visits out!!!'
         write(fout,*) ' !!! mass',fjout,' visits out!!!'
      endif
      
      if (acc_count .ne. prev_acc) then 
           write(3,'(100i3)') (boot_inds(i), i=1,nobs )    
      end if 
      
      call wr_date(fout,.true.,.true.)
      
c
ccccccccccccccccccccccccccc
ccc  write chain info   ccc
ccccccccccccccccccccccccccc
c
      write(*,*) '... calling wrchainfo() ...' 

      
      call wrchainfo(regname,kreg,bayesf,midx,bstar,gj,kj,
     &      imax,nvout,ibm,icut,idx,fout, .false.,beta)
      
      ! add line for variable inclusion
      write(foutv,5193) acc_count,  (beta(i), i=1,kreg)  
      ! add line for top model 
      call check_top_mod_post(bayesf,midx,imax,top_num,top_midx,
     & top_mod_post)
      write(foutm,5193) acc_count, (top_mod_post(i), i=1,top_num) 
      
5193  format(i3,t10,100f12.5)     
      call wr_date(fout,.true.,.true.)
      
      
      
      
c
ccccccccccccccccccccccccccc
ccc     jointness       ccc
ccccccccccccccccccccccccccc
c
            ! this part is not using in our study, sdg
      if (dojoint) then
         call jointness(kreg,bayesf,midx,
     &         imax,fout,foutbj,fouttj,fail)
         call wr_date(fout,.true.,.true.)
      endif
c
ccccccccccccccccccccccccccc
ccc     writepost       ccc
ccccccccccccccccccccccccccc
c
      if (wrpost) then
         write(*,*) '... calling wrpostinfo() ...'      
         call wrpostinfo(regname,nobs,kreg,ztz,z,
     &         bayesf,gj,dstar,bstar,midx,imax,fout,fail)
         if (fail) write(*,*) '>>> fail in wrchainfo().1'
         call wr_date(fout,.true.,.true.)
      endif
c
ccccccccccccccccccccccccccccccc
c  within-smaple prediction   c
ccccccccccccccccccccccccccccccc
c   
c 
      
      ! this part is not using in our study, sdg
      if (lpsloop) then
         write(*,*) '... looping LPS ...'   
         write(fout,*)
         call barra('L',75,fout)
         write(fout,*) 'within-of-sample prediction; nf =',nf
         write(fout,*)         
         call lps(iprior,imax,ibm,nobs,kreg,
     &         bayesf,gj,dstar,bstar,midx,
     &         avey,ssqyn,ztz,z,y,nf,yf,zf,fout,
     &         avelps(1),bestlps(1),fulllps(1),nulllps(1),
     &         fail)
         
         if (fail) write(fout,*) '>>> fail in lps()'
         call wr_date(fout,.true.,.true.)
         do 55 i=2,ilps
            call barra("L",75,fout)
            write(fout,*)
            write(fout,*) i
            write(*,*) '@@@', i,' of ',ilps
            write(fout,*)       
            
            call splitdata(fout,regname,idum,ntot,nobs,nf,kreg,
     &            xdata,z,ztz,ydata,y,avey,ssqyn,lpsloop,yf,zf,fail)
            
            call runchain(fix_start,iprior,idum,initrep,mnumrep,nobs,
     &       kreg,
     &            ztz,z,y,ssqyn,
     &            freq,bayesf,gj,dstar,bstar,midx,kj,
     &            imax,nvout,fout,fjout,fjsum,fail)
            if (fail) write(fout,*) '>>> fail runchain in inner loop'
            call wr_date(fout,.true.,.true.)
            
            call lps(iprior,imax,ibm,nobs,kreg,
     &            bayesf,gj,dstar,bstar,midx,
     &            avey,ssqyn,ztz,z,y,nf,yf,zf,fout,
     &            avelps(i),bestlps(i),fulllps(i),nulllps(i),
     &            fail)
            if (fail) write(fout,*) '>>> fail lps in inner loop'
            call wr_date(fout,.true.,.true.)
            write(fout,*)
            write(*,*)'***********************************************'
            write(*,*)'     BMA         Best        Full        Null'
            write(*,*)'-----------------------------------------------'
            do 54 ii=1,i           
               write(*,'(4f12.4)') avelps(i),bestlps(i),
     &               fulllps(i),nulllps(i)
   54       enddo
            write(*,*)'***********************************************'            
   55    enddo
         write(fout,*)
         call barra('L',75,fout)
         write(fout,*)'     BMA         Best        Full        Null'
         call barra('-',75,fout)
         do 56 i=1,ilps
            write(fout,'(4f12.4)') avelps(i),bestlps(i),
     &            fulllps(i),nulllps(i)
   56    enddo
         call barra('L',75,fout)
         write(fout,*)
         endif
c 
c
ccccccccccccccccccccccccccc
ccc    convergence      ccc
ccccccccccccccccccccccccccc
c  
      ! this part is not using in our study, sdg
      ! keep it simple   
      if (dogm) then
         call convergence(regname,
     &         iprior,idum,initrep,mnumrep,nobs,kreg,
     &         ztz,z,y,ssqyn,
     &         freq,bayesf,gj,dstar,bstar,midx,kj,
     &         imax,nvout,icut,fout,fail,fjout,mmodel,fjsum,fix_start)
      endif
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c          
     

 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc   
5189  end do
      ! Boostrap ends from here, repeating the entire process finished         
 !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
      write(foutv,*) '_________________________________________________'
      write(fout,*)  'Total number of failed bootstrap replicates :', 
     & fail_count   
        call wr_date(fout,.true.,.true.)
        call wr_date(foutv,.true.,.true.)
        call wr_date(foutm,.true.,.true.)
      
      close(foutm)
      close(foutv)
      close(3)
      close(13)
      stop
      end
      
      
      
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccc                                                                 ccc
cccc  subroutines                                                    ccc
cccc                                                                 ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine setup                                                ccc
ccc                                                                  ccc
ccc  modified:                                                       ccc
ccc    06.07.25 readdata is an independent routine                   ccc
ccc                                                                  ccc 
ccc    06.07.19 misc minutia                                         ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine setup(
     &      iprior,idum,initrep,mnumrep,bootrep,top_num,
     &      fdatname,standard,
     &      lpsloop,ilps,split,
     &      wrpost,dogm,dojoint,
     &      fout,foutv,foutm,foutbj,fouttj,maxboot,max_top,fix_start)   
      !sdg
c
      implicit real*8(a-h,o-z)
      
      integer fout,foutv,foutm,fpar,foutbj,fouttj,maxboot,top_num
      parameter(fpar=19)    
      logical wrpost,lpsloop,dogm,dojoint,standard,fix_start
      character*10 nombre
      character*12 fdatname
      integer initrep,mnumrep,bootrep,max_top    !sdg
c
c  read paremeters
c
      open(unit=fpar,file='fls.par')     

      read(fpar,'(a10)') nombre
      !read(fpar,'(a10)') var_incl_boot
      ! output for variable selction 
      open(unit=foutv, file=trim(adjustl(nombre))//'_var_incl.out')
      ! output for top model selection 
      open(unit=foutm, file=trim(adjustl(nombre))//'_mod_post.out')
      open(unit=fout, file=trim(adjustl(nombre))//'.out')
      
      read(fpar,*) fdatname
         
      call wr_date(fout,.true.,.true.)
      write(fout,*)
      write(fout,'(" This file is ",a15)') trim(adjustl(nombre))//'.out' 
      write(fout,'(" The variable out file is   ",a30)') 
     & nombre//'_var_incl.out' 
      write(fout,'(" The top model out file is   ",a30)') 
     & nombre//'_mod_post.out' 
      
      write(foutv,'(" This file is ",a20)') trim(adjustl(nombre))
     & //'_var_incl.out'     
      write(*,   '(" The all records file  is  ",a15)') 
     & trim(adjustl(nombre))//'.out'
      write(*,   '(" The variable out file is  ",a30)') 
     & trim(adjustl(nombre))//'_var_incl.out'
      write(*,   '(" The top model out file is  ",a30)') 
     & trim(adjustl(nombre))//'_mod_post.out'
      
      write(foutv,*)
      write(foutv,'("The all records file is ",a15)') 
     & trim(adjustl(nombre))//'.out'
      write(foutm,*)
      write(foutm,'("The all records file is ",a15)') 
     & trim(adjustl(nombre))//'.out' 
      write(fout,'(" Data file is ",a15)') fdatname
      write(foutv,'(" Data file is ",a15)') fdatname
      write(foutm,'(" Data file is ",a15)') fdatname
      write(*,   '(" Data file is ",a15)') fdatname
      
      ! Header for output of top models

      
      
      write(fout,*)
      write(foutv,*)
      write(foutm,*)
      call barra('-',55,fout)
            call barra('-',55,foutv)
                  call barra('-',55,foutm)
      read(fpar,*) idum
      if (idum.gt.0) idum=-idum
      write(fout,*) '.. random seed ......................',idum
      write(foutv,*) '.. random seed ......................',idum
      write(foutm,*) '.. random seed ......................',idum
      write(*,*)    '.. random seed ......................',idum
      read(fpar,*) iprior
      if ((iprior.lt.1).or.(iprior.gt.9)) then
      write(*,*)    'error: iprior must be between 1 and 9!!!'
      write(fout,*) 'error: iprior must be between 1 and 9!!!'
      write(foutv,*) 'error: iprior must be between 1 and 9!!!'
      close(fout)
      stop
      endif
      write(fout,*) '.. prior ............................',iprior
      write(foutv,*) '.. prior ............................',iprior
      write(foutm,*) '.. prior ............................',iprior      
      write(*,*)    '.. prior ............................',iprior
      read(fpar,*) initrep
      write(fout,*) '.. burn-in draws ....................',initrep
      write(foutv,*) '.. burn-in draws ....................',initrep
      write(foutm,*) '.. burn-in draws ....................',initrep
      write(*,*)    '.. burn-in draws ....................',initrep
      read(fpar,*) mnumrep
      write(fout,*) '.. mc3 draws.........................',mnumrep
      write(foutv,*) '.. mc3 draws ............................',mnumrep
      write(foutm,*) '.. mc3 draws ............................',mnumrep
      write(*,*)    '.. mc3 draws.........................',mnumrep
      read(fpar,*) bootrep 
      
      if (bootrep .gt. maxboot) then
         write(fout,*) 'Bootstrap replicates is over the limits 10000!!'
        write(foutv,*) 'Bootstrap replicates is over the limits 10000!!'
        write(foutm,*) 'Bootstrap replicates is over the limits 10000!!'
         write(*,*)    'Bootstrap replicates is over the limits 10000!!'
         close(fout)
         close(foutv)
         close(foutm)
         stop
      end if
      
      read(fpar,*) top_num  
      if (top_num .gt. max_top) then
         write(fout,*) 'Top number of models is over the limits 100!!'
        write(foutm,*) 'Top number of models is over the limits 100!!'
         write(*,*)    'Top number of models is over the limits 100!!'
         close(fout)
         close(foutv)
         close(foutm)
         stop
      end if
      
      write(fout,*) '.. bootstrap data replicates ........',bootrep
      write(foutv,*) '.. bootstrap data replicates ........',bootrep
      write(foutm,*) '.. bootstrap data replicates ........',bootrep
      write(*,*)    '.. bootstrap data replicates ........',bootrep
      
      write(fout,*) '.. The number of top models to present ..',top_num
      write(foutv,*) '.. The number of top models to present ..',top_num
      write(foutm,*) '.. The number of top models to present ..',top_num
      write(*,*) '.. The number of top models to present ..',top_num
      
      read(fpar,18) standard
      if (standard) then 
        write(fout,*) '.. Xs will be standardised'
        write(foutv,*) '.. Xs will be standardised'
        write(*,*)    '.. Xs will be standardised'
      endif
      read(fpar,18) lpsloop
      read(fpar,*) ilps
      ilps = min(ilps,25)
      if (lpsloop) then
        write(fout,*) '.. Do ',ilps,' sample splits'
        write(*,*) '.. Do ',ilps,' sample splits'
       endif
      read(fpar,*) split
      read(fpar,18) wrpost     
      read(fpar,18) dogm
      if (dogm) then
      write(fout,*) '.. will do G&M.......................'
      write(*,*)    '.. will do G&M.......................'
      endif
      
       
      read(fpar,18) dojoint
      if (dojoint) then
        write(fout,*) '.. will do Jointness Stuff...........'
        write(*,*)    '.. will do Jointness Stuff...........'
        open(unit=foutbj, file=nombre//'bj.dat')
        open(unit=fouttj, file=nombre//'tj.dat')
        write(fout,17) nombre//'.bj.out'
        write(fout,17) nombre//'.tj.out'
 17     format(' ... accompanying file: ',a20)
      endif
      
      read(fpar,18) fix_start
      if (fix_start) then
        write(fout,*) '.. The starting chain is always fixed...........'
        write(*,*)    '.. The starting chain is always fixed...........'
      end if
      
      write(foutv,*) '_________________________________________________'
      
  18   format(l1)
      close(fpar)
c 
      return
      end
c
      
      
            
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   ! Write a function to general simple random sampling with 
   ! replacement         
c       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine gen_srs_with_replace (num_to_draw, max_ind, gens)
            
!      implicit real*8(a-h,o-z)
      

!      integer num_to_draw, max_ind, gens(num_to_draw)   
!      integer method,brng, seed, i
!      INTEGER*4 errcode
!      common method, stream
!      !INTEGER*4 stream(2)
!      real a,b, xr(num_to_draw), xrt
      
      
      ! use the MKL library to sample
!       a = 0d0
!       b = 1d0

!     ***** Generating *****
!      errcode=vsrnguniform( method, stream, num_to_draw, xr ,a ,b )

       
 !     do 5181 i = 1,num_to_draw
 !        xrt = xr(i)*real(max_ind)      
 !        gens(i) = ceiling(xrt)
!5181  end do
      
            
!      end subroutine 
      
      
     
      
c      ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ran_srs(idum,num_to_draw,max_ind,gens)
            
      implicit real*8(a-h,o-z)
      

      integer num_to_draw,max_ind, gens(num_to_draw)
      integer idum,i
      real*8  rngs

!     ***** Generating *****

      do 5151 i = 1,num_to_draw
         rngs = ran2(idum)*dble(max_ind)      
         gens(i) = ceiling(rngs)
5151  end do
      
            
      end subroutine
      
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine readdata(fout,fdatname,
     &           lpsloop,ntot,split,nobs,nf,
     &           standard,
     &           kreg,regname,x,y)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8(a-h,o-z)
      
      integer fout,fdat
      parameter (fdat=20)
      parameter (maxk=104,maxn=88)
      character*12 regname(maxk),fdatname
      real*8 y(maxn),x(maxn,maxk),xsum(maxk),xssq(maxk)
      logical lpsloop,standard
c      
      open(unit=fdat, file=fdatname)
c
c  read kreg and nobs from *.dat file
c      
      read(fdat,*) ntot
      write(fout,*) '.. total nobs.........................',ntot
      write(*,*)    '.. total nobs.........................',ntot
      read(fdat,*) kreg
      write(fout,*) '.. total kreg.........................',kreg
      write(*,*)    '.. total kreg.........................',kreg
 
      if (kreg.gt.104) then
         write(fout,*) '!!!!! ERROR: kreg > 104 !!!!'
         write(*,*)    '!!!!! ERROR: kreg > 104 !!!!'  
         stop
      endif

      if (lpsloop) then 
         nf = int((1.0d0 - split)*ntot)
         nobs = ntot - nf
        write(fout,*) '.. LPS: will use nf .................',nf
        write(*,*) '.. LPS: will use nf .................',nf
      else
         nf = 0
         nobs = ntot
      endif           
c
c  read regressors names
c      
      do 20 i=1,kreg
       read(fdat,'(a12)') regname(i)
 20   enddo
c
c  init stuff
c 
      sumy = 0.0d0
      ssqy = 0.0d0
      do 50 j=1,kreg
         xsum(j) = 0.0d0
         xssq(j) = 0.0d0
 50   enddo
c
c  read data
c   
      do 100 i=1,ntot
         read(fdat,*) y(i),(x(i,j), j=1,kreg)
            sumy = sumy + y(i)
            ssqy = ssqy + y(i)**2
            do 70 j=1,kreg
               xsum(j) = xsum(j) + x(i,j)
               xssq(j) = xssq(j) + x(i,j)**2
 70         enddo
 100  enddo
      close(fdat)      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute means and vars
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      avey = sumy / dble(ntot)
      vary = ssqy / dble(ntot) - avey**2
      
      ssqy = ssqy - dble(ntot)*avey**2
      
      do 110 j=1,kreg
         xsum(j) = xsum(j)/dble(ntot)
         xssq(j) = xssq(j)/dble(ntot) - xsum(j)**2
110   enddo
 
      write(fout,*)
      write(fout,114) 
 114  format(t35,'Mean        St.Dev.')
      call barra('_',55,fout)
      write(fout,120) avey,dsqrt(vary)
      write(fout,*)
      do 115 j=1,kreg
         write(fout,122) j,regname(j),xsum(j),dsqrt(xssq(j))
 115  enddo
 120  format('Dep variable is Growth',t31,2g12.4)
 122  format('X(',i2,'): ',a12,t31,2g12.4)
      write(fout,*)
      call barra('_',55,fout)
      write(fout,*)
      write(*,*) '... Done w/ readdata ... '
      
      if (standard) then
      write(fout,*) 
      write(fout,*) '*******************************'
      write(fout,*) '**  Xs will be standardised  **'
      write(fout,*) '*******************************'
      write(fout,*) 
      do 200 i=1,nobs
      do 200 j=1,kreg
       x(i,j) = (x(i,j) - xsum(j))/dsqrt(xssq(j))
 200  enddo
      endif
c
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine splitdata(fout,regname,idum,ntot,nobs,nf,kreg,
     &     x,z,ztz,yin,yout,avey,ssqyn,lpsloop,yf,zf,fail)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit real*8(a-h,o-z)
      
      integer fout,fdat
      parameter (fdat=20)
      parameter (maxk=104,maxn=88,maxnf=50)
      integer sample(maxnf)
      character*12 regname(maxk)
      logical fail,lpsloop,inference(maxn)
c
c  note that zf is transposed of what it would normally be because later
c  this is more convenient in cholesky in LPS
c
      real*8 yin(maxn),yout(maxn),x(maxn,maxk),
     & z(maxn,maxk),ztz(maxk,maxk),avez(maxk),ssqz(maxk),
     & yf(maxnf),zf(maxk,maxnf)
c
c  init stuff
c 
      fail = .false.
      inobs = 0
      nnf = 0 
      
      sumy = 0.0d0
      ssqy = 0.0d0 
          
      do 10 j=1,kreg
         avez(j) = 0.0d0
         ssqz(j) = 0.0d0
  10  enddo
      do 20 i=1,nf
  20  sample(i) = 0    
c
c  if LPS draw nf out of nobs
c   
      if (lpsloop) then
	      call drawsample(idum,ntot,nf,sample,fail)
	      if (fail) then
		      write(*,*) ' Fail in drawsample!'
		      write(fout,*) ' Fail in drawsample!'
		      return
	      endif
	      do 50 i=1,ntot
	        inference(i) = .true.
  50          enddo
	      do 55 i=1,nf
	        inference(sample(i)) = .false.
  55          enddo 	       
      endif     
c
c  split sample
c  
      if (lpsloop) then
         do 101 i=1,ntot
          if (inference(i)) then          
c    if case i is in inference sample:          
               inobs = inobs +1
               yout(inobs) = yin(i)
               do 91 j=1,kreg
                     z(inobs,j) = x(i,j)
                     avez(j) = avez(j) + z(inobs,j)
                     ssqz(j) = ssqz(j) + z(inobs,j)**2
  91            enddo
                sumy = sumy + yout(inobs)
                ssqy = ssqy + yout(inobs)**2
             else
c   otherwise case i is part of prediction sample             
                nnf = nnf + 1
                yf(nnf) = yin(i)
c
c  note that zf is transposed of what it would normally be because later
c  this is more convenient in cholesky in LPS
c
                do 95 j=1,kreg
                     zf(j,nnf) = x(i,j)
 95             enddo
           endif
 101     continue 
c if no sample split:   
      else     
         do 125 i=1,ntot
            yout(i) = yin(i)
            sumy = sumy + yout(i)
            ssqy = ssqy + yout(i)**2
            do 111 j=1,kreg
               z(i,j) = x(i,j)
               avez(j) = avez(j) + z(i,j)
               ssqz(j) = ssqz(j) + z(i,j)**2
  111       enddo
  125     enddo
      endif
      
      close(fdat)     
           
      if (lpsloop) then
         call barra('.',55,fout)
         write(fout,*) ' Cases used for prediction: '
         write(fout,*)  (sample(i),i=1,nf)
         write(fout,*)
         call barra('.',55,fout)
      endif
      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  compute means of nobs = (nobs - nf) inference sample
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      avey = sumy/dble(nobs)
      vary = ssqy/dble(nobs) - avey**2
      ssqyn = ssqy - dble(nobs)*avey**2
      
      do 110 j=1,kreg
         avez(j) = avez(j)/dble(nobs)
         ssqz(j) = ssqz(j)/dble(nobs) - avez(j)**2
110   enddo
c
c  de-mean z's
c
      do 191 i=1,nobs
         do 190 j=1,kreg
            z(i,j) = z(i,j) - avez(j)
 190     enddo 
 191  enddo 
c
c  de-mean z_f's
c
      do 193 i=1,nf
         do 192 j=1,kreg
            zf(j,i) = zf(j,i) - avez(j)
 192     enddo
 193  enddo
c
c  construct z'z -> ztz
c      
      do 500 i=1,kreg
         do 400 j=1,kreg
            ztz(i,j) = 0.0d0
            do 300 n=1,nobs
               ztz(i,j) = ztz(i,j) + z(n,i)*z(n,j)
 300        enddo
 400     enddo
 500  enddo
c
c  if sample was split re-write means and StDev
c 
      if (lpsloop) then
       write(fout,*)  
       write(fout,*) ' LPS will be called'
       write(fout,*) ' Sub-sample used inference: '
       write(fout,*) '  Nobs: ',nobs

       write(fout,114) 
 114  format(t35,'mean        St.Dev.')
        call barra('_',55,fout)
       write(fout,120) avey,dsqrt(vary)
       write(fout,*)
       do 115 j=1,kreg
         write(fout,122) j,regname(j),avez(j),dsqrt(ssqz(j))
 115   enddo
 120   format('Dep variable is Growth',t31,2g12.4)
 122   format('X(',i2,'): ',a12,t31,2g12.4)
       write(fout,*)
       call barra('_',55,fout)
       write(fout,*)
      endif
      
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine check_top_mod_post                                   ccc
ccc                                                                  ccc
ccc  modified: Shutong Ding, 2015-06-01                              ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    
      subroutine check_top_mod_post(visits,midx,imax,top_num,top_midx,
     & top_mod_post )
       
      implicit real*8(a-h,o-z)
      integer maxm, max_top, top_num, imax, ix, jx
      parameter(maxm = 200000, max_top = 100)
      real*8 visits(maxm),midx(maxm,2),top_midx(max_top,2), 
     & top_mod_post(max_top)
      
      ! initializing, if not visited, default set to 0
      top_mod_post = 0d0
      
      do 5195 jx = 1,top_num
      do 5191 ix = 1,imax
          ! search over the all visited models to find matches
          if ((midx(ix,1) .eq. top_midx(jx,1)) .and. (midx(ix,2) .eq. 
     &    top_midx(jx,2)) ) then   ! both indexes matches
             top_mod_post(jx) = visits(ix)      
          end if    
5191  continue    
5195  end do
      
      end subroutine check_top_mod_post
      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine wrchainfo                                            ccc
ccc                                                                  ccc
ccc  modified: july 06                                               ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     modified to check the top listed models, sdg
       
      subroutine wrchainfo(regname,kreg,visits,midx,bstar,g0j,kj,
     &      imax,nvout,ibm,icut,indx,fout,justcut,beta)

      implicit real*8(a-h,o-z)
      parameter(maxm=200000,thres=0.25d0,maxk=104) 
      logical justcut
      integer model(maxk),ir(maxk),fout,kjj,
     &      indxx(maxm),indx(maxm),kj(maxm)
      real*8 visits(maxm),beta(maxk),bb(maxk),bbb(maxk),
     &      bstar(maxk,maxm),g0j(maxm),midx(maxm,2),
     &      vkj(0 : maxk)
     
      character*12 regname(maxk)
      character*25 height
c 
      height='XXXXXXXXXXXXXXXXXXXXXXXXX'
c      
      write(fout,*)
      write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%%'   
      write(fout,*) '%%%   WrChainInfo    %%%'
      write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%%'  
      write(fout,*)      
c
c  initialize ---careful with data statements as opt can screw things up
c   
      vkj(0) = 0.0d0    
       do 10 i=1,kreg
      vkj(i) = 0.0d0
  10  beta(i)= 0.0d0
c
c here we go...
c 
      write(*,*) '... wrchainfo: visited models ...'      
c     
c  write chain info: visited models
c
      write(*,*)    '... Number of models visited is:', imax
      write(fout,*) 'Number of models visited is:', imax
      write(fout,*)
      if (nvout.gt.0) write(fout,*) nvout,' visits out!!!'            
c
c  sort according to number of visits w/ auxiliary array indxx
c  which contains indexes of models in _ascending_ order
c  
      write(*,*) '... wrchainfo: sorting ...'      
      write(*,*) '... imax',imax
      call indexx(imax,visits,indxx)
      tot=0.0d0
c
c reset indx to descending order
c
      do 200 i=1,imax
 200  indx(i)=indxx(imax+1-i) 
c
c  write out post mass of best models
c     
      write(fout,*)
      call barra('-',70,fout)
      cumass_1 = 0.0d0

      do 420 i=1,imax
         cumass = cumass_1 + visits(indx(i))
	 if ((i.eq.5).or.(i.eq.25).or.(i.eq.50).or.
     &       (i.eq.75).or.(i.eq.100)) then
            write(fout,421) i,cumass*100.0d0
	 endif
	 if ((cumass_1.lt.0.25d0).and.(cumass.ge.0.25d0)) then
            i25=i
	    write(fout,421) i25,cumass*100.0d0
	 endif
	 if ((cumass_1.lt.0.5d0).and.(cumass.ge.0.5d0)) then
            i50=i
	    write(fout,421) i50,cumass*100.0d0
	 endif
	 if ((cumass_1.lt.0.75d0).and.(cumass.ge.0.75d0)) then
            i75=i
	    write(fout,421) i75,cumass*100.0d0
	 endif
	 if ((cumass_1.lt.0.9d0).and.(cumass.ge.0.9d0)) then
            i90=i
	    write(fout,421) i90,cumass*100.0d0
	 endif
	 if ((cumass_1.lt.0.95d0).and.(cumass.ge.0.95d0)) then
            i95=i
	    write(fout,421) i95,cumass*100.0d0
	 endif
	 if ((cumass_1.lt.0.99d0).and.(cumass.ge.0.99d0)) then
            i99=i
	    write(fout,421) i99,cumass*100.0d0
	 endif
	 cumass_1=cumass
 420  enddo       
      write(fout,421) imax,cumass*100.0d0
 421  format('best ',i6,' models account for ',f8.1,' % of mass')
      icut=i90
 
      call barra('-',70,fout)
      write(fout,*)
      if (justcut) then
        return
      endif
c-----------------------------------------------------------------------
c  loop over visited models to compute average number of regressors
c-----------------------------------------------------------------------
      write(fout,*)
      call barra('=',55,fout) 
      write(fout,*) ' number of variables in models'
      call barra('_',55,fout) 
      write(fout,*) '      nmod      postprob      ave         std'
      call barra('_',55,fout) 
      sumvisit = 0.0d0
      sumreg   = 0.0d0
      ssqreg   = 0.0d0
      
      do 450 i=1,imax
c
c 1. find # regressors in each model
c         
         imd = indx(i)
         call g2model(midx(imd,1),midx(imd,2), kreg, model)
         nreg = kj(imd)
c
c  vkj() is an array 0..kreg where we store the visits to models with
c  0, 1, 2, ...kreg regressors
c         
         vkj(nreg) = vkj(nreg) + visits(imd)         
c
c 2. accumulate stats
c 
         sumvisit = sumvisit + visits(imd)
         sumreg=sumreg + visits(imd)* dble(nreg)
         ssqreg=ssqreg + visits(imd)* dble(nreg**2)
         
         if ((i.eq.5).or.(i.eq.10).or.
     &         (i.eq.25).or.(i.eq.50).or.
     &         (i.eq.75).or.(i.eq.100).or.
     &         (i.eq.i25).or.(i.eq.i50).or.
     &         (i.eq.i75).or.(i.eq.i90).or.
     &         (i.eq.i95).or.(i.eq.i99)) 
     &         then
         ave = sumreg/sumvisit
         std = dsqrt(ssqreg/sumvisit - ave**2)        
         write(fout,451) i,sumvisit*100.0d0,ave, std     
      endif
 450  enddo
      ave = sumreg/sumvisit
      std = dsqrt(ssqreg/sumvisit - ave**2)         
      write(fout,451) imax,sumvisit*100.0d0,ave, std
 451  format(i10,3f12.1)
c------------------------------------------
c end of loop over models
c------------------------------------------
      write(fout,*)
      call barra('-',55,fout)
      
         write(fout,*) ' Mass by model size '
         write(fout,*) '--------------------'  
      
      xvkj = 0
      
      do 666 i=0,kreg
         if (vkj(i).gt.xvkj) xvkj=vkj(i)
 666  enddo
      
      s=0.0d0
      sum = 0.0d0
      ssq = 0.0d0
      
      do 710 i=0,kreg
         sum=sum+dble(i)*vkj(i)
         ssq=ssq+(dble(i)**2)*vkj(i)      
         if (vkj(i).ge.1.d-4) then      
            write(fout,711) height(1:int(25.0d0*vkj(i)/xvkj)),
     &            i,vkj(i)*100
         endif
         s=s+vkj(i)
  710 enddo
      
      
  711 format(a25,t27,i3,f8.3)

      call barra('-',55,fout)
      write(fout,'(t5,a24,f8.3)') 'Average model size.....',sum
      write(fout,'(t5,a24,f8.3)') 'st.dev. of model size..',
     &      dsqrt(ssq-sum**2)
      
      
      
      write(fout,'(t5,a7)') 'Prior:'
         theta=0.5d0
         write(fout,'(t7,a24,f8.3)') 
     &         'Expected model size....',theta*kreg
         write(fout,'(t7,a24,f8.3)') 
     &         'st.dev.  model size....',
     &         dsqrt(theta*kreg*(1.0d0-theta))      
      
      write(fout,*)    
      write(fout,'(t5,a24,f8.3)') 'This should be one.....', s
      call barra('-',55,fout)
      write(fout,*)
c
c prior prob
c                   
      write(fout,'(a33,t45,e12.4,x,a1)')
     &      'prior prob for a single model is',
     &      (1.0d0/dble(2.0d0**kreg))*100.0d0,'%'
c
c  write models with post prob larger than thres 
c                     
      write(fout,'(a41,t45,f8.4,x,a1)') 
     &      'models with post probability larger than',thres,'%'
      write(fout,*)
      call barra('-',70,fout)
      write(fout,*)
      write(fout,*) '   postprob        regressors'
      call barra('-',70,fout)
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop over visited models
ccc       
      do 500 i=1,imax   
         percvisits = visits(indx(i))*100.0d0
         call g2model2(midx(indx(i),1),midx(indx(i),2),
     &         kreg,model,kjj,ir)
c
c  accumulate post mass for the betas
c               
         do 460 ib=1,kjj
 460     beta(ir(ib)) = beta(ir(ib)) + percvisits
c
c  write out model post prob and regressors for all
c  models w/ post prob >= thres
c               
         if (percvisits.ge.thres) then
            tot = tot + percvisits
            write(fout,467) i, percvisits,(ir(k),k=1,kjj)
 467        format(i3,t5,f6.2,'%',x,'|',100(:i3))
         endif
 500  enddo
ccc
c end of loop over visited models
cccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c write post prob accounted for the models above
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(fout,*)
      write(fout,510) tot
 510  format(t4,f6.2,'%')
      write(fout,*)
      call barra('-',70,fout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c best model
cccccccccccccccccccccccccccccccccccccccccccccccccccccc       
      indice = indx(1)
      ibm = indice     
      call g2model2(midx(indice,1),midx(indice,2), 
     &      kreg, model,kjj,ir)   
      do 550 j=1,kreg
 550  bb(j) = 0.0d0
      do 555 j=1,kjj
         bb(ir(j)) = bstar(j,indice)/(1.0d0 + g0j(indice) )
 555  continue

c
c store next best-model betas in bbb()
c  
      if (imax.gt.1) then
         indice = indx(2)     
         call g2model2(midx(indice,1),midx(indice,2), 
     &         kreg, model,kjj,ir)      
         do 556 j=1,kreg
            bbb(j) = 0.0d0
 556     continue
         do 557 j=1,kjj
 557     bbb(ir(j)) = bstar(j,indice)/(1.0d0 + g0j(indice) )
      endif
      
c
c  post probability of incl: regressors
c            
      write(fout,*)
      write(fout,590) 
 590  format('post prob of incl and 2 best models m_j''s')
      call barra('_',70,fout)
      write(fout,*)
      do 600 i=1,kreg
 600  write(fout,610) i,regname(i),beta(i),bb(i),bbb(i)
 610  format(i2,1x,a12,2x,f6.2,' %',2x,2g14.4)
      write(fout,*)
      call barra('=',70,fout)
c     
                
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc  subroutine wrpostinfo
ccc
ccc  modified 
ccc    06.07.19 post variance of betas fixed
ccc    06.05.11 to compute ave # regressors in models
ccc    06.04.14 to simply write out means and sds of 
ccc             continuous part of posterior distribution of betas
ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine wrpostinfo(regname,nobs,kreg,ztz,z,
     &      prob,g0j,dstar,bstar,midx, 
     &      imax,fout,fail)

      implicit real*8(a-h,o-z)
      parameter(maxm=200000,maxk=104,maxn=88)  
      logical fail
      integer model(maxk),ir(maxk),fout,ir2(maxk)     
      character*12 regname(maxk)
      
      real*8 prob(maxm),beta(maxk),ztz(maxk,maxk),
     &      ztzj(maxk,maxk),z(maxn,maxk),zj(maxn,maxk),
     &      zj2(maxn,maxk),ztzj2(maxk,maxk),
     &      c(maxn),cz(maxk),cz2(maxk),p(maxk),
     &      bstar(maxk,maxm),dstar(maxm),midx(maxm,2),
     &      g0j(maxm),
     &      bvar(maxk),bmean(maxk),bssq(maxk),
     &      bvarc(maxk),bmeanc(maxk),bsign(maxk)            
c
      write(fout,*)
      write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%'   
      write(fout,*) '%%%   wrpostinfo    %%%'
      write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%' 
c
c Initialize stuff
c
      fail = .false.
   
      do 455 i=1,kreg
       bmean(i) = 0.0d0
       bvar(i) = 0.0d0
       bssq(i) = 0.0d0
       bmeanc(i) = 0.0d0
       bvarc(i) = 0.0d0
       beta(i) = 0.0d0
       bsign(i) = 0.0d0
 455  continue
c
c Accumulate Post Prob of Including Betas
c            
      do 500 imd=1,imax
         call g2model2(midx(imd,1),midx(imd,2),kreg,model,kj,ir)
         do 460 ib=1,kj
  460      beta(ir(ib)) = beta(ir(ib)) + prob(imd)
  500  enddo 
c
c  Posterior moments of the betas
c  ===============================
c  See FEDEA WP www.fedea.es/pub/Papers/1997/dt97-25.pdf
c  eq (4.5) on page 13 for the expression of the student-t.
c  Complication arises because we must construct c'(I - F(F'F)^{-1}F')c
c  where c is the l-th column of Z, and F is Z without such column...
c  
      nu = nobs-1
ccccccccccccccccccccccccccccccccccccc
c  Loop over visited models
ccccccccccccccccccccccccccccccccccccc         
      do 1010 imd=1,imax            
         vb = dstar(imd)/(dble(nu)*(g0j(imd)+1.0d0))
         call g2model2(midx(imd,1),midx(imd,2),kreg,model,kj,ir)         
c
c  for each visited model construct (z_j'z_j) and z_j
c         
         do 730 j=1,kj
            do 720 i=1,j
 720        ztzj(i,j) = ztz(ir(i),ir(j))
            do 725 i=1,nobs
 725        zj(i,j) = z(i,ir(j))
 730     enddo
         
         k = kj-1
c
c compute the variance for each regressor jreg
c         
         do 1000 jreg=1,kj
c
c store in ir2() indices 1..jreg-1, jreg+1,...kj
c         	                 
            do 805 j=1,jreg-1
 805        ir2(j) = j	    
            do 806 j=jreg,kj
 806        ir2(j) = j+1	    
c
c construct submatrices of ztzj and zj in ztzj2 and zj2
c by dropping the column (row) corresponding to jreg
c 	     
            do 820 j=1,k
               do 810 i=1,j
 810           ztzj2(i,j) = ztzj(ir2(i),ir2(j))
               do 815 i=1,nobs
 815           zj2(i,j) = zj(i,ir2(j))
 820        enddo
c
c  c: extract corresponding column from z
c  c2 = c'c
c 
            c2 = 0.0d0
            do 840 i=1,nobs
               c(i) = zj(i,jreg)
               c2 = c2 + c(i)**2
 840        continue
c
c cz: is a (k x 1) column vector given by (Zj2')(c)
c
            do 880 j=1,k
               cz(j) = 0.0d0
               do 850 i=1,nobs
                 cz(j) = cz(j) + c(i)*zj2(i,j)
 850           enddo          
 880        enddo
c
c  Must construct c'[I - F(F'F)^{-1}F']c = c'c - c'F(F'F)^{-1}F'c
c  with: 
c      F'c := cz = (Zj2')c
c      F'F := ztzj2
c
c  choldc + cholsc will return 
c     cz2 = (F'F)^{-1}F'c = (ztzj2)^{-1}cz
c
            call choldc(ztzj2,k,maxk,p,fail)
            if (fail) then
             write(fout,*) jreg, 'fail choldc!!!'
             write(fout,*) k,p 
             return
            endif
            call cholsl(ztzj2,k,maxk,p,cz,cz2)
c
c  now we compute x:= c'F(F'F)^{-1}F'c = c'F [cz2]
c            
            x = 0.0d0
            do 900 j=1,k
 900        x = x + cz(j)*cz2(j)
c
c   finally  x:= c'c - c'F(F'F)^{-1}F'c
c  (added: 06.07.19)
c            
            x =  c2 - x     
c       
            if (x.ne.0) then
            vvb = (vb/x)*(dble(nu)/(dble(nu)-2.0d0))
            else
             call barra('!',50,fout) 
             write(fout,*) 'x eq 0!!!! regr #', ir(jreg),
     &                     'prob:',prob(imd) 
              write(fout,*) 'imd:',imd
              write(fout,*) (model(i),i=1,kreg)
              call barra('!',50,fout)
             goto 999
            endif            
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c #1 is the \sum_i p_i \sigma^2_{j | M_i}
c    (weighted average of the post variances of betaj for ea Model Mi)
c	    
	    bvar(ir(jreg)) = bvar(ir(jreg)) + vvb*prob(imd)
c
c #2 is the \sum_i p_i \mu^2_{j | M_i}
c    (weighted average of the post means**2 of betaj for ea Model Mi)
c	         
            bssq(ir(jreg)) = bssq(ir(jreg))
     &  + ((bstar(jreg,imd)/(g0j(imd)+1.0d0))**2)*prob(imd)         
c
c #3 is the \sum_i p_i \mu_{j | M_i}
c    (weighted average of the post means of betaj for ea Model Mi)
c	    
            bmean(ir(jreg)) = bmean(ir(jreg))
     &      +(bstar(jreg,imd)/(g0j(imd)+1.0d0))*prob(imd)
     
            if (bstar(jreg,imd).ge.0.0d0) then
              bsign(ir(jreg))=bsign(ir(jreg))+prob(imd)
            endif
  999    continue
 1000    enddo
        
 1010 enddo
c 
c..............................................................
c  var(beta_i) =  \sum_i p_i \sigma^2_i
c               + \sum_i p_i (\mu_i - \mu)^2
c               + (1-\sum_i p_i) *\mu^2  [this corresponds to the mass at 0]
c
c  the first part is already stored in varb2, the 2nd equals:
c
c        \sum_i p_i \mu_i^2 - \mu^2
c
c...............................................................
c
c first we write out the parts and only then we store the
c variance in bvar(i)
c 
         write(fout,*) 'Beta: Post Prob of Inclusion'
         do 1500 i=1,kreg  
            write(fout,'(i3,2x,e20.10)') i,beta(i)       
            bvar(i) = bvar(i) + bssq(i) - bmean(i)**2 
c
c if we want the moments of the continuous part only,
c must rescale since beta(i) may be < 1.
c             
            if (beta(i).gt.0) then
               bsign(i) = bsign(i)/beta(i)        
               bsign(i) = 100.0d0*max(bsign(i),1.0d0-bsign(i))               
               bmeanc(i) = bmean(i)/beta(i)
               bvarc(i) = bvar(i)/beta(i) +
     &               (beta(i)-1.0d0)*bmeanc(i)**2
            else
               bsign(i)  = -99.9d0
               bmeanc(i) = -99.9d0
               bvarc(i)  = -99.9d0
            endif       
 1500    enddo
c
c  post probability of incl: regressors
c            
       write(fout,*)
       write(fout,*) ' Betas: Posterior Moments (Unconditional and Condi
     &tional on Inclusion)'
      write(fout,*)
       write(fout,1550)
 1550 format(t13,'p[i]',3x,'Mean[i]',5x,'sd[i]',5x,'|m/sd|',
     & 2x,'Mean[i|I]',2x,'sd[i|I]',3x,'|m/sd|',1x,'Sign')       
       call barra('_',85,fout) 
       do 3000 i=1,kreg
          if (bvar(i).gt.0) then
             x1 = dsqrt(bvar(i))
             x2 = abs(bmean(i)/dsqrt(bvar(i)))
          else
             x1 = -99.9d0
             x2 = -99.9d0
          endif
          if (bvarc(i).gt.0) then
             z1 = dsqrt(bvarc(i))
             z2 = abs(bmeanc(i)/dsqrt(bvarc(i)))
          else
             z1 = -99.9d0
             z2 = -99.9d0
          endif          
          write(fout,3400) i,regname(i),beta(i)*100.0d0, 
     &          bmean(i), x1, x2, bmeanc(i), z1, z2, bsign(i)                
 3000  enddo
       call barra('_',85,fout) 
       write(fout,*)
 3400 format(i2,1x,a8,1x,f5.1,2g12.4,f5.1,2g12.4,f5.1,1x,f5.1)         
c
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine convergence                                          ccc
ccc                                                                  ccc
ccc  created: 06.07.19                                               ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
       subroutine convergence(regname,
     &   iprior,idum,initrep,mnumrep,nobs,kreg,
     &   ztz,z,y,ssqyn,
     &   freq,bayesf,gj,dstar,bstar,midx,kj,
     &   imax,nvout,icut,fout,fail,fjout,mmodel,fjsum,fix_start)   

       implicit real*8(a-h,o-z)
       parameter(maxm=200000,maxk=104,maxn=88)
       
       integer fout     
       integer iprior,idum,initrep,mnumrep,nobs     
       logical fail,visited,allvis,fix_start
       integer kj(maxm),mmodel(maxk),idx(maxm),ir(maxk)
       character*12 regname(maxk)

       real*8 z(maxn,maxk), y(maxn), ztz(maxk,maxk),beta(maxk),
     & midx(maxm,2),midx2(maxm,2),fjout,
     & bayesf(maxm),freq(maxm),bayesf2(maxm),freq2(maxm),
     & gj(maxm),dstar(maxm),bstar(maxk,maxm)

       initrep=min(100000,initrep/3)
       mnumrep=min(500000,mnumrep/3)
             
       write(fout,*)
       write(fout,*) ' 2nd G&M chain: ',initrep,mnumrep
       write(fout,*)
       write(*,*)
       write(*,*)    ' ... 2nd G&M chain: ',initrep,mnumrep
       write(*,*)
c
c  things we don't want to overwrite:
c
c     freq2, bayesf2, midx2, imax2, fjsum2
c     
       call runchain(fix_start,iprior,idum,initrep,mnumrep,nobs,kreg,
     &       ztz,z,y,ssqyn,
     &       freq2,bayesf2,gj,dstar,bstar,midx2,kj,
     &       imax2,nvout,fout,fjout,fjsum2,fail)     

       write(*,*) '... G&M: done w/ chain ...'
       
       call wr_date(fout,.true.,.true.)
       
       if (fail) then
          write(fout,*) ' runchain failed in G&Mstuff!!!'
          write(*,*)    ' runchain failed in G&Mstuff!!!'
          return
       endif
      
       call wrchainfo(regname,kreg,bayesf2,midx2,bstar,gj,kj,
     &      imax2,nvout,ibm,icut,idx,fout, .false.,beta)
          
      write(*,*) '... G&M: done w/ wrchaininfo ...'
      call wr_date(fout,.true.,.true.)
      write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'   
      write(fout,*) '%%%        Convergence         %%%'
      write(fout,'(a15,i5,a15)') ' %%%   Will use ',icut,
     &                           '    models  %%%'
      write(fout,*) '%%%  90% of mass of 1st chain  %%%'
      write(fout,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'  
         
      if (fail) then
         write(fout,*) ' wrchainfo failed in G&Mstuff!!!'
         write(*,*)    ' wrchainfo failed in G&Mstuff!!!'
         return
      endif
c
c  icut will tell us how many models in this 2nd chain are needed
c  to account for 90% of posterior mass
c
c  idx contains the sorted array-indices of the top icut models
c    
      fja  = 0.0d0
      fjaa = 0.0d0
      fjaa1= 0.0d0
      visa = 0.0d0
c
c  now we just look at the icut best models (whose ind are the top
c  icut entries in idx()) of 2nd chain
c
         allvis =.true.
         inv=0
         xmass=0.0d0
         do 100 i=1,icut
          visited = .false.
c
c  accumulate marginal likelihoods (need to multiply by ml(m) since 
c  bayesf are normalized numbers
c
           fja= fja + bayesf2(idx(i))*fjsum2
c
c  loop over visited models in chain #1, if it was visited there,
c  accumulate bayesf in fjaa and prob in visa
c
          do 8625 imd=1,imax                                           
            if ((midx2(idx(i),1).eq.midx(imd,1)).and.
     &          (midx2(idx(i),2).eq.midx(imd,2))) then                           
                 fjaa = fjaa  + bayesf2(idx(i))*fjsum2
                 fjaa1= fjaa1 + bayesf(imd)*fjsum
                 visa = visa + freq(imd)
                 visited = .true.
                 goto 8626
            endif
 8625    enddo
 8626    continue
         if (.not.visited) then
          inv=inv+1
          xmass = xmass + bayesf2(idx(i))*100.0d0  
	  allvis = .false.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                              c
c             call g2model2(midx2(idx(i),1),midx2(idx(i),2),   c
c        &                  kreg,mmodel,kkjj,ir)               c
c             write(fout,101) inv,bayesf2(idx(i))*100.0d0      c
c             write(fout,*) 'idx1: ', midx2(idx(i),1)          c
c             write(fout,*) 'idx2: ', midx2(idx(i),2)          c
c             write(fout,*) '# reg: ',kkjj                     c
c             write(fout,*) (ir(j),j=1,kkjj)                   c
c                                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         endif
  100    enddo
  101    format(i3,'. Model not visited in 1st chain, mass: ',e12.4,'%') 
         write(fout,*) 
         write(fout,*) 'Mass of not visited %: ', xmass
         if (allvis) write(fout,*)'all visited !'
c
c  write results
c
      write(fout,*)
      call barra('=',75,fout)
      write(fout,*)
      write(fout,*) 'George & McCulloch: direct estimation'
      write(fout,*) 'post model prob of all visited models'
      write(fout,*)
      write(fout,20) visa
  20  format(' vis(A)................',e12.4)
      write(fout,30) fja
  30  format(' ml(A).................',e12.4)
c
c  we want fjsum of the 1st chain here!
c
      write(fout,40) fjsum
  40  format(' ml(M).................',e12.4)
      if (fja.ne.0) then
      write(fout,50) (visa*(fjsum/fja))*100.0d0
      else
      write(fout,50) -99.99
      endif
  50  format(' pmp(M)................',e12.4,' %') 
      if (fja.ne.0) then
      write(fout,60) (fjaa/fja)*100.0d0
      write(fout,61) (fjaa/fjaa1)*100.0d0
      else
      write(fout,50) -99.99
      endif
  60  format(' %A visited by chain...',e12.4,' %') 
  61  format(' this should be 100....',e12.4,' %') 
      write(fout,*)
      call barra('=',75,fout)
      write(fout,*)
c
      call wr_date(fout,.true.,.true.)
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine jointness                                            ccc
ccc                                                                  ccc
ccc  created: 06.07.19                                               ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        
      subroutine jointness(kreg,prob,midx,
     &      imax,fout,foutbj,fouttj,fail)

      implicit real*8(a-h,o-z)
      parameter(maxm=200000,maxk=104) 
      logical fail
      integer model(maxk),ir(maxk),fout,foutbj,fouttj,
     &      indxx(maxm),indx(maxm)
      real*8 prob(maxm),beta(maxk),
     &      midx(maxm,2),
     &      bjoint(maxk,maxk),bjoint1(maxk,maxk),bjoint2(maxk,maxk),
     &      tjoint(maxk,maxk,maxk),tjoint1(maxk,maxk,maxk),
     &      tjoint2(maxk,maxk,maxk)

      write(fout,*)
      write(fout,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%'   
      write(fout,*) ' %%%     Jointness     %%%'
      write(fout,*) ' %%%%%%%%%%%%%%%%%%%%%%%%%'  
      write(fout,*)      
c
c  initialize ---careful with data statements as opt can screw things up
c   
       do 11 i=1,kreg 
            beta(i)=0.0d0      
        do 11 j=1,kreg
            bjoint(i,j) =0.0d0
            bjoint1(i,j)=0.0d0
            bjoint2(i,j)=0.0d0 
            do 11 k=1,kreg
            tjoint(i,j,k) =0.0d0
            tjoint1(i,j,k)=0.0d0
            tjoint2(i,j,k)=0.0d0   
  11  enddo
      qjoint=0.0d0
      qunion=0.0d0
      q5joint=0.0d0
      q5union=0.0d0
c
c here we go...
c 
      fail=.false.
      write(*,*) '... jointness: visited models ...'      
c
c  sort according to number of prob w/ auxiliary array indxx
c  which contains indexes of models in _ascending_ order
c  
      write(*,*) '... wrchainfo: sorting ...'      
      write(*,*) '... imax',imax
      call indexx(imax,prob,indxx)
c
c reset indx to descending order
c
      do 200 i=1,imax
 200  indx(i)=indxx(imax+1-i) 
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c loop over visited models
cccccccccccccccccccccccccccccccccccccccccccccccccccccc 
          
      do 500 i=1,imax   
         percprob = prob(indx(i))*100.0d0
         call g2model2(midx(indx(i),1),midx(indx(i),2),
     &                 kreg,model,kjj,ir)


         do 431 ib=1,kjj-1
          do 431 jb=ib+1,kjj
          bjoint(ir(ib),ir(jb))=bjoint(ir(ib),ir(jb))
     &          + percprob
          do 431 kb=jb+1,kjj
          tjoint(ir(ib),ir(jb),ir(kb))=tjoint(ir(ib),ir(jb),ir(kb))
     &          + percprob   
 431     enddo
c
c if the 1st 4/5 vars are included then ir(4/5) must be 4/5
c 
         if (ir(4).eq.4) qjoint = qjoint + percprob
         if (ir(5).eq.5) q5joint = q5joint + percprob
c
c now the union of the first 4/5
c    
         do 440 ib=1,5
           if (ir(ib).le.5) then 
             qunion = qunion + percprob
             q5union = q5union + percprob
             goto 441
           endif
           if (ir(ib).le.4) then 
             qunion = qunion + percprob
             goto 441
           endif
 440     enddo 
 441     continue    

c
c  accumulate post mass for the betas
c               
         do 460 ib=1,kjj
 460        beta(ir(ib)) = beta(ir(ib)) + percprob
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c end of loop over visited models
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
 500  enddo
c.....................................................
c bi - jointness stuff
c.....................................................
c
       b1ave  = 0.0d0
       b1sdv  = 0.0d0
       b2ave  = 0.0d0
       b2sdv  = 0.0d0
       bcorr  = 0.0d0
       nobs   = 0
       b1ave1 = 0.0d0
       b1sdv1 = 0.0d0
       b2ave1 = 0.0d0
       b2sdv1 = 0.0d0
       bcorr1 = 0.0d0
       nobs1  = 0
       b1ave2 = 0.0d0
       b1sdv2 = 0.0d0
       b2ave2 = 0.0d0
       b2sdv2 = 0.0d0
       bcorr2 = 0.0d0
       nobs2  = 0
       idec = 0
       ist  = 0
       ivst = 0
       ifav = 0
c
c  loop over "connections"
c
      write(foutbj,*) ' i j pri prj prij bj bf'
      l=0
      lpos=0
      do 615 i=1,kreg
      do 615 j=i+1,kreg
        l=l+1
        if (bjoint(i,j).gt.0) then
           lpos=lpos+1
	          if ((beta(i)+beta(j)-bjoint(i,j)).gt.0d0) then
	           bjoint1(i,j)= 100.0d0*bjoint(i,j)            
     &             /(beta(i)+beta(j)-bjoint(i,j))
	          else
	           write(fout,*) 'error in bi-joint'
	           fail=.true.
	           return
	          endif
	          
	          if ((beta(i)+beta(j)-2.0d0*bjoint(i,j)).gt.0) then
	             bjoint2(i,j)= bjoint(i,j)
     &                  /(beta(i)+beta(j)-2.0d0*bjoint(i,j))
	          else
	             bjoint2(i,j)=9.99d9
	             write(fout,*) 'bi-joint ',i,j, 'denom = 0!'
	          endif
	     else
	          bjoint1(i,j)=0.0d0
	          bjoint2(i,j)=0.0d0
       endif	
c
c  write to ????bj.dat, note that r read.table likes e not d format!
c	
         write(foutbj,'(i6,2i3,5e15.8)') l,i,j,beta(i),beta(j),
     &                  bjoint(i,j),bjoint1(i,j),bjoint2(i,j)

cendif0	
	bj1 = bjoint1(i,j)
	bj2 = bjoint2(i,j)
	      b1ave=b1ave+bj1
	      b1sdv=b1sdv+bj1**2
	      b2ave=b2ave+bj2
	      b2sdv=b2sdv+bj2**2
	      bcorr=bcorr+bj1*bj2
	      nobs=nobs+1
       if (bjoint2(i,j).gt.1.0d0)  then
c  jointness bf > 1	      
	      b1ave1=b1ave1+bj1
	      b1sdv1=b1sdv1+bj1**2
	      b2ave1=b2ave1+bj2
	      b2sdv1=b2sdv1+bj2**2
	      bcorr1=bcorr1+bj1*bj2
	      nobs1=nobs1+1
        else
c disjointness bf < 1	      
	      b1ave2=b1ave2+bj1
	      b1sdv2=b1sdv2+bj1**2
	      b2ave2=b2ave2+bj2
	      b2sdv2=b2sdv2+bj2**2
	      bcorr2=bcorr2+bj1*bj2
	      nobs2=nobs2+1	
c for disjointness accumm	
	       if (bj2.lt.(1.0d-2)) then
	           idec = idec + 1
	          else
	              if (bj2.lt.(1.0d0/30.0d0)) then
	                  ivst = ivst+1
	               else
	                 if (bj2.lt.(1.d-1)) then
	                    ist = ist + 1
	                  else
	                    if (bj2.lt.(1.0d0/3.0d0)) ifav = ifav+1
	                 endif
	              endif
	        endif
                
         endif
 615  enddo
       write(fout,'(a11,i4,a8,i5,a20,f6.2,a1)') ' There are ',lpos,
     &' out of ', l, 
     &' pij greater than 0:',100.*real(lpos)/real(l),'%'
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
       b1ave = b1ave/dble(nobs)
       b1sdv = dsqrt(b1sdv/dble(nobs) - b1ave**2)
       b2ave = b2ave/dble(nobs)
       b2sdv = dsqrt(b2sdv/dble(nobs) - b2ave**2)
       bcorr = dsqrt((bcorr/dble(nobs) - b1ave*b2ave)/(b1sdv*b2sdv))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       b1ave1 = b1ave1/dble(nobs1)
       b1sdv1 = dsqrt(b1sdv1/dble(nobs1) - b1ave1**2)
       b2ave1 = b2ave1/dble(nobs1)
       b2sdv1 = dsqrt(b2sdv1/dble(nobs1) - b2ave1**2)
       bcorr1 = dsqrt((bcorr1/dble(nobs1)-b1ave1*b2ave1)       
     &                 /(b1sdv1*b2sdv1))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
       b1ave2 = b1ave2/dble(nobs2)
       b1sdv2 = dsqrt(b1sdv2/dble(nobs2) - b1ave2**2)
       b2ave2 = b2ave2/dble(nobs2)
       b2sdv2 = dsqrt(b2sdv2/dble(nobs2) - b2ave2**2)
       bcorr2 = dsqrt((bcorr2/dble(nobs2)-b1ave2*b2ave2)
     &                 /(b1sdv2*b2sdv2))       
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc       
      
      
      
      write(fout,*) 
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*) ' x   Bi-jointness measures    x'
      write(fout,*) ' x   descriptive statistics   x'
      write(fout,*) ' x   .....................    x'  
      write(fout,*) ' x   +/- means PO .gt/lt. 1   x' 
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*)     
      write(fout,*)  '            bj1         bj1+        bj1-'
      write(fout,616) ' mean',b1ave,b1ave1,b1ave2
      write(fout,616) ' sd  ',b1sdv,b1sdv1,b1sdv2
      write(fout,*)     
      write(fout,*)  '            bj2         bj2+        bj2-'
      write(fout,616) ' mean',b2ave,b2ave1,b2ave2
      write(fout,616) ' sd  ',b2sdv,b2sdv1,b2sdv2
      write(fout,*)
      write(fout,616) ' corr',bcorr,bcorr1,bcorr2
      write(fout,617) ' nobs',nobs,nobs1,nobs2
      write(fout,*)
 616  format(a6,3f12.3)
 617  format(a6,3i12)
      itot = ifav + ist + ivst + idec
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*) ' x   Bi- disjointness by BF ratio   x'
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*)
      write(fout,'(a22,i4,f8.2,a2)') ' favourable    (1/3) ', ifav, 
     &			100.0*real(ifav)/real(itot),'%'               
      write(fout,'(a22,i4,f8.2,a2)') ' strong       (1/10) ', ist, 
     &			100.0*real(ist)/real(itot),'%' 
      write(fout,'(a22,i4,f8.2,a2)') ' very strong  (1/30) ', ivst, 
     &			100.0*real(ivst)/real(itot),'%' 
      write(fout,'(a22,i4,f8.2,a2)') ' decisive    (1/100) ', idec, 
     &			100.0*real(idec)/real(itot) ,'%' 
      write(fout,'(a22,i4,f8.2,a2)') ' total               ', itot, 
     &			100.0*real(itot)/real(nobs),'%' 
      write(fout,*)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c
c  write out bi - jointness /  a selection to fout and all to fout2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write(fout,*)
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*) ' x   Bi-Jointness measures    x'
      write(fout,*) ' x   bj > 75%  or bf > 3      x' 
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*)
      write(fout,*) '    i  j    pri     prj    prij     bj      bf'
      call barra('-',52,fout)
      l = 0
      do 804 i=1,kreg-1
        do 804 j=i+1,kreg 
           if ((bjoint1(i,j).ge.(75.0d0)).or.
     &        (bjoint2(i,j).ge.(3.0d0))) then
             l=l+1
             write(fout,'(3i3,4f8.2,e12.3)') l, i,j,
     &               beta(i),beta(j),
     &               bjoint(i,j),bjoint1(i,j),bjoint2(i,j)
          endif
 804    enddo
c
c  write out bi - dis-jointness 
c
!      write(fout,*)
!      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
!      write(fout,*) ' x   bi dis-jointness measures    x'
!      write(fout,*) ' x   bf < 1/100                   x'
!      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
!      write(fout,*)
!      write(fout,*) '    i  j    pri     prj    prij     rbj      bf'
!      call barra('-',52,fout)
!      l = 0     
!      do 806 i=1,kreg-1
!        do 806 j=i+1,kreg        
!          if    (bjoint2(i,j).le.(1/(1000.0)) ) then
!             l=l+1
!             write(fout,'(i5,2i3,4f8.2,e12.3)') l,i,j, 
!     &             beta(i),beta(j),
!     &             bjoint(i,j),bjoint1(i,j),bjoint2(i,j)
!           endif
! 806  enddo
      write(fout,*)
c.....................................................
c end of bi - jointness stuff
c.....................................................                         
c.....................................................
c tri - jointness stuff
c.....................................................    

      write(fouttj,*) ' i j k pri prj prk prijk tj bf '
      l=0
      lpos=0
      ifav=0
      ist=0
      ivst=0
      idec=0
      do 816 i=1,kreg
      do 816 j=i+1,kreg
      do 816 k=j+1,kreg
	      l = l+1
	      if (tjoint(i,j,k).gt.0) then
	        lpos=lpos+1
	        x1 = beta(i)+beta(j)+beta(k)-
     &               (bjoint(i,j)+bjoint(j,k)+bjoint(i,k))+
     &               tjoint(i,j,k)
	        if (x1.gt.0) then
	             tjoint1(i,j,k)=100.0d0*tjoint(i,j,k)/x1
	        else
	 	           write(fout,*) 'error in tri-joint'
	 	           fail=.true.
		           return
	       endif
	        x1 = x1 - tjoint(i,j,k)
	       if (x1.gt.0) then
	             tjoint2(i,j,k)=tjoint(i,j,k)/x1
	        else
	 	           write(fout,*) 'error in tri-joint'
	 	           fail=.true.
		           return
	        endif
	        else
	        tjoint1(i,j,k)=0.0d0
	        tjoint2(i,j,k)=0.0d0
	       endif
c for disjointness accumm
               tj2=tjoint2(i,j,k)	
	       if (tj2.lt.(1.d-2)) then
	           idec = idec + 1
	          else
	              if (tj2.lt.(1.0d0/30.0d0)) then
	                  ivst = ivst+1
	               else
	                 if (tj2.lt.(1.0d-1)) then
	                    ist = ist + 1
	                  else
	                    if (tj2.lt.(1.0d0/3.0d0)) ifav = ifav+1
	                 endif
	              endif
	        endif
c
c write tri-joint to fouttj
c    	         
	         write(fouttj,'(i5,3i3,6e13.6)') l,i,j,k,
     &                beta(i),beta(j),beta(k),
     &                tjoint(i,j,k),tjoint1(i,j,k),tjoint2(i,j,k)
 816  enddo
       write(fout,*)
       write(fout,'(a10,i5,a8,i5,a21,f6.2,a1)') 'there are ',lpos,
     &' out of ', l, 
     &' pijk greater than 0;',100.*real(lpos)/real(l),'%'
       write(fout,*) 
c
      write(fout,*)
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 
      write(fout,*) ' x   Tri-Jointness Measures   x'
      write(fout,*) ' x   tj1 > 70% or tj2 > 1     x'
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'  
      write(fout,*) 
      write(fout,820) 'i','j','k','pi','pj','pk','pijk','tj1','tj2'
 820  format(t8,a1,t11,a1,t14,a1,t19,a2,t27,a2,t35,a2,t42,a4,t51,a3,t60,
     &a3)
      call barra('-',65,fout)

      l=0   
      do 904 i=1,kreg-1
       do 904 j=i+1,kreg
        do 904 k=j+1,kreg
          if ((tjoint1(i,j,k).ge.(70.0d0)).or.
     &        (tjoint2(i,j,k).ge.(1.0d0))) then
             l=l+1
             write(fout,'(i6,3i3,5f8.2,e12.3)') l,i,j,k, 
     &             beta(i),beta(j),beta(k),
     &             tjoint(i,j,k),tjoint1(i,j,k),tjoint2(i,j,k)
          endif
 904    enddo

      itot = ifav + ist + ivst + idec
      write(fout,*)
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*) ' x  Tri- disjointness by bf ratio   x'
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*)
      write(fout,'(a22,i6,f8.2,a2)') ' favourable    (1/3) ', ifav, 
     &			100.0*real(ifav)/real(itot),'%'               
      write(fout,'(a22,i6,f8.2,a2)') ' strong       (1/10) ', ist, 
     &			100.0*real(ist)/real(itot),'%' 
      write(fout,'(a22,i6,f8.2,a2)') ' very strong  (1/30) ', ivst, 
     &			100.0*real(ivst)/real(itot),'%' 
      write(fout,'(a22,i6,f8.2,a2)') ' decisive    (1/100) ', idec, 
     &			100.0*real(idec)/real(itot) ,'%' 
      write(fout,'(a22,i6)') ' total               ', itot 
!     &			100.0*real(itot)/real(l),'%' 
      write(fout,*)

!      write(fout,*)
!      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx' 
!      write(fout,*) ' x   tri dis-jointness measures   x'
!      write(fout,*) ' x   tj < 1/1000   | pijk > 0     x'
!      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'  
!      write(fout,*) 
!      write(fout,820) 'i','j','k','pi','pj','pk','pijk','tj1','tj2'
!      call barra('-',65,fout)
!     
!      l=0   
!      do 906 i=1,kreg-1
!        do 906 j=i+1,kreg
!         do 906 k=j+1,kreg
!          if ((tjoint(i,j,k).gt.0).and.
!     &        (tjoint2(i,j,k).le.(1.0/1000.0))) then
!                   l=l+1
!             write(fout,'(i5,3i3,5f8.2,e12.3)') l,i,j,k, 
!     &             beta(i),beta(j),beta(k),
!     &             tjoint(i,j,k),tjoint1(i,j,k),tjoint2(i,j,k)
!          endif
! 906  enddo
      
      write(fout,*)
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*) ' x   Quad-jointness measure   x'
      write(fout,*) ' x   Look only at vars 1:4    x'
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*)
                  
      write(fout,'(a12, g12.6,a2)') 'pr[1:4]   =', qjoint,' %'
      write(fout,'(a12, g12.6,a2)') 'pr[union] =', qunion,' %'
      

      write(fout,'(a12, e12.6,a2)') 'QJ  =', 100.0d0*qjoint/qunion,' %'
      write(fout,'(a12, g12.6)')    'PO  =', qjoint/(qunion-qjoint)
      write(fout,*)
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*) ' x   Quint-jointness measure  x'
      write(fout,*) ' x   Look only at vars 1:5    x'
      write(fout,*) ' xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
      write(fout,*)
                  
      write(fout,'(a12, g12.6,a2)') 'pr[1:5]   =', q5joint,' %'
      write(fout,'(a12, g12.6,a2)') 'pr[union] =', q5union,' %'
      

      write(fout,'(a12, e12.6,a2)')'QJ  =',100.0d0*q5joint/q5union,' %'
      write(fout,'(a12, g12.6)')   'PO  =',q5joint/(q5union-q5joint)
      
c.....................................................
c end of jointness stuff
c.....................................................                         
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine computefj                                            ccc
ccc                                                                  ccc
ccc  modified:                                                       ccc
ccc       06.05.20 to accommodate non-uniform prior on model size    ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine computefj(iprior,model,nobs,kreg,ztz,z,y,ssqyn,
     &      kj,g0j,dstar,bstar,fj,fail)
c
c  for the choice of prior (1-9), see:
c   
c  c. fernandez, e. ley and mark f.j. steel (2001):
c "benchmark priors for bayesian model averaging" 
c  journal of econometrics, 100:2 (february), 381-427.
c
      implicit real*8(a-h,o-z)
      parameter (maxk=104,maxn=88)
c dummy args:      
      real*8 ztz(maxk,maxk),z(maxn,maxk),y(maxn),bstar(maxk)
      integer iprior
      logical fail
c local vars:      
      real*8 p(maxk),ztzj(maxk,maxk),zj(maxn,maxk),zy(maxk)
      integer model(maxk),ir(maxk)
c      
      fail = .false.

      dn  = dble(nobs)
      dk  = dble(kreg)     
      call ireg(kreg,model,ir,kj)
      dkj = dble(kj)
      
      if (kj.eq.0) then
         fj =  - 0.5d0*(dn-1.d0)*dlog(ssqyn)
         g0j = 1.d250
         goto 999
      endif
      
c
c  kj is the total number of regressors included (excl intercept)
c  ir() contains the list of regressors (excluding intercept)
c
c  select the kj regressors specified in ir()
c  we only need the upper part of ztzj
c
      do 20 j=1,kj
         do 10 i=1,j
 10      ztzj(i,j) = ztz(ir(i),ir(j))
         do 15 i=1,nobs
 15      zj(i,j) = z(i,ir(j))
 20   enddo

c since the z's have mean 0: 
c
c y'x (x'x)^{-1} x'y = ((\sum y)^2)/n + y'z bstar
c
c  compute y'z_j (z_j'z_j)^{-1} z_j'y by computing the cholesky
c  decomposition of (z'_{(j)} z'_{(j)}) and solving 
c  bstar = (z_j'z_j)^{-1} z_j'y         by backsubstitution
c
      call choldc(ztzj,kj,maxk,p,fail) 
      if (fail) return

      do 30 j=1,kj
         zy(j) = 0.0d0
         do 30 i=1,nobs
 30   zy(j) = zy(j) + zj(i,j)*y(i)
c
c       (x_j'x_j)^{-1} x'_j y
c       = (\bar y, (z_j'z_j)^{-1} z'_j y) =
c       = (\bar y, \beta^*_j)'
c 
      call cholsl(ztzj,kj,maxk,p,zy,bstar)
c
c ssqy = nobs * var(y) = y'y - ((y'1)**2)/nobs
c      
      ymzjy = ssqyn
      do 40 j=1,kj
 40   ymzjy = ymzjy - zy(j)*bstar(j)
c
c different priors and stuff...
c
      
      if (iprior.eq.1) then
         g0j = 1.0d0/dn
       elseif (iprior.eq.2) then
         g0j = dkj/dn
       elseif (iprior.eq.3) then
         g0j = (dk**(1.0d0/dkj))/dn
       elseif (iprior.eq.4) then
          g0j = dsqrt(dkj/dn)	 
       elseif (iprior.eq.5) then
         g0j = (dlog(dn))**(-3)
       elseif (iprior.eq.6) then
         g0j = dlog(dkj+1.0d0)/dlog(dn)
       elseif (iprior.eq.7) then
         aux = 0.15411d0*0.64889d0**(1.0d0/dkj)
         g0j = aux/(1.0d0 - aux)
       elseif (iprior.eq.8) then
         g0j = 1.0d0/dsqrt(dn)
       elseif (iprior.eq.9) then
         g0j = 1.0d0/(dk**2)
      endif
      
      aux1 = 1.0d0/(1.0d0 + g0j)
      aux2 = 1.0d0 - aux1
      aux3 = 0.5d0*dkj*dlog(aux2)
         
      dstar = aux1*ymzjy + aux2*ssqyn
c
c finally log(f(j))
c
      fj = aux3 - 0.5d0*(dn-1.d0)*dlog(dstar)
c
c for handling informative priors on model size
c  
 999  continue
                 
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc  subroutine runchain                                             ccc
ccc                                                                  ccc
ccc  modified:                                                       ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine runchain(fix_start,iprior,idum,initrep,mnumrep,nobs,
     &      kreg,ztz,z,y,ssqyn,
     &      visits,fjlog,g0j,dstar,bstar,
     &      midx,kj,
     &      imax,nvout,fout,fjout,fjsum,fail)
c
      implicit real*8(a-h,o-z)
      integer fout
      parameter(maxm=200000,maxm2=5000,maxk=104,maxn=88)
c  dummy arguments:
      real*8 z(maxn,maxk),y(maxn),ztz(maxk,maxk),bstar(maxk,maxm)
      logical fail, fix_start 
      integer kj(maxm)
c local vars:      
      real*8 visits(maxm),fjlog(maxm),g0j(maxm),dstar(maxm),
     &      bstarfj(maxk),bstarfjold(maxk),
     &      fjstack(maxm2),g0jstack(maxm2),dstarstack(maxm2),
     &      bstarstack(maxk,maxm2)
      logical move, visited, evaluated    
      real*8 midx(maxm,2),idxstack(maxm2,2),idxold(2),idxnew(2)
      integer model(maxk),mc(maxk),kjstack(maxm2),model_start(maxk)
      
c................................................................
      
      fail = .false.
      fjout = 0.0d0
      nvout = 0
      model_start = 0
      model_start(1:10) = 1
c     
c-----------------------------------------------------------------------
c  run mc3
c-----------------------------------------------------------------------
c 1. start the chain
c-----------------------------------------------------------------------
c
c  generate one model at random
c 
c **k67**
      
      
      if (kreg.gt.52) then
         if (fix_start) then
               call get2modidx(kreg, model_start,idxold(1),idxold(2)) 
         else  
         idxold(1) =  (ran2(idum)*(2.0d0**52) + 1.0d0)
         k2 = kreg - 52
         idxold(2) =  (ran2(idum)*(2.0d0**k2) + 1.0d0)
         end if
      else
         if (fix_start) then
             idxold(1) = getmodidx(kreg, model_start)
         else    
         idxold(1) =  (ran2(idum)*(2.0d0**kreg) + 1.0d0)
         idxold(2) = 0
         end if
      endif
      call g2model(idxold(1),idxold(2),kreg,model)
      write(fout,*) 
      write(fout,*)' Starting model is:'
      write(fout,*) ' idxold(1): ', idxold(1)
      if (kreg.gt.52) write(fout,*) ' idxold(2): ', idxold(2)
c **k67**
      write(fout,*) ' model(1:kreg): '
      write(fout,*) (model(i),i=1,kreg)
      write(fout,*)
      
      call computefj(iprior,model,nobs,kreg,ztz,z,y,ssqyn,
     &      kjfjold,g0jfjold,dstarfjold,bstarfjold,fjold,fail)
      if (fail) then
      write(fout,*) 'start of chain: fail in computefj!'
      return
      endif
c
c start storing from the bottom of the stack
c     
c **k67**
      idxstack(maxm2,1) = idxold(1)
      idxstack(maxm2,2) = idxold(2)
c **k67**
      fjstack(maxm2) = fjold
      kjstack(maxm2) = kjfjold
      g0jstack(maxm2) = g0jfjold
      dstarstack(maxm2) = dstarfjold
      do 260 j=1,kjstack(maxm2)
 260  bstarstack(j,maxm2) = bstarfjold(j)
c
c  initialize the chain ---run for initrep repetitions---
c 
      write(*,*) '... initializing chain ... '
c      
c  flush all values in stack!
c		  
      nmv2 = 1
      iwrite=initrep/5

      do 8090 init=1,initrep

      if (mod(init,iwrite).eq.0) then 
       write(*,'("Warmup Draw:",i7,":",i7,x,f4.0,"%")') 
     & init,initrep,
     & 100.0d0*(dble(init)/dble(initrep+mnumrep))
      endif
                        
         call mcand(kreg,model,mc,idum)
c     
c   get model index
c 
         evaluated = .false.
c **k67**         
         call get2modidx(kreg,mc,idxnew(1),idxnew(2))
         if ((idxnew(1).eq.idxold(1)).and.
     &       (idxnew(2).eq.idxold(2))) goto 8090
c **k67**
c
c  see if it is stored in the stack
c
          imin = maxm2-(nmv2-1)
          do 7010 i = imin, maxm2
c **k67**
               if ((idxnew(1).eq.idxstack(i,1)).and. 
     &             (idxnew(2).eq.idxstack(i,2))) then
c **k67**     
                  evaluated = .true.
                  fjnew = fjstack(i)
                  index = i
                  goto 7299
               endif
 7010     enddo
 
 7299    continue
c
c  if it was not stored, then compute fj and other stuff
c
         if (.not. evaluated) then
	   
            call computefj(iprior,mc,nobs,kreg,ztz,z,y,ssqyn,
     &            kjfj,g0jfj,dstarfj,bstarfj,fjnew,fail)
     
      if (fail) then
      write(fout,*) init,'warm-up: fail in computefj!'
      return
      endif
c     
c  store it for future use.  
c  start storing from the bottom of stack to the top,
c  when full, flush the bottom
c
            if (nmv2.lt.maxm2) then
               nmv2 = nmv2+1
               index = maxm2-(nmv2-1)
            else
               do 7350 i=maxm2,2,-1
c **k67**
                  idxstack(i,1) = idxstack(i-1,1)
                  idxstack(i,2) = idxstack(i-1,2)
c **k67**
                  fjstack(i) = fjstack(i-1)
                  kjstack(i) = kjstack(i-1)
                  g0jstack(i) = g0jstack(i-1)
                  dstarstack(i) = dstarstack(i-1)
                  do 7340 j=1,kjstack(i-1)
 7340             bstarstack(j,i) = bstarstack(j,i-1)
 7350          continue
               index = 1
            endif
c **k67**	    
            idxstack(index,1) = idxnew(1)
            idxstack(index,2) = idxnew(2)
c **k67**            
            fjstack(index) = fjnew
            kjstack(index) = kjfj
            g0jstack(index) = g0jfj
            dstarstack(index) = dstarfj
            do 7360 j=1,kjstack(index)
 7360       bstarstack(j,index) = bstarfj(j)
         endif
c     
c end of stuff when evaluated was false
c


c
c see if chain moves to new model
c
         dif = (fjnew-fjold) 
	 
         if (dif .ge. 0.d0) then
            move = .true.
         elseif (dif.ge.dlog(ran2(idum)))  then
            move = .true.
         else
            move = .false.
         endif
	 
         if (move) then
	    do 8000 i=1,kreg
 8000       model(i) = mc(i)
            fjold = fjnew
c **k67**
            idxold(1) = idxnew(1)
            idxold(2) = idxnew(2)
c **k67**            
            indexstack = index
         endif
	 
 8090 continue              
c
c here is the first model visited:
c
      nmv  = 1
      imd = 1

      visits(imd) = 1.0d0
c **k67**
      midx(imd,1) = idxold(1)
      midx(imd,2) = idxold(2)
c **k67**
      fjlog(imd) = fjold
      kj(imd) = kjstack(indexstack)
      g0j(imd) = g0jstack(indexstack)
      dstar(imd) = dstarstack(indexstack)
      do 8092 j=1,kj(imd)
 8092 bstar(j,imd) = bstarstack(j,indexstack)
      
      do 8095 i=2,maxm
 8095 visits(i) = 0.0d0

c
c-----------------------------------------------------------------------
c  2. run the chain
c-----------------------------------------------------------------------
c             
      call wr_date(fout,.false.,.true.)
      write(*,*) '... running chain ... '

      numrep = 1
      fjlogmean = fjold
      
      iwrite=mnumrep/10
      
 8100 numrep = numrep + 1
c
c  keep the terminal informed of each 1/10 draws
c
       if (mod(numrep,iwrite).eq.0) then
        write(*,'("Chain Draw:",i7,":",i7,x,f4.0,"%")') 
     & numrep,mnumrep,
     & 100.0d0*(dble(numrep+initrep)/dble(initrep+mnumrep))
       endif
c     
c  pick candidate
c 
      call mcand(kreg,model,mc,idum)
c     
c  get model index
c 
c **k67**
      call get2modidx(kreg,mc,idxnew(1),idxnew(2))
      if ( (idxnew(1).eq.idxold(1)).and.
     &     (idxnew(2).eq.idxold(2)) ) goto 8500
               
      visited   = .false.
      evaluated = .false.
c **k67**
c     
c see if it is an already _visited_ model 
c                  
      do 8200 i=1, min(nmv,maxm)
c **k67**
         if ( (idxnew(1).eq.midx(i,1)).and.
     &        (idxnew(2).eq.midx(i,2)) ) then
c **k67**
            visited = .true.
            evaluated = .true.
            fjnew = fjlog(i)
            index = i
            goto 8299
         endif
 8200 enddo
c     
c otherwise see if it is an already _evaluated_ model
c
      if (nmv2.ge.1) then
         do 8210 i = maxm2-(nmv2-1), maxm2
c **k67**
            if ( (idxnew(1).eq.idxstack(i,1)).and.
     &           (idxnew(2).eq.idxstack(i,2)) ) then 
c **k67**
               evaluated =.true.
               fjnew = fjstack(i)
               index = i
               goto 8299
            endif
 8210    enddo
      endif
            
 8299 continue
c     
c  if we don't have logfj stored, then compute it
c
      if (.not. evaluated) then
         call computefj(iprior,mc,nobs,kreg,ztz,z,y,ssqyn,
     &         kjfj,g0jfj,dstarfj,bstarfj,fjnew,fail)
         if (fail) then
            write(fout,*) numrep,'chain: fail in computefj!'
            return
         endif
      endif

c
c see if chain moves to new model
c
         dif = (fjnew-fjold) 
      
         if (dif .ge. 0.0d0) then
            move = .true.
         elseif (dif.ge.dlog(ran2(idum)))  then
            move = .true.
         else
            move = .false.
         endif
c
c
c
         if (move) then
            
            do 8300 i=1,kreg
 8300       model(i) = mc(i)
            fjold  = fjnew
c **k67**
            idxold(1) = idxnew(1)
            idxold(2) = idxnew(2)
c **k67**
            if (visited) then
               imd = index
               goto 8500
            else
               nmv = nmv+1
               if (nmv.le.maxm) then
                  imd = nmv
c **k67**
                  midx(imd,1) = idxnew(1)
                  midx(imd,2) = idxnew(2)
c **k67**
                  fjlog(imd) = fjnew
                  fjlogmean = fjlogmean + fjnew
               else
                  nvout = nvout + 1
                  fjout =fjout + dexp(fjnew)
                  goto 8501
               endif
            endif
c
c  we get to the next if when we wanna move _and_ there is 'room'
c  to move
c
         if (evaluated) then
c
c   move stuff from "evaluated" stack to "visited" stack
c      
            kj(imd) = kjstack(index)
            g0j(imd) = g0jstack(index)
            dstar(imd) = dstarstack(index)
            do 8310 j=1,kj(imd)
 8310       bstar(j,imd) = bstarstack(j,index)
         else
c
c   if (.not. evaluated) computefj was called above and now
c   we can store stuff in "visited" stack
c 
            kj(imd) = kjfj
            g0j(imd) = g0jfj
            dstar(imd) = dstarfj
            do 8312 j=1,kj(imd)
 8312       bstar(j,imd) = bstarfj(j)
         endif
c
c  don't move and model not evaluated before
c            
      elseif (.not. evaluated) then
c
c  store in "evaluated" stack
c 
         if (nmv2.lt.maxm2) then
            nmv2 = nmv2+1
            index = maxm2-(nmv2-1)
         else
c     
c when stack is full, disregard value at bottom, move everyone
c 1 position downwards and store new value at top
c
            do 8350 i=maxm2,2,-1
c **k67**
               idxstack(i,1) = idxstack(i-1,1)
               idxstack(i,2) = idxstack(i-1,2)
c **k67**
               fjstack(i) = fjstack(i-1)
               kjstack(i) = kjstack(i-1)
               g0jstack(i) = g0jstack(i-1)
               dstarstack(i) = dstarstack(i-1)
               do 8340 j=1,kjstack(i-1)
 8340          bstarstack(j,i) = bstarstack(j,i-1)
 8350       continue
            index = 1
         endif
c **k67**
         idxstack(index,1) = idxnew(1)
         idxstack(index,2) = idxnew(2)
c **k67**         
         fjstack(index) = fjnew
         kjstack(index) = kjfj
         g0jstack(index) = g0jfj
         dstarstack(index) = dstarfj
         do 8360 j=1,kjstack(index)
 8360    bstarstack(j,index) = bstarfj(j)
      endif
c     
 8500 visits(imd) = visits(imd) + 1.0d0
c
c we get directly to 8501 when the visited-models stack was full
c
 8501 continue
      if (numrep.lt.mnumrep) goto 8100
c
c-----------------------------------------------------------------------
c  end of chain loop
c-----------------------------------------------------------------------
c                                 
      imax = min(maxm,nmv)
      fjlogmean = fjlogmean/dble(imax)
      write(fout,*)
      write(fout,*) 'fjlogmean...', fjlogmean
      write(fout,*)
c
c-----------------------------------------------------------------------
c  let's do the accounting...
c-----------------------------------------------------------------------
c      
      fjsum = 0.0d0
c
c  when we visit models for which we have no room in the
c  visited-models stack we still want [sum fjlog(i)] = 1;
c  we thus comment out the following line
c
c      if (nvout.gt.0) fjsum = fjout/dexp(fjlogmean)
c
c fjout = sum dexp() already
c
      fjout = fjout/exp(fjlogmean)

      do 9000 i=1,imax
         fjlog(i) = dexp(fjlog(i) - fjlogmean)
         fjsum = fjsum + fjlog(i)
         visits(i) = visits(i)/dble(mnumrep-nvout)
 9000 enddo
      
      fjsum2=0.0d0
      fjsum3=0.0d0
      do 9010 i=1,imax
        fjlog(i) = fjlog(i)/fjsum
        fjsum2 = fjsum2 + fjlog(i)
        fjsum3 = fjsum3 + visits(i)
 9010 enddo 
 
c      write(*,*)
c      write(*,   '("Sum of Normalised BFs:   ",e20.10)') fjsum2
c      write(*,   '("Sum of Normalised Visits:",e20.10)') fjsum3
c      write(fout,*)
c      write(fout,'("Sum of Normalised BFs:   ",e20.10)') fjsum2
c      write(fout,'("Sum of Normalised Visits:",e20.10)') fjsum3
c      write(fout,*)
      
      fjout = fjout/(fjsum+fjout)
      
      if (nvout .gt. 0) write(fout,9020) maxm,fjout*100.0d0
 9020 format(' mass not in the stored ',i7,' models is ',g6.4,' %')
c
c  compute corr coefficient
c 
      sumb = 0.0d0
      sumf = 0.0d0
      ssqb = 0.0d0
      ssqf = 0.0d0
      corr = 0.0d0
      
      do 9100 i=1,imax
         sumb = sumb +  fjlog(i)
         ssqb = ssqb +  fjlog(i)**2
         sumf = sumf +  visits(i)
         ssqf = ssqf +  visits(i)**2
 9100 enddo
 
      sumb = sumb/dble(imax)
      sumf = sumf/dble(imax)
      ssqb = sqrt(ssqb/dble(imax) - sumb**2)
      ssqf = sqrt(ssqf/dble(imax) - sumf**2)
      
      do 9200 i=1,imax
 9200 corr = corr+(fjlog(i)-sumb)*(visits(i)-sumf)
      corr = corr/(ssqb*ssqf*dble(imax))
      
      write(fout,*)
      write(fout,9209) sumb,sumf
 9209 format(2x,'means of fj and freq are    ',2e15.4)
      write(fout,9210) ssqb,ssqf
 9210 format(2x,' stds of fj and freq are    ',2e15.4)
      write(fout,9215) corr
 9215 format(2x,' corr coef of fj w/ freq is ',e15.4)
      write(fout,*)

      fjsum = fjsum*dexp(fjlogmean)
     
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc                                                                  ccc
ccc   function: gammln                                               ccc
ccc                                                                  ccc
ccc this function computes ln(gamma(x)) for x > 0.                   ccc
ccc taken from pascal numerical recipes p 177                        ccc
ccc                                                                  ccc
ccc modified to improve accuracy when x<1 using:                     ccc
ccc                                                                  ccc
ccc                     pi*e                                         ccc
ccc gamma(1-e) = ---------------------                               ccc
ccc               gamma(1+e)*sin(pi*e)                               ccc
ccc                                                                  ccc
ccc                                                                  ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function gammln(xx)
c
      real*8 xx
      real*8 cof(6),stp,half,one,fpf,x,y,tmp,ser,g,dpi
c
      data cof/76.18009173d0,-86.50532033d0,24.01409822d0, 
     &         -1.231739516d0,0.120858003d-2,-.536382d-5/
      data stp/2.50662827465d0/
      data dpi/3.14159265358979d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
c                 
      if (xx.lt.1) then
         x = one - xx
         tmp =  x + fpf
         tmp = (x + half)*dlog(tmp)-tmp
         ser = one
         do 1000 j=1,6
            x=x+one
            ser = ser+cof(j)/x
 1000    continue
         g = tmp + dlog(stp*ser)
         y = dpi*(one-xx)
         gammln = dlog(y/sin(y)) - g
c
       else
c
         x = xx - one
         tmp =  x + fpf
         tmp = (x + half)*dlog(tmp)-tmp
         ser = one
         do 1100 j=1,6
            x=x+one
            ser = ser+cof(j)/x
 1100    continue
         gammln = tmp + dlog(stp*ser)
      endif
c
      return
      end
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  basic mc3m stuff
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c**********************************************
      subroutine ireg(n,m,ir,kj)
c**********************************************
c      
c inputs:
c     n.......is the number of all possible regressors
c     m(n)...m(i) contains 0s or 1s in each coordinate
c output:
c     ir.....is a kj-dim vector with the indeces of the included regressors
c
      integer  n,m(n),ir(n)
c
      ii=1
      do 10 i=1,n
         if (m(i).eq.1) then
            ir(ii)=i
            ii=ii+1
         endif
 10   continue
      kj = ii -1
      return
      end
c
c**********************************************
      real*8 function getmodidx(k,m)
c**********************************************
c      
c inputs:
c     k.le.52.......is the number of all possible regressors
c     m(k)....contains 0s or 1s in each coordinate
c output:
c     getmodidx..is the real*8 associated with the state vector
c
      real*8 idx,z
      integer  k,i,m(k)
      
        z = 2.0d0**(k-1)
        idx = 1      
        do 10 i=1,k
         idx = idx + z*m(i)
         z = z/2.0d0
 10     continue
        getmodidx = idx
        return
 
      end
c
c
c**********************************************
      subroutine get2modidx(k,m,idx1,idx2)
c**********************************************
c      
c inputs:
c     k>52.......is the number of all possible regressors
c     m(1:k)....contains 0s or 1s in each coordinate
c output:
c     idx1..is the real*8 associated with the 1st 52-var state vector
c     idx2..is the real*8 associated with the 2nd 52-var state vector
c   
      integer  i,k, m(k), m1(52), m2(52) 
      real*8 idx1,idx2,getmodidx
      
      k1=min(52,k)     
      do 10 i=1,k1
         m1(i)=m(i)
 10   continue
      idx1 = getmodidx(k1,m1)
      
      k2 = k-52
      if (k2.ge.1) then
         do 15 i=1,k2
            m2(i)=m(i+52)
 15      continue
         idx2 = getmodidx(k2,m2)
      else
         idx2=0.0d0
      endif
      
      return
      end
c
c**********************************************
      subroutine gmodel(idx,k,m)
c**********************************************
c
c  given the model index idx and the number of all possible regressors 
c  k<53
c  it returns a binary array, m, of dimension k.  if the ith coordinate is 1
c  the ith regressor is included in the model.  otherwise it is excluded.
c      
c inputs:
c     idx.....is the model index (real*8)
c     k.......is the number of all possible regressors k=1..52
c
c output:
c     m(k)...contains 0s or 1s in each coordinate
c
      integer i, k, m(k)
      real*8 idx,x,z

      x = 2.0d0**(k-1)
      z = idx
 
      do 20 i=1,k
         if (z.gt.x) then
            m(i) = 1
            z = z - x
           else
            m(i) = 0
         endif
         x = x / 2.0d0
 20   continue
      return
      end
c
c
c**********************************************
      subroutine g2model(idx1,idx2,k,m)
c**********************************************
c
c  given the model index idx and the number of all possible regressors k
c  it returns a binary array, m, of dimension k.  if the ith coordinate is 1
c  the ith regressor is included in the model.  otherwise it is excluded.
c      
c inputs:
c     idx.....is the model index (real*8)
c     k>52....is the number of all possible regressors
c
c output:
c     m(k)...contains 0s or 1s in each coordinate
c
      integer  i,k,m1(52),m2(52),m(k)
      real*8 idx1,idx2
      
      k1 = min(k,52)
      call gmodel(idx1,k1,m1)
      do 10 i=1,k1
         m(i)=m1(i)
 10   continue
      
      k2 = k-52
      if (k2.ge.1) then
         call gmodel(idx2,k2,m2)
         do 15 i=1,k2
            m(52+i)=m2(i)
 15      continue
      endif
      
      return
      end
c
c**********************************************
      subroutine g2model2(idx1,idx2,k,m,kj,ir)
c**********************************************
c
c  given the model index idx and the number of all possible regressors k
c  it returns a binary array, m, of dimension k.  if the ith coordinate is 1
c  the ith regressor is included in the model.  otherwise it is excluded.
c      
c inputs:
c     idx.....is the model index (real*8)
c     k>52....is the number of all possible regressors
c
c outputs:
c     m(k)...contains 0s or 1s in each coordinate
c     kj is the number of included regressors
c     if() has in positions 1..kj the indexes of included regressors      
c
      integer k,m(k),ir(k),kj
      real*8 idx1,idx2
      
      call g2model(idx1,idx2,k,m)
      call ireg(k,m,ir,kj)
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                               c
c   c**********************************************                             c
c         subroutine gmodel2(j,n,m,kj,ir)                                       c
c   c**********************************************                             c
c   c                                                                           c
c   c  given the model index j (and the number of all possible regressors n)    c
c   c  it returns a binary array of dimension n.  if the kth coordinate is 1    c
c   c  the kth regressor is included in the model.  otherwise it is excluded.   c
c   c                                                                           c
c   c inputs:                                                                   c
c   c     j.......is the model index                                            c
c   c     n.......is the number of all possible regressors                      c
c   c                                                                           c
c   c output:                                                                   c
c   c     m(n)...m(i) contains 0s or 1s in each coordinate                      c
c   c     kj.....number of included regressors                                  c
c   c     ir(n)..array w/ indexes of incl regressors                            c
c   c                                                                           c
c         integer  n,m(n),ir(n)                                                 c
c         real*8 j,ni, k                                                        c
c                                                                               c
c         ni = 2.0d0**(n-1)                                                     c
c         k = j                                                                 c
c         kj = 0                                                                c
c         ii = 1                                                                c
c                                                                               c
c         do 20 i=1,n                                                           c
c            if (k.gt.ni) then                                                  c
c               m(i) = 1                                                        c
c               kj = kj + 1                                                     c
c               k = k-ni                                                        c
c               ir(ii) = i                                                      c
c               ii = ii + 1                                                     c
c            else                                                               c
c               m(i) = 0                                                        c
c            endif                                                              c
c            ni = ni/2.0d0                                                      c
c    20   continue                                                              c
c         return                                                                c
c         end                                                                   c
c                                                                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                       c
c   c**********************************************     c
c         integer function icard(n,m)                   c
c   c**********************************************     c
c   c                                                   c
c   c  returns the number of 1's in binary array m(n)   c
c   c                                                   c
c         integer i,j,n,m(n)                            c
c         j=0                                           c
c         do 10 i=1,n                                   c
c           j=j+m(i)                                    c
c    10   continue                                      c
c         icard=j                                       c
c         return                                        c
c         end                                           c
c                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c**********************************************
      subroutine mcand(n,m,mc,idum)
c**********************************************
c
c  given the model represented by the binary array m(n)
c  randomly choose candidate model, mc, in nbd(m)
c
      implicit real*8 (a-h,o-z)
      integer n, m(n), mc(n), idum
	        
      do 100 i=1,n
 100  mc(i) = m(i)
c
c
c |--|----------|
c 0  1          n
c     
c draw from 0,1,...,n
c
      i1 = int(ran2(idum)*(n+1.0d0))
c
c if we get 0 we stay w/ same model
c otherwise we just switch the value of the i1th coordinate
c
      if (i1.gt.0) mc(i1) = abs(mc(i1)-1)
	  return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  misc aux routines
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------------------------------------------------
c
      subroutine barra(car,n,fout)
c
c-----------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
c
c
      integer i, n, fout
      character*1 car
      character*1 bar(85)
      n2=min(85,n)
      do 1000 i=1,n2
 1000 bar(i) = car
      write(fout,*) (bar(i),i=1,n2)
      return
      end
c-----------------------------------------------------------------------
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc  random 
ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c
c-----------------------------------------------------------------------
c
      real*8 function ran2(idum)
c
c   this function contains a portable random number generator 
c   uniform [0,1]
c   (numerical recipes p. 273)
c
      implicit real*8 (a-h,o-z)
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*8 am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1.d0/im1,imm1=im1-1,
     &      ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     &      ntab=32,ndiv=1+imm1/ntab,eps=1.2d-7,rnmx=1.0d0-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do 11 j=ntab+8,1,-1
            k=idum/iq1
            idum=ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum=idum+im1
            if (j.le.ntab) iv(j)=idum
 11      continue
         iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=dmin1(am*iy,rnmx)
      return
      end
c
c-----------------------------------------------------------------------

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc
ccc  matrix inversion, determinant, etc routines
ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine choldc(a,n,np,p,fail)
c
c a(1:n,1:n); n =< np (physical dimension)
c      
      logical fail
      integer n,np
      real*8 a(np,np),p(n)
      integer i,j,k
      real*8 sum
      fail = .false.
      do 13 i=1,n
         do 12 j=i,n
           sum=a(i,j)
          do 11 k=i-1,1,-1
             sum=sum-a(i,k)*a(j,k)
 11       continue
          if(i.eq.j)then
             if(sum.le.0.d0) then
                fail = .true.
                return
             endif
             p(i)=sqrt(sum)
         else
             a(j,i)=sum/p(i)
         endif
 12     continue
 13   continue
      return
      end

      subroutine cholsl(a,n,np,p,b,x)
      integer n,np
      real*8 a(np,np),b(n),p(n),x(n)
      integer i,k
      real*8 sum
      do 12 i=1,n
        sum=b(i)
        do 11 k=i-1,1,-1
          sum=sum-a(i,k)*x(k)
11      continue
        x(i)=sum/p(i)
12    continue
      do 14 i=n,1,-1
        sum=x(i)
        do 13 k=i+1,n
          sum=sum-a(k,i)*x(k)
13      continue
        x(i)=sum/p(i)
14    continue
      return
      end
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine indexx(n,arr,indx)
      integer n,indx(n),m,nstack
      real*8 arr(n)
      parameter (m=7,nstack=50)
      integer i,indxt,ir,itemp,j,jstack,k,l,istack(nstack)
      real*8 a
      do 11 j=1,n
         indx(j)=j
 11   continue
      jstack=0
      l=1
      ir=n
 1    if(ir-l.lt.m)then
         do 13 j=l+1,ir
            indxt=indx(j)
            a=arr(indxt)
            do 12 i=j-1,1,-1
               if(arr(indx(i)).le.a)goto 2
               indx(i+1)=indx(i)
 12         enddo
            i=0
 2          indx(i+1)=indxt
 13      enddo
         if(jstack.eq.0)return
         ir=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+ir)/2
         itemp=indx(k)
         indx(k)=indx(l+1)
         indx(l+1)=itemp
         if(arr(indx(l+1)).gt.arr(indx(ir)))then
            itemp=indx(l+1)
            indx(l+1)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l)).gt.arr(indx(ir)))then
            itemp=indx(l)
            indx(l)=indx(ir)
            indx(ir)=itemp
         endif
         if(arr(indx(l+1)).gt.arr(indx(l)))then
            itemp=indx(l+1)
            indx(l+1)=indx(l)
            indx(l)=itemp
         endif
         i=l+1
         j=ir
         indxt=indx(l)
         a=arr(indxt)
 3       continue
         i=i+1
         if(arr(indx(i)).lt.a)goto 3
 4       continue
         j=j-1
         if(arr(indx(j)).gt.a)goto 4
         if(j.lt.i)goto 5
         itemp=indx(i)
         indx(i)=indx(j)
         indx(j)=itemp
         goto 3
 5       indx(l)=indx(j)
         indx(j)=indxt
         jstack=jstack+2
         if(jstack.gt.nstack)pause 'nstack too small in indexx'
         if(ir-i+1.ge.j-l)then
            istack(jstack)=ir
            istack(jstack-1)=i
            ir=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         endif
      endif
      goto 1
      end
c
c-----------------------------------------------------------------------
c
      subroutine wr_date(filename,wfile,wscreen)
c
c-----------------------------------------------------------------------
c
      integer filename
      character*24 the_date
      logical wfile,wscreen

      call fdate(the_date)
      
      if (wfile) then
         write(filename,101)
         write(filename,100) the_date
         write(filename,101)
       endif
         
      if (wscreen) then
          write(*,101)
          write(*,100) the_date
          write(*,101)
      endif
          
  100 format('|',4x,a24,4x,'|')
  101 format('*--------------------------------*')
      return
      end                                                                      
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      real*8 function studt1(zst,nust,bst,ast)
      
      implicit real*8 (a-h,o-z)
      real*8 zst,bst,ast
      integer nust
      
      studt1 = 
     & (dexp(gammln(0.5d0*(1.d0+nust))-gammln(0.5d0*nust))
     & /dsqrt(3.14159265358979d0*nust*ast))/
     & (1.d0+((zst-bst)**2)/(nust*ast))**(0.5d0*(1.d0+nust))
      
      return
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine lps(iprior,imax,ibm,nobs,kreg,
     &      visits,g0j,dstar,bstar,midx, 
     &      ymean,ssqyn,ztz,z,y,nf,yf,zf,fout,
     &      avelps,bestlps,fulllps,nulllps,
     &      fail)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      
      implicit real*8(a-h,o-z)
      integer fout
      parameter(maxm=200000,maxk=104,maxn=88,maxnf=50)
c dummy args:      
      real*8 y(maxn),yf(maxnf),
     & zf(maxk,maxnf),z(maxn,maxk),ztz(maxk,maxk),
     & bstar(maxk,maxm),dstar(maxm),g0j(maxm),visits(maxm),
     & midx(maxm,2)
c local vars:  
      real*8 p(maxk),zfj(maxk,maxnf),
     & zfdummy(maxk),yvarfull(maxnf),ymfull(maxnf),ztzj(maxk,maxk),
     & ym(maxnf,maxm),yvar(maxnf,maxm),bstarfull(maxk),
     & fulllps,nulllps
      integer ir(maxk), model(maxk) 
      logical fail
    
      fail = .false.           
      nu = nobs-1
      yvar2=ssqyn/dble(nu)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  full model 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c full model, compute ymfull and yvarfull
c
      do 100 j=1,kreg
         model(j) = 1
 100  enddo
      
      call computefj(iprior,model,nobs,kreg,ztz,z,y,ssqyn,
     &      kjfull,gjfull,dstarfull,bstarfull,fjfull,fail)
      if (fail) then
         write(fout,*) 'fail in lps()...full...computefj'
         return
      endif
      
      call choldc(ztz,kreg,maxk,p,fail)
      if (fail) then
         write(fout,*) 'fail in lps()...full...choldc'
         return
      endif
      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                        c
c   loop over observations to be predicted: i_f=1..nf    c
c                                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      do 150 i_f=1,nf
         call cholsl(ztz,kreg,maxk,p,zf(1,i_f),zfdummy)
         zzz = 0.0d0
         do 125 j=1,kreg
            zzz = zzz + zf(j,i_f)*zfdummy(j)
 125     enddo
         zzz = 1.d0+1.d0/nobs+ zzz/(1.d0+gjfull)
         yvarfull(i_f) = zzz*dstarfull/(nobs-1.d0)
         
         yyy = 0.0d0
         do 140 j=1,kreg
            yyy = yyy + zf(j,i_f)*bstarfull(j)
 140     enddo
         ymfull(i_f) = ymean + yyy/(1.d0+gjfull)
 150  enddo

cccccccccccccccccccccccccccccccc
c                              c
c   Now for all other models   c
c                              c
cccccccccccccccccccccccccccccccc
c     loop over models         c
cccccccccccccccccccccccccccccccc
      
      write(*,*) '... looping over models ...',imax
      
      do 1000 imd = 1,imax      
c
c first, generate model, get kj and ir()
c            
         call g2model2(midx(imd,1),midx(imd,2),kreg,model,kj,ir)
c
c  select the kj regressors specified in ir() from z'z and z_f
c  also substract sample averages of z_j's to the z_{f,j}'s.
c
         do 20 j=1,kj
            do 15 i=1,kj
               ztzj(i,j) = ztz(ir(i),ir(j))
 15         enddo
            do 17 i=1,nf
c
c  zmean must be substracted from zf() if not done above
c
               zfj(j,i)  = zf(ir(j),i)
 17         enddo
 20      enddo
c
c  compute cholesky dcmp of (Z_j'Z_j)
c 
         call choldc(ztzj,kj,maxk,p,fail)
         if (fail) then
            write(fout,*) (ir(j),j=1,kj)
            do 21 i=1,kj
               write(fout,22) (ztzj(i,j),j=1,i)
 21         enddo
 22         format(29(:f10.2))
            write(16,*) '...fail choldc @ predict...'
            return
         endif
c      
c  loop over columns of z'_{f,j} w/ cholsl to obtain
c  (z'_j z_j)^{-1} z'_{f,j}
c
c
c loop over observations to be predicted: i_f=1..nf
c
         do 200 i_f=1,nf
            call cholsl(ztzj,kj,maxk,p,zfj(1,i_f),zfdummy)
            zzz = 0.0d0
c
c  premultiply by  z_{f,j} to get
c  zzz <-- z_{f,j} (z'_j z_j)^{-1} z'_{f,j}
c
            do 25 j=1,kj
               zzz = zzz + zfj(j,i_f)*zfdummy(j)
 25         enddo
            zzz = 1.d0+1.d0/nobs+ zzz/(1.d0+g0j(imd))
            
            yvar(i_f,imd) = zzz*dstar(imd)/dble(nu)
            
            yyy = 0.0d0
            do 40 j=1,kj
               yyy = yyy + zfj(j,i_f)*bstar(j,imd)
  40        enddo
            ym(i_f,imd) = ymean + yyy/(1.d0+g0j(imd))
            
            
 200     enddo
 1000 continue
c     
c loop over yf's for bma and best model's predictive densities
c
      write(*,*) '... looping over nf: ',nf
      
      avelps = 0.0d0
      bestlps = 0.0d0
      fulllps = 0.0d0
      nulllps=0.0d0
      
      do 3400 i_f = 1,nf
         
         xlps = 0.0d0     
         do 1200 imd=1,imax      
            xlps = xlps + visits(imd)*
     &            studt1(yf(i_f),nu,ym(i_f,imd),yvar(i_f,imd))            
 1200    enddo
         
         bestlps = bestlps - 
     &         dlog(studt1(yf(i_f),nu,ym(i_f,ibm),yvar(i_f,ibm)))  
         fulllps =fulllps-
     &         dlog(studt1(yf(i_f),nu,ymfull(i_f),yvarfull(i_f)))  
         
         nulllps = nulllps - dlog(studt1(yf(i_f),nu,ymean,yvar2))         
         avelps = avelps - dlog(xlps) 
 3400 enddo
      
      avelps  = avelps/dble(nf)
      bestlps = bestlps/dble(nf)
      fulllps = fulllps/dble(nf)
      nulllps = nulllps/dble(nf)
      
      write(fout,*)
      write(fout,*) 'lps with ',nf,' observations'
      write(fout,*)
      write(fout,3410) avelps
      write(fout,3411) bestlps
      write(fout,3413) fulllps
      write(fout,3414) nulllps
      write(fout,*)
 3410 format('bma lps............',g18.8)
 3411 format('best-model lps.....',g18.8)
 3413 format('full-model lps.....',g18.8)
 3414 format('null-model lps.....',g18.8)
c
c      
      return
      end
c

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  subroutine drawsample                                           ccc
ccc  created July 2006                                               ccc
ccc    Draws a non-repl sample of k elements out of n                ccc
ccc    Based on                                                      ccc
ccc    McLeod and Bellhouse (1983) Appl Statist 32(2): 182-184       ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c    
      subroutine drawsample(idum,n,k,sample,fail)
c     
      implicit real*8(a-h,o-z)
      
      logical fail
      integer idum,n,k,mxk
      parameter(mxk=50)
      integer sample(mxk)
c     
      if (n.lt.k) then
         fail = .true.
         return
      else
         fail = .false.
      endif
c
c start assigning first k elements to the sample
c  
      do 10 i=1,k
         sample(i) = i
 10   enddo  
      if (n.eq.k) return
c
c then for i>k, draw u = U[0,1] and make
c j <- 1 + int(u*i)
c if j is less than k, then i replaces the jth sample element
c     
      do 20 i=k+1,n
         u = ran2(idum)
         j = 1 + int(u*i)
         if (j.le.k) sample(j) = i
 20   enddo
      
      return
      end          
c   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccc  subroutine ols                                                  ccc
ccc  created July 2006                                               ccc
ccc    Computes OLS estimates, z is demeaned                         ccc
ccc    Here we already know kreg, so no need for maxk                ccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c    
      subroutine ols(fout,regname,nobs,kreg,z,ztz,y,avey,fail)
c      
      implicit real*8(a-h,o-z)
      parameter(maxk=104,maxn=88)
c dummy args:      
      real*8 z(maxn,maxk),ztz(maxk,maxk),y(maxn)  
      character*12 regname(maxk)
      integer fout
      logical fail
c local vars:      
      real*8 betahat(maxk),betastd(maxk),
     &      zy(maxk),zzi(maxk,maxk),p(maxk)
c start:      
      fail = .false.
c
c  OLS estimates
c
      do 531 j=1,kreg
         zy(j) = 0.0d0
         do 530 i=1,nobs
            zy(j) = zy(j) + z(i,j)*y(i)
  530    enddo
  531 enddo
c
      call choldc(ztz,kreg,maxk,p,fail)
      if (fail) then
         write(fout,*) 'OLS: ztz is singular!!'
         write(*,*) 'OLS: ztz is singular!!'
         return
      endif
      call cholsl(ztz,kreg,maxk,p,zy,betahat)
      
      s2 = 0.0d0
      do 550 i=1,nobs
         zbeta = 0.0d0
         do 540 j=1,kreg
            zbeta = zbeta + z(i,j)*betahat(j)
  540    enddo
         s2 = s2 + (y(i) - avey - zbeta)**2
 550  enddo
      s2 = s2/dble(nobs-kreg)
      
      do 601 i=1,kreg
         do 600 j=1,kreg
            zzi(i,j) = 0.0d0
 600     enddo        
         zzi(i,i) = 1.0d0
 601  enddo
      
      do 610 j=1,kreg
         call cholsl(ztz,kreg,maxk,p,zzi(1,j),zzi(1,j))
 610  enddo
      
      write(fout,*)   
      write(fout,*) 'Full-model OLS estimates'
      write(fout,*) 
      write(fout,612)
 612  format(t33, 'Betahat      St.Dev.    ratio')
      call barra('_',65,fout)
      write(fout,620) avey, s2
      write(fout,*) 
      do 615 j=1,kreg
         betastd(j) = dsqrt(s2*zzi(j,j))
         write(fout,622) j,regname(j),betahat(j),betastd(j),
     &         abs(betahat(j)/betastd(j))
 615  enddo
 620  format(' Intercept',      t31,2g12.4)
 622  format(' X(',i2,'): ',a12,t31,2g12.4,f6.1)
      call barra('_',65,fout)
      
      return
      end    
      
cEOF 7/19/2006
                                                                      
                                                                       
