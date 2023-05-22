!Note: I put the shrug (¯\_(ツ)_/¯) emoji everywhere I'm not certain. 
!In a proper academia sense: "The correctness of these is left as an exercise to the reader."
!TODO:
!  Parallelism in OMP (probably will not use MPI)
!     You will need to first list all the transitions without calculating them and then iterate over them
!

module constants

      double precision, parameter :: dblTol = 1d-10

!     the one to rule them all
      double complex, parameter :: iu = (0d0,1d0)
      
!      
      !custom physical constant values, not finished
      integer, parameter :: const_count = 100 !TODO
      character(10) :: const_chars(const_count)
      
     ! data const_chars//
 
!     all in SI   
      double precision, parameter :: h = 6.62607004d-34
      double precision, parameter :: h1953 = 6.6252d-34
      double precision, parameter :: pi = 3.14159265359
      double precision, parameter :: hbar = h/2d0/pi
      double precision, parameter :: el = 1.60217662d-19
      double precision, parameter :: el_AU = 1d0
      double precision, parameter :: mne = 1.67492749804d-27
      double precision, parameter :: mpr = 1.67262192369d-27
      double precision, parameter :: mpr_AU = 1.83615D3
      double precision, parameter :: mpr_amu = 1.007276466621
      double precision, parameter :: mel = 9.1093837015d-31
      double precision, parameter :: cc = 299792458d0
      double precision, parameter :: cc_AU=137.036d0
      double precision, parameter :: NA = 6.02214076D23
      double precision, parameter :: BohrM = hbar/2.0d0/mel*el
      double precision, parameter :: BohrM_AU = 0.5d0
      double precision, parameter :: BohrM_CGS = BohrM_AU/cc_AU
      double precision, parameter :: gs = 2.002318
      double precision, parameter :: NucM = hbar/2.0d0/mpr*el
      double precision, parameter :: NucM_AU=0.5d0/(mpr_AU)
      double precision, parameter :: NucM_CGS=0.5d0/(mpr_AU)/137.0359991 ! ¯\_(ツ)_/¯
      double precision, parameter :: BohrR = 0.52917721092d0 !in Angstrom
      double precision, parameter :: permi0 = 8.8541878128D-12
      double precision, parameter :: perme0 = 1.25663706212D-6
      double precision, parameter :: kb = 1.380649d-23
      double precision, parameter :: Rgas = 8.31446261815324
      double precision, parameter :: SI_2_Debye=1d0/3.33564D-30
      double precision, parameter :: DEBYE_2_SI=3.33564D-30
      
!     energy conversions in SI
      double precision, parameter :: Hz_2_cm = 3.33564d-11
      double precision, parameter :: GHz_2_cm = 1d9*Hz_2_cm
      double precision, parameter :: cm_2_Hz = 1d0/hz_2_cm
      double precision, parameter :: cm_2_GHz = 1d0/GHz_2_cm
      
      double precision, parameter :: amu_2_au = 1d0/1.82289d0
      
!     conversions to AU
      double precision, parameter :: SI2AU_debye = 1.0d0/2.541746913d0
      double precision, parameter :: SI2AU_BohrMagneton = BohrM_AU/BohrM
      double precision, parameter :: SI2AU_Hz    = 1.51983D-16
      double precision, parameter :: SI2AU_Kg    = 1d0/9.10939D-31
      double precision, parameter :: SI2AU_m     = 1d0/5.29177d-11 !bohr radius
      
!     conversions to "SI"
      double precision, parameter :: AU2SI_debye = 2.541746913d0
      double precision, parameter :: AU2SI_BohrMagneton = 2.0d0
      double precision, parameter :: AU2SI_Kg = SI2AU_Kg**(-1)
     
      character(1), parameter :: axes(0:4) = ['N','X','Y','Z','N']
      
      contains
      
      function SI2SI(from,to,val)result(res)
         character(*) from,to
         double precision val,res
      end function SI2SI
      
      function SI2AU(from,val)result(res)
         character(*) from
         double precision val,res
      end function SI2AU
      
      function AU2SI(from,val)result(res)
         character(*) from
         double precision val,res
      end function AU2SI
      
      subroutine setConstant(const,val,info)
         character(*),intent(in) :: const
         integer,intent(out) :: info
         double precision,intent(in) :: val
         
         
      end subroutine setConstant

      function Levi()result(res)
         integer res(0:4,0:4,0:4)
         
         res=0
         do i = 1,3
            do j = 1,3
               do k = 1,3
                  if(i/=j .and. j/=k .and. k/=i)then
                     if(i-j==-1 .or. i-j==2)then
                        res(i,j,k)=1
                     else if(i-j==-2 .or. i-j==1)then
                        res(i,j,k)=-1
                     end if
                  end if
               end do
            end do
         end do
         
      end function Levi
      
end module constants

module molecule
   use constants
   implicit none
   
   private atsy,atnn,atsy_iso,atnn_iso,symbol_len,iso_count,symbol_iso_len
   public molecule_spec,mol
   
   integer,parameter :: symbol_iso_len=6,symbol_len=2,iso_count=4
   
   character(symbol_len) :: atsy(118)
   integer :: atnn(118),atnn_iso(iso_count),atan_iso(iso_count)
   character(symbol_iso_len) :: atsy_iso(iso_count)
   

   type Mol
      character(:),allocatable :: mol_name
      double precision dip(3),gt(3,3),quad(3,3),abc(3),eps(3,3)
      double precision, allocatable :: R(:,:)
   end type Mol
   
   type(mol) :: molecule_spec
   
   data atNN/1,4,6,10,11,12,14,16,19,20,23,24,27,28,31,32,35,40,&
    39,40,45,48,51,52,55,56,59,59,63,65,70,73,75,79,80,84,85,88,&
    89,91,93,96,97,101,103,106,108,112,115,119,122,128,127,131,133,&
    137,139,140,141,144,145,150,152,157,159,162,165,167,169,173,175,&
    178,181,184,186,190,192,195,197,201,204,207,209,209,210,222,223,&
    226,227,232,231,238,237,244,243,247,247,251,254,257,256,254,256,&
    261,262,263,262,265,266,0,0,0,0,0,0,0,0,0/
   
   data atsy/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',&
   'Na','Mg','Al','Si',' P',' S','Cl','Ar',&
   ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
            'Ga','Ge','As','Se','Br','Kr',&
   'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',&
            'In','Sn','Sb','Te',' I','Xe',&
   'Cs','Ba','La',&
                 'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',&
                 'Er','Tm','Yb','Lu',&
   'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',&
            'Tl','Pb','Bi','Po','At','Rn',&
   'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf',&
   'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg',&
   'Cn','Nh','Fl','Mc','Lv','Ts','Og'/
   
   data atsy_iso/'     D','     T','   C13','   C14'/
   data atan_iso/1,1,6,6/
   data atnn_iso/2,3,13,14/
   
   contains
   
   function kdelta(a,b)result(res)
      INTEGER a,b,res
      
      
      if(a==b)then
         res=1
      else
         res=0
      end if
   end function kdelta
   
   function ROTMAT_2_ANGLES(rot) result(arr)
      double precision, intent(in) :: rot(3,3)
      double precision :: arr(3)
      
      arr(1)=atan2(rot(3,2),rot(3,3))
      arr(2)=atan2(-rot(3,1),sqrt(rot(3,2)**2+rot(3,3)**2))
      arr(3)=atan2(rot(2,1),rot(1,1))
      arr=arr/pi*180.0d0
   end function ROTMAT_2_ANGLES
   
   function CALC_Q(n,z,r)result(q)
      INTEGER :: i,j,k,n,z(n)
      double precision :: r(3,n),q(3,3)
      
      q=0d0
      do i = 1,n
         do j=1,3
            do k=1,3
               q(j,k)=q(j,k)+dble(z(i))*(dble(kdelta(j,k))*DOT_PRODUCT(r(:,i),r(:,i))-r(j,i)*r(k,i))
            end do
         end do
      end do
   end function CALC_Q
   
   FUNCTION CALC_INERT(n,m,r) RESULT(t)
      INTEGER :: i,n
      double precision,intent(in) :: m(n),r(3,n)
      double precision :: t(3,3)
      
      t=0d0
      do i = 1,n
         t(1,1)=t(1,1) + m(i)*(r(2,i)**2+r(3,i)**2)
         t(2,2)=t(2,2) + m(i)*(r(1,i)**2+r(3,i)**2)
         t(3,3)=t(3,3) + m(i)*(r(2,i)**2+r(1,i)**2)
         t(1,2)=t(1,2) - m(i)*r(1,i)*r(2,i)
         t(1,3)=t(1,3) - m(i)*r(1,i)*r(3,i)
         t(2,3)=t(2,3) - m(i)*r(2,i)*r(3,i)
      end do
      t(3,2)=t(2,3)
      t(3,1)=t(1,3)
      t(2,1)=t(1,2)
      
   END FUNCTION CALC_INERT
   
   FUNCTION CALC_COM(n,m,r) RESULT (COM)
      INTEGER :: i,n
      double precision,intent(in) :: m(n),r(3,n)
      double precision :: COM(3)
      
      COM=0d0
      do i = 1,n
         COM=COM+m(i)*r(:,i)
      end do
      
      COM=COM/SUM(m)
   END FUNCTION CALC_COM

   SUBROUTINE SHIFT_COORD(n,r,c)
      INTEGER n,i
      double precision,intent(inout) :: r(3,n)
      double precision,intent(in) :: c(3)
      
      do i = 1,n
         r(:,i)=r(:,i)-c
      end do
   END SUBROUTINE SHIFT_COORD
   
   
   function Shift_Dipole_do(dip,cm,charge)result(res)
      double precision :: dip(3),cm(3),res(3),charge
      res=dip-charge*CM
   end function Shift_Dipole_do
   
   function Shift_GT_do(gt,dip,cm,ama,n,r0,r1)result(res)
      double precision :: gt(3,3),cm(3),res(3,3),r0(3,n),r1(3,n),t0(3,3),t1(3,3),ama(n),dip(3),dp
      integer :: i,j,n
      dp = DOT_PRODUCT(dip,CM)
      t0 = calc_inert(n,ama,r0)
      t1 = calc_inert(n,ama,r1)
      do i = 1,3
         do j = 1,3
            res(i,j)=gt(i,j)*t0(j,j)-mpr_amu*(2d0*kdelta(i,j)*dp-dip(i)*CM(j)-CM(i)*dip(j))
            res(i,j)=res(i,j)/t1(j,j)
         end do
      end do
   end function Shift_GT_do
   
   function SHIFT_QUAD_TR_do(quad,dip,cm,charge)result(res)
      double precision :: quad(3,3),dip(3),cm(3),res(3,3),charge
      integer :: i,j
      
      do i = 1,3
         do j = 1,3
            res(i,j)=quad(i,j)-CM(i)*dip(j)-CM(j)*dip(i)
         end do
      end do
   end function SHIFT_QUAD_TR_do
   
   function shift_eps_do(eps,cm)result(res)
      double precision :: eps(3,3),cm(3),res(3,3)
      res=eps
   end function shift_eps_do
   
   function findSymbol(arr,n,m,str)result(res)
      integer res,i,n,m
      
      character(m) arr(n)
      character(*) :: str
      
      res=0
      do i = 1,n
         if(str==trim(adjustl(arr(i))))then
            res=i
            exit
         end if
      end do
   end function findSymbol
   
   subroutine Symbol2ANNN(symbol,AN,NN)
      integer, intent(out) :: AN,NN
      character(*),intent(in) :: symbol
      
      integer :: idx
      
      AN=-1
      NN=-1
      idx = findSymbol(atsy,size(atsy,dim=1),symbol_len,symbol)
      if(idx/=0)then
         AN=idx
         NN=atnn(idx)
         return
      else
         idx = findSymbol(atsy_iso,size(atsy_iso,dim=1),symbol_iso_len,symbol)
         if(idx/=0)then
            AN=atan_iso(idx)
            NN=atnn_iso(idx)
         else
            !write(6,'(A,A)')'Unknown atom symbol ',symbol
         end if
      end if
      
   end subroutine Symbol2ANNN
   
   function getIsotopeMass(AN,NN)result(m)
      integer,intent(in) :: AN,NN
      INTEGER :: n,atomic,nucleon,length,atomic2,i
      double precision :: m,mass
      double precision :: masses(30)
      double precision, parameter :: carbon12 = 12.011
      character(4) :: symbol
      open(99,file='isotopes.par')
!19       continue
19       read(99,1996)atomic,symbol,nucleon,mass
         if(atomic==AN)then
            if(nucleon==NN)then
               m=mass
               goto 20
            end if
            read(99,1996)atomic2,symbol,nucleon,mass
            do while(atomic2==atomic)
               if(nucleon==NN)then
                  m=mass
                  goto 20
               end if
               read(99,1996)atomic2,symbol,nucleon,mass
            end do
            
            write(6,'(A,I3,A,I3)') 'Could not find atom with AN/NN ',AN,' ',NN
            stop
         else
            goto 19
         end if
         
!21       if(atomic==n .and.abs(mass/m-1d0)<1d-1)goto 20
      !goto 19
20    close(99)
1996  format(BN,I4,A4,I5,F15.11)
      m = m/(carbon12/12d0)
   end function getIsotopeMass
   
end module molecule

module arrayOperations
   
   implicit none
   contains
   
   pure subroutine SwapEl(arr,n,a,b)
      Integer,intent(in) :: n,a,b
      double precision,intent(inout) :: arr(n)
      double precision :: buf
      
      if(a==b)return
      
      buf=arr(a)
      arr(a)=arr(b)
      arr(b)=buf
   end subroutine SwapEl
   
   pure subroutine SwapEl_int(arr,n,a,b)
      Integer,intent(in) :: n,a,b
      Integer,intent(inout) :: arr(n)
      Integer :: buf
      
      if(a==b)return
      
      buf=arr(a)
      arr(a)=arr(b)
      arr(b)=buf
   end subroutine SwapEl_int
   
   subroutine Reverse_arr_int(arr,n)
      INTEGER n,i
      INTEGER arr(n),buf
      
      if(MODULO(n,2)==0)then
         do i = 1,n/2
            buf=arr(i)
            arr(i)=arr(n-i+1)
            arr(n-i+1)=buf
         end do
      else
         do i = 1,(n-1)/2
            buf=arr(i)
            arr(i)=arr(n-i+1)
            arr(n-i+1)=buf
         end do
      end if
   end subroutine Reverse_arr_int
   
   subroutine Reorder_arr_2d(arr,n,order,column)
      INTEGER n,order(n),i
      double precision arr(n,n),buf(n,n)
      logical column
      
      buf=arr
      
      if(column)then
         do i = 1,n
            arr(:,i)=buf(:,order(i))
         end do
      else
         do i = 1,n
            arr(i,:)=buf(order(i),:)
         end do
      end if
   end subroutine Reorder_arr_2d
   
   subroutine Reorder_arr_col(arr,colSize,n,order)
      INTEGER ColSize,n,order(colSize),i
      double precision arr(colSize,n),buf(colSize,n)
      
      buf=arr
      do i = 1,colSize
         arr(i,:)=buf(order(i),:)
      end do
   end subroutine Reorder_arr_col
   
end module arrayOperations

module wigner
   use constants

   implicit none
   INTEGER, parameter :: factSize_8=169,factsize_4=34
   INTEGER W3Jsize
   double precision,allocatable :: factorials(:)
   integer W3J_count
   !double complex,parameter :: iu=(0d0,1d0)
   logical factFlag
   
   public factorials,factFlag
   
   contains
   
   subroutine COMPUTE_FACTORIALS(kin)
      INTEGER i,top
      double precision kin
      
      if(kind(kin)==8)then
         top=factSize_8
      else if(kind(kin)==4)then
         top=factsize_4
      end if
      allocate(factorials(top))
      
      do i = 1,top
         factorials(i)=fact(i)
      end do
      
      
      ! do i = 171,1754
         ! facts_10(i) = FACT_10(i)
      ! end do
      
      factFlag=.true.

   end subroutine COMPUTE_FACTORIALS
   
   !"adapted" for half-integers
   pure function W6j_R(j1,j2,j3,j4,j5,j6)result(res)
      double precision,intent(in) :: j1,j2,j3,j4,j5,j6
      double precision :: res
      
      
      if(.not.(Triangle_rule_r(j1,j2,j3).and.Triangle_rule_r(j1,j5,j6).and.Triangle_rule_r(j2,j4,j6).and.Triangle_rule_r(j3,j4,j5)))then
         res=0d0
         return
      end if
      
      !if(j2==0.5 .or. j2==1.0)then
      !if(.not.isInt)
      res=(-1.0)**(j1+j2+j4+j5)*RacahR(j1,j2,j5,j4,j3,j6)
      !end if
   end function W6j_R
   
   !watchout for precision
   !this might fail later
   !I am relying on the observation that all half-integers can be EXACTLY represented by double-precisions
   !ISBN: 978-0-306-47123-0, Devanathan, Angular Momentum Techniques ...
   !TODO: check correctness....at some point....in my life
   pure function Racah_exact(a,b,c,d,e,f) result(res)
      double precision,intent(in) :: a,b,c,d,e,f
      double precision res
      !yuhp, totally safe, will not be indeterminate whatsoever
      if(.not.(Triangle_rule_r(a,b,e).and.Triangle_rule_r(c,d,e).and.Triangle_rule_r(a,c,f).and.Triangle_rule_r(b,d,f)))then
         res=0d0
         return
      end if
      
      if(e==0.5)then
         if(c==d+0.5)then
            if(a==b+0.5)then
               res=(-1.0)**(b+d-f)*sqrt(((b+d+f+2.0)*(b+d-f+1.0))/((2.0*b+1.0)*(2.0*b+2.0)*(2.0*d+1.0)*(2.0*d+2.0)))
            else if(a==b-0.5)then
               res=(-1.0)**(b+d-f)*sqrt(((f-b+d+1.0)*(f+b-d))/(2.0*b*(2.0*b+1.0)*(2.0*d+1.0)*(2.0*d+2.0)))
            end if
         else if(c==d-0.5)then
            if(a==b+0.5)then
               res=(-1.0)**(b+d-f)*sqrt(((f+b-d+1.0)*(f-b+d))/((2.0*b+1.0)*(2.0*b+2.0)*2.0*d*(2.0*d+1.0)))
            else if(a==b-0.5)then
               res=(-1.0)**(b+d-f-1.0)*sqrt(((b+d+f+1.0)*(-f+b+d))/(2.0*b*(2.0*b+1)*2.0*d*(2.0*d+1.0)))
            end if
         end if
      else if(e==1.0)then
         if(c==d+1.0)then
            if(a==b+1.0)then
               res=(-1.0)**(b+d-f)*sqrt(((b+d+f+3.0)*(b+d+f+2.0)*(b+d-f+2.0)*(b+d-f+1.0))/(4.0*(2.0*b+1.0)*(2.0*b+3)*(b+1.0)*(2.0*d+1.0)*(2.0*d+3.0)*(d+1.0)))
            else if(a==b)then
               res=(-1.0)**(b+d-f)*sqrt(((b+d+f+2.0)*(b+d-f+1.0)*(-b+d+f+1.0)*(b-d+f))/(4.0*b*(2.0*b+1.0)*(b+1.0)*(2.0*d+1.0)*(d+1.0)*(2.0*d+3.0)))
            else if(a==b-1.0)then
               res=(-1.0)**(b+d-f)*sqrt(((b-d+f)*(b-d+f-1.0)*(-b+d+f+2.0)*(-b+d+f+1.0))/(4.0*(2.0*b+1)*(2.0*b-1)*b*(d+1.0)*(2.0*d+1)*(2.0*d+3.0)))
            end if
         else if(c==d)then
            if(a==b+1.0)then
               res=(-1.0)**(b+d-f)*sqrt(((b+d+f+2)*(b-d+f+1)*(b+d-f+1)*(-b+d+f))/(4.0*(2.0*b+1)*(2.0*b+3)*(b+1.0)*(2.0*d+1)*d*(d+1.0)))
            else if(a==b)then
               !in Devanathan this one is wrong
               !see Angular Momentum in Quantum Physics, Biedenharn, Louck
               res=(-1.0)**(b+d-f-1.0)*(b*(b+1)+d*(d+1)-f*(f+1))/sqrt((4.0*b*(2.0*b+1)*(b+1)*(2.0*d+1)*d*(d+1)))
            else if(a==b-1.0)then
               res=(-1.0)**(b+d-f-1.0)*sqrt(((b+d+f+1)*(b+d-f)*(b-d+f)*(-b+d+f+1))/(4.0*(2.0*b+1)*(2.0*b-1)*b*(d+1)*d*(2.0*d+1)))
            end if
         else if(c==d-1.0)then
            if(a==b+1.0)then
               res=(-1.0)**(b+d-f)*sqrt(((-b+d+f)*(-b+d+f-1)*(b-d+f+2)*(b-d+f+1))/(4.0*d*(2.0*b+1)*(2.0*b+3)*(b+1)*(2.0*d-1)*(2.0*d+1)))
            else if(a==b)then
               res=(-1.0)**(b+d-f-1.0)*sqrt(((b+d+f+1)*(b-d+f+1)*(-b+d+f)*(b+d-f))/(4.0*b*d*(2.0*b+1)*(b+1)*(2.0*d+1)*(2.0*d-1)))
            else if(a==b-1.0)then
               res=(-1.0)**(b+d-f)*sqrt(((b+d+f+1)*(b+d+f)*(b+d-f)*(b+d-f-1))/((4.0*b*d*(2.0*b+1)*(2.0*b-1)*(2.0*d+1)*(2.0*d-1))))
            end if
         end if
      endif
   end function Racah_exact
      
   pure function Racah_exact_int(a,b,c,d,e,f)result(res)
      Integer,intent(in) :: a,b,c,d,e,f
      double precision res
      
      res=Racah_exact(real(a,kind(res)),real(b,kind(res)),real(c,kind(res)),real(d,kind(res)),real(e,kind(res)),real(f,kind(res)))
   end function Racah_exact_int
      
   pure function Racah_tr_int(a,b,c)result(res)
      Integer,intent(in) :: a,b,c
      double precision :: res
      INTEGER one,two,three,four
      one=a+b-c
      two=a-b+c
      three=-a+b+c
      four=a+b+c+1
      
      if(one<0 .or. two<0 .or. three<0)then
         res=0d0
      else
         if(one>170 .or. two>170 .or. three>170 .or. four>170)then
            res = sqrt((FACT_10((one))*FACT_10((two))*FACT_10((three)))/FACT_10((four)))
         else
            res = sqrt((FACT((one))*FACT((two))*FACT((three)))/FACT((four)))
         end if
      end if
   end function Racah_tr_int
   
   pure function W6j(j1,j2,j3,j4,j5,j6)result(res)
      Integer,intent(in) :: j1,j2,j3,j4,j5,j6
      double precision :: res
      INTEGER m1,m2,m3,m4,m5,m6
      
      res=(-1)**(j1+j2+j4+j5)*Racah(j1,j2,j5,j4,j3,j6)
   end function W6j
   
   !RACAH COEFFICIENT
   pure function Racah(a,b,c,d,e,f) result(res)
      Integer,intent(in) :: a,b,c,d,e,f
      INTEGER :: x,lmax,lmin,buf1,buf2,buf3,buf4,buf5,buf6,buf7,buf8
      double precision :: res,buf
      
      if((d==1) .and. e /= 1)then
         res=(-1d0)**(a+d-e-f)*Racah_exact_int(e,c,b,f,d,a)
         return
      else if((f==1) .and. e /= 1)then
         res=Racah_exact_int(a,c,b,d,f,e)
         return
      else if(e==1)then
         res=Racah_exact_int(a,b,c,d,e,f)
         return
      end if
      buf=0d0
      
      res=racah_tr_int(a,b,e)*racah_tr_int(c,d,e)*racah_tr_int(a,c,f)*racah_tr_int(b,d,f)
      if(res==0d0)return
      
      lMin=max(a+b+e,c+d+e,a+c+f,b+d+f)
      lMax=min(a+b+c+d,a+d+e+f,b+c+e+f)
      
      x=lmin
      do while(x<=lmax)
         buf1=x+1
         buf2=x-a-b-e
         buf3=x-c-d-e
         buf4=x-a-c-f
         buf5=x-b-d-f
         buf6=a+b+c+d-x
         buf7=a+d+e+f-x
         buf8=b+c+e+f-x
         if(buf1>170 .or. buf2>170 .or. buf3>170 .or. buf4 > 170 .or. buf5>170 .or. buf6>170 .or. buf7>170 .or. buf8>170)then
            buf=buf + (-1)**(x+a+b+c+d)*FACT_10(buf1)/(FACT_10(buf2)*FACT_10(buf3)*FACT_10(buf4)*FACT_10(buf5)*FACT_10(buf6)*FACT_10(buf7)*FACT_10(buf8)) 
         else
            buf=buf + (-1)**(x+a+b+c+d)*FACT(buf1)/(FACT(buf2)*FACT(buf3)*FACT(buf4)*FACT(buf5)*FACT(buf6)*FACT(buf7)*FACT(buf8)) 
         end if
         x=x+1
      end do
      
      res=res*buf
   end function Racah
   
   pure function Racah_tr(a,b,c)result(res)
      double precision,intent(in) :: a,b,c
      double precision res,one,two,three,four
      
      if(a+b-c<0 .or. a-b+c<0 .or. -a+b+c<0)then
         res=0d0
         return
      end if
      
      one=a+b-c
      two=a-b+c
      three=-a+b+c
      four=a+b+c+1
      ! if(.not.(isInt(one).and.isInt(two).and.isInt(three).and.isInt(four)))then
         ! res=Huge(1.0d0)+1d0
         ! return
      ! end if
      if(one>170 .or. two>170 .or. three>170 .or. four>170)then
         res = sqrt((FACT_10(int(one))*FACT_10(int(two))*FACT_10(int(three)))/FACT_10(int(four)))
      else
         res = sqrt((FACT(int(one))*FACT(int(two))*FACT(int(three)))/FACT(int(four)))
      end if
   end function Racah_tr  
     
   !RACAH COEFFICIENT FOR HALF-INTEGERS  
   pure function RacahR(a,b,c,d,e,f) result(res)
      double precision,intent(in) :: a,b,c,d,e,f
      double precision :: res,buf,x,lmax,lmin
      
      if(e==1.0 .or. e==0.5)then
         res=Racah_exact(a,b,c,d,e,f)
         return
      else if((f==1.0 .or. f==0.5))then
         res=Racah_exact(a,c,b,d,f,e)
         return
      else if((d==1.0 .or. d==0.5))then
         res=(-1.0)**(a+d-e-f)*Racah_exact(e,c,b,f,d,a)
         return
      else if((a==1.0 .or. a==0.5))then
         res=(-1.0)**(a+d-e-f)*Racah_exact(e,b,c,f,a,d)
         return
      else if((b==1.0 .or. b==0.5))then
         res=(-1.0)**(b+c-e-f)*Racah_exact(a,e,f,d,b,c)
         return
      else if((c==1.0 .or. c==0.5))then
         res=(-1.0)**(b+c-e-f)*Racah_exact(a,f,e,d,c,b)
         return
      end if
      buf=0d0
      
      res=racah_tr(a,b,e)*racah_tr(c,d,e)*racah_tr(a,c,f)*racah_tr(b,d,f)
      if(res==0d0)return
      
      lMin=max(a+b+e,c+d+e,a+c+f,b+d+f)
      lMax=min(a+b+c+d,a+d+e+f,b+c+e+f)
      
      x=lmin
      do while(x<=lmax)
         buf=buf + (-1)**(x+a+b+c+d)*FACTR(x+1)/(FACTR(x-a-b-e)*FACTR(x-c-d-e)*FACTR(x-a-c-f)*FACTR(x-b-d-f)*FACTR(a+b+c+d-x)*FACTR(a+d+e+f-x)*FACTR(b+c+e+f-x)) 
         x=x+1
      end do
      
      res=res*buf
   end function RacahR
   
   pure function FactR(n)result(res)
      double precision, intent(in) :: n
      double precision res,x
      
      if(n==0.0)then
         res=1.0
      else if(factFlag)then
         res=factorials(NINT(n))
      else
         res=1
         x=2
         do while(x<=n)
            res=x*res
         end do
      end if
   end function FactR
      
   pure function IsHalfInt(num)result(res)
      double precision,intent(in) :: num
      double precision bufR
      INTEGER numI
      logical res
      
      bufr=num*2.0
      if(isInt(bufr))then
         res=mod(NINT(bufr),2)==1
      else
         res=.false.
      end if
   end function IsHalfInt
   
   pure function IsInt(num)result(res)
      double precision,intent(in) :: num
      INTEGER numI
      logical res
      
      numI=NINT(num)
      res=numI-num==0.0
   end function IsInt
   
   pure function CG_SR(N,S,J,MN,MS,MJ)result(res)
      Integer,INTENT(IN) :: N,MN
      double precision,INTENT(IN) :: S,J,MJ,MS
      double precision :: res
      
      res=0d0
      if(.not. Triangle_rule_sr(N,S,J) .or. Mn+Ms/=Mj .or. abs(mn)>n .or. abs(ms)>s .or. abs(mj)>j)return
      if(S==0.5)then
         if(j==N+S)then
            if(ms==0.5)then
               res=sqrt((N+mj+0.5)/(2d0*N+1))
            else if(ms==-0.5)then
               res=sqrt((N-mj+0.5)/(2d0*N+1))
            end if
         else if(j==N-S)then
            if(ms==0.5)then
               res=-sqrt((N-mj+0.5)/(2d0*N+1))
            else if(ms==-0.5)then
               res=sqrt((N+mj+0.5)/(2d0*N+1))
            end if
         end if
      else if(S==1)then
         if(J==dble(N+1))then
            if(ms==1d0)then
               res=sqrt((N+mj)*(N+mj+1d0)/((2d0*N+1d0)*(2d0*N+2d0)))
            else if(ms==0d0)then
               res=sqrt((N-mj+1d0)*(N+mj+1d0)/((2d0*N+1d0)*(N+1d0)))
            else if(ms==-1d0)then
               res=sqrt((N-mj)*(N-mj+1d0)/((2d0*N+1d0)*(2d0*N+2d0)))
            end if
         else if(J==dble(N))then
            if(ms==1d0)then
               res=-sqrt((N+mj)*(N-mj+1d0)/(2d0*N*(N+1d0)))
            else if(ms==0d0)then
               res=dble(mj)/sqrt(N*(N+1d0))
            else if(ms==-1d0)then
               res=sqrt((N-mj)*(N+mj+1d0)/(2d0*N*(N+1d0)))
            end if
         else if(J==dble(N-1))then
            if(ms==1d0)then
               res=sqrt((N-mj)*(N-mj+1d0)/(2d0*N*(2d0*N+1d0)))
            else if(ms==0d0)then
               res=-sqrt((N-mj)*(N+mj)/(N*(2d0*N+1d0)))
            else if(ms==-1d0)then
               res=sqrt((N+mj)*(N+mj+1d0)/(2d0*N*(2d0*N+1)))
            end if
         end if
         return
      end if
   end function CG_SR
   
   !too lazy to rewrite everything
   !N is actually (=>) J
   !J => N
   !Mn => Mj
   !Mj => Mn
   !CLEBSCH GORDAN COEFFICIENT FOR UNCOUPLING INTO MOLECULE FIXED AXES
   !|N K S J M > = SUM(Ms)[CG_SR_REV() |J Kj M > | S -Ks> ]
   !probably unreliably, I dont know how to get this uncoupling
   !di Lauro might be wrong here
   pure function CG_SR_REV(N,S,J,MN,MS,MJ)result(res)
      Integer,INTENT(IN) :: J,MJ
      double precision,INTENT(IN) :: S,N,MN,MS
      double precision :: res
      
      res=0d0
      if(.not. Triangle_rule_sr(J,S,N) .or. Mn+Ms/=Mj .or. abs(mn)>n .or. abs(ms)>s .or. abs(mj)>j)return
      if(S==0.5)then
         if(j==N+S)then
            if(ms==0.5)then
               res=sqrt((N+mj+0.5)/(2d0*N+1))
            else if(ms==-0.5)then
               res=sqrt((N-mj+0.5)/(2d0*N+1))
            end if
         else if(j==N-S)then
            if(ms==0.5)then
               res=-sqrt((N-mj+0.5)/(2d0*N+1))
            else if(ms==-0.5)then
               res=sqrt((N+mj+0.5)/(2d0*N+1))
            end if
         end if
      else if(S==1)then
         if(J==dble(N+1))then
            if(ms==1d0)then
               res=sqrt((N+mj)*(N+mj+1d0)/((2d0*N+1d0)*(2d0*N+2d0)))
            else if(ms==0d0)then
               res=sqrt((N-mj+1d0)*(N+mj+1d0)/((2d0*N+1d0)*(N+1d0)))
            else if(ms==-1d0)then
               res=sqrt((N-mj)*(N-mj+1d0)/((2d0*N+1d0)*(2d0*N+2d0)))
            end if
         else if(J==dble(N))then
            if(ms==1d0)then
               res=-sqrt((N+mj)*(N-mj+1d0)/(2d0*N*(N+1d0)))
            else if(ms==0d0)then
               res=dble(mj)/sqrt(N*(N+1d0))
            else if(ms==-1d0)then
               res=sqrt((N-mj)*(N+mj+1d0)/(2d0*N*(N+1d0)))
            end if
         else if(J==dble(N-1))then
            if(ms==1d0)then
               res=sqrt((N-mj)*(N-mj+1d0)/(2d0*N*(2d0*N+1d0)))
            else if(ms==0d0)then
               res=-sqrt((N-mj)*(N+mj)/(N*(2d0*N+1d0)))
            else if(ms==-1d0)then
               res=sqrt((N+mj)*(N+mj+1d0)/(2d0*N*(2d0*N+1)))
            end if
         end if
         return
      end if
   end function CG_SR_REV
   
   !CLEBSCH GORDAN COEFFICIENT FOR HALF-INTEGERS
   pure function CG_R(j1,j2,j3,m1,m2,m3)result(res)
      double precision A,B,C,buf
      double precision,INTENT(IN) :: j1,j2,j3,m1,m2,m3
      double precision res,i
      
      if((m1+m2)/=m3 .or. .not.TRIANGLE_RULE_R(j1,j2,j3)) then
         res = 0d0
         return 
      end if
      
      
      if(j2==0d0)then
         if(j1==j3.and.m1==m3)then
            res=1
         else
            res=0
         end if
         return
      else if(j2==1d0)then
         if(j3==j1+1)then
            if(m2==1d0)then
               res=sqrt((j1+m3)*(j1+m3+1d0)/((2d0*j1+1d0)*(2d0*j1+2d0)))
            else if(m2==0)then
               res=sqrt((j1-m3+1d0)*(j1+m3+1d0)/((2d0*j1+1d0)*(j1+1d0)))
            else if(m2==-1)then
               res=sqrt((j1-m3)*(j1-m3+1d0)/((2d0*j1+1d0)*(2d0*j1+2d0)))
            end if
         else if(j3==j1)then
            if(m2==1)then
               res=-sqrt((j1+m3)*(j1-m3+1d0)/(2d0*j1*(j1+1d0)))
            else if(m2==0)then
               res=dble(m3)/sqrt(j1*(j1+1d0))
            else if(m2==-1)then
               res=sqrt((j1-m3)*(j1+m3+1d0)/(2d0*j1*(j1+1d0)))
            end if
         else if(j3==j1-1)then
            if(m2==1)then
               res=sqrt((j1-m3)*(j1-m3+1d0)/(2d0*j1*(2d0*j1+1d0)))
            else if(m2==0)then
               res=-sqrt((j1-m3)*(j1+m3)/(j1*(2d0*j1+1d0)))
            else if(m2==-1)then
               res=sqrt((j1+m3)*(j1+m3+1d0)/(2d0*j1*(2d0*j1+1)))
            end if
         end if
         return
      else if(j2==2d0)then
         if(m2==2d0)then
            if(j3==j1+2)then
               res=sqrt(((j1+m3-1d0)*(j1+m3)*(j1+m3+1)*(j1+m3+2))/((2*j1+1)*(2*j1+2)*(2*j1+3)*(2*j1+4)))
            else if(j3==j1+1)then
               res=-sqrt(((j1+m3-1d0)*(j1+m3)*(j1+m3+1)*(j1-m3+2))/(2*j1*(j1+1)*(j1+2)*(2*j1+1)))
            else if(j3==j1)then
               res=sqrt((3d0*(j1+m3-1)*(j1+m3)*(j1-m3+1)*(j1-m3+2))/((2*j1-1)*2*j1*(j1+1)*(2*j1+3)))
            else if(j3==j1-1)then
               res=-sqrt(((j1+m3-1d0)*(j1-m3)*(j1-m3+1)*(j1-m3+2))/(2*(j1-1)*j1*(j1+1)*(2*j1+1)))
            else if(j3==j1-2)then
               res=sqrt(((j1-m3-1d0)*(j1-m3)*(j1-m3+1)*(j1-m3+2))/((2*j1-2)*(2*j1-1)*2*j1*(2*j1+1)))
            end if
         else if(m2==1d0)then
            if(j3==j1+2)then
               res=sqrt(((j1-m3+2d0)*(j1+m3+2)*(j1+m3+1)*(j1+m3))/((2*j1+1)*(j1+1)*(2*j1+3)*(j1+2)))
            else if(j3==j1+1)then
               res=-(j1-2d0*m3+2)*sqrt(((j1+m3+1d0)*(j1+m3))/(2*j1*(2*j1+1)*(j1+1)*(j1+2)))
            else if(j3==j1)then
               res=(1d0-2*m3)*sqrt(((3*(j1-m3+1d0)*(j1+m3)))/((2*j1-1)*j1*(2*j1+2)*(2*j1+3)))
            else if(j3==j1-1)then
               res=(j1+2d0*m3-1)*sqrt(((j1-m3+1d0)*(j1-m3))/((j1-1)*j1*(2*j1+1)*(2*j1+2)))
            else if(j3==j1-2)then
               res=-sqrt(((j1-m3+1d0)*(j1-m3)*(j1-m3-1)*(j1+m3-1))/((j1-1)*(2*j1-1)*j1*(2*j1+1)))
            end if
         else if(m2==0d0)then
            if(j3==j1+2)then
               res=sqrt(((3d0*(j1-m3+2)*(j1-m3+1)*(j1+m3+2)*(j1+m3+1)))/((2*j1+1)*(2*j1+2)*(2*j1+3)*(j1+2)))
            else if(j3==j1+1)then
               res=m3*sqrt(((3d0*(j1-m3+1)*(j1+m3+1)))/(j1*(2*j1+1)*(j1+1)*(j1+2)))
            else if(j3==j1)then
               res=(3d0*m3**2-j1*(j1+1))/sqrt((2d0*j1-1)*j1*(j1+1)*(2*j1+3))
            else if(j3==j1-1)then
               res=-m3*sqrt(((3d0*(j1-m3)*(j1+m3)))/((j1-1)*j1*(2*j1+1)*(j1+1)))
            else if(j3==j1-2)then
               res=sqrt((3d0*(j1-m3)*(j1-m3-1)*(j1+m3)*(j1+m3-1))/((2*j1-2)*(2*j1-1)*j1*(2*j1+1)))
            end if
         else if(m2==-1d0)then
            if(j3==j1+2)then
               res=sqrt(((j1-m3+2d0)*(j1-m3+1)*(j1-m3)*(j1+m3+2))/((2*j1+1)*(j1+1)*(2*j1+3)*(j1+2)))
            else if(j3==j1+1)then
               res=(j1+2d0*m3+2)*sqrt(((j1-m3+1d0)*(j1-m3))/(j1*(2*j1+1)*(2*j1+2)*(j1+2)))
            else if(j3==j1)then
               res=(2d0*m3+1)*sqrt((3d0*(j1-m3)*(j1+m3+1))/((2*j1-1)*j1*(2*j1+2)*(2*j1+3)))
            else if(j3==j1-1)then
               res=-(j1-2d0*m3-1)*sqrt(((j1+m3+1d0)*(j1+m3))/((j1-1)*j1*(2*j1+1)*(2*j1+2)))
            else if(j3==j1-2)then
               res=-sqrt(((j1-m3-1d0)*(j1+m3+1)*(j1+m1)*(j1+m3-1))/((j1-1)*(2*j1-1)*j1*(2*j1+1)))
            end if
         else if(m2==-2d0)then
            if(j3==j1+2)then
               res=sqrt(((j1-m3-1d0)*(j1-m3)*(j1-m3+1)*(j1-m3+2))/((2*j1+1)*(2*j1+2)*(2*j1+3)*(2*j1+4)))
            else if(j3==j1+1)then
               res=sqrt(((j1-m3-1d0)*(j1-m3)*(j1-m3+1)*(j1+m3+2))/(j1*(2*j1+1)*(j1+1)*(2*j1+4)))
            else if(j3==j1)then
               res=sqrt((3d0*(j1-m3-1)*(j1-m3)*(j1+m3+1)*(j1+m3+2))/((2*j1-1)*(2*j1+2)*(2*j1+3)))
            else if(j3==j1-1)then
               res=sqrt(((j1-m3-1d0)*(j1+m3)*(j1+m3+1)*(j1+m3+2))/((j1-1)*j1*(2*j1+1)*(2*j1+2)))
            else if(j3==j1-2)then
               res=sqrt(((j1+m3-1d0)*(j1+m3)*(j1+m3+1)*(j1+m3+2))/((2*j1-2)*(2*j1-1)*2*j1*(2*j1+1)))
            end if
         end if
         return
      end if
      
      A = (FACTN(NINT(j1+j2-j3))*FACTN(NINT(j3+j1-j2))*FACTN(NINT(j2+j3-j1)))/(FACTN(NINT(j1+j2+j3+1)))
      B = FACTN(NINT(j1+m1))*FACTN(NINT(j1-m1))*FACTN(NINT(j2+m2))*FACTN(NINT(j2-m2))*FACTN(NINT(j3+m3))*FACTN(NINT(j3-m3))
      
      if((A==0d0).OR.(B==0d0))then
         res=0d0
         return
      end if
      
      C = 0
      buf = 0
      i = max(-j3+j2-m1,-j3+j1+m2,0d0)
      do while (i <= min(j1+j2-j3,j1-m1,j2+m2))
         C = (FACTN(NINT(j1+j2-j3-i))*FACTN(NINT(j1-m1-i))*FACTN(NINT(j2+m2-i))*FACTN(NINT(j3-j2+m1+i))*FACTN(NINT(j3-j1-m2+i)))
         buf = buf + (-1)**i/FACT(NINT(i))*C**(-1)
         i=i+1
      end do
      
      res = sqrt(dble(2*j3+1)*A*B)*buf
   end function CG_R
   
   pure function TRIANGLE_RULE(a,b,c) result(res)
      Integer, INTENT(IN) :: a,b,c
      logical res
      
      if(a+b>=c.AND.a+c>=b.AND.b+c>=a)then
         res = .true.
      else
         res = .false.
      end if
   end function TRIANGLE_RULE

   pure function TRIANGLE_RULE_R(a,b,c) result(res)
      double precision, intent(in) :: a,b,c
      logical res
      
      if(a+b>=c.AND.a+c>=b.AND.b+c>=a)then
         res = .true.
      else
         res = .false.
      end if
   end function TRIANGLE_RULE_R
   
   pure function TRIANGLE_RULE_SR(a,b,c) result(res)
      Integer, INTENT(IN) :: a
      double precision, intent(in) :: b,c
      logical res
      
      if(a+b>=c.AND.a+c>=b.AND.b+c>=a)then
         res = .true.
      else
         res = .false.
      end if
   end function TRIANGLE_RULE_SR
   
   !Clebsch Gordan coefficient defined by 3 Js and 3Ms, returns double
   pure function CG(j1,j2,j3,m1,m2,m3) result(res)
      double precision A,B,C,buf
      Integer,INTENT(IN) :: j1,j2,j3,m1,m2,m3
      double precision res
      INTEGER i
      
      if((m1+m2)/=m3 .or. .not.TRIANGLE_RULE(j1,j2,j3) .or. abs(m1) > J1 .or. abs(m2) > J2 .or. abs(m3) > J3) then
         res = 0d0
         return 
      end if
      
      
      if(j2==0)then
         if(j1==j3.and.m1==m3)then
            res=1
         else
            res=0
         end if
         return
      else if(j2==1)then
         if(j3==j1+1)then
            select case(m2)
               case(1)
                  res=sqrt((j1+m3)*(j1+m3+1d0)/((2d0*j1+1d0)*(2d0*j1+2d0)))
               case(0)
                  res=sqrt((j1-m3+1d0)*(j1+m3+1d0)/((2d0*j1+1d0)*(j1+1d0)))
               case(-1)
                  res=sqrt((j1-m3)*(j1-m3+1d0)/((2d0*j1+1d0)*(2d0*j1+2d0)))
               ! case default
                  ! STOP 'WRONG CG: J2=1 M=?'
            end select
         else if(j3==j1)then
            select case(m2)
               case(1)
                  res=-sqrt((j1+m3)*(j1-m3+1d0)/(2d0*j1*(j1+1d0)))
               case(0)
                  res=dble(m3)/sqrt(j1*(j1+1d0))
               case(-1)
                  res=sqrt((j1-m3)*(j1+m3+1d0)/(2d0*j1*(j1+1d0)))
               ! case default
                  ! STOP 'WRONG CG: J2=1 M=?'
            end select
         else if(j3==j1-1)then
            select case(m2)
               case(1)
                  res=sqrt((j1-m3)*(j1-m3+1d0)/(2d0*j1*(2d0*j1+1d0)))
               case(0)
                  res=-sqrt((j1-m3)*(j1+m3)/(j1*(2d0*j1+1d0)))
               case(-1)
                  res=sqrt((j1+m3)*(j1+m3+1d0)/(2d0*j1*(2d0*j1+1)))
               ! case default
                  ! STOP 'WRONG CG: J2=1 M=?'
            end select
         end if
         return
      else if(j2==2)then
         if(m2==2)then
            if(j3==j1+2)then
               res=sqrt(((j1+m3-1d0)*(j1+m3)*(j1+m3+1)*(j1+m3+2))/((2*j1+1)*(2*j1+2)*(2*j1+3)*(2*j1+4)))
            else if(j3==j1+1)then
               res=-sqrt(((j1+m3-1d0)*(j1+m3)*(j1+m3+1)*(j1-m3+2))/(2*j1*(j1+1)*(j1+2)*(2*j1+1)))
            else if(j3==j1)then
               res=sqrt((3d0*(j1+m3-1)*(j1+m3)*(j1-m3+1)*(j1-m3+2))/((2*j1-1)*2*j1*(j1+1)*(2*j1+3)))
            else if(j3==j1-1)then
               res=-sqrt(((j1+m3-1d0)*(j1-m3)*(j1-m3+1)*(j1-m3+2))/(2*(j1-1)*j1*(j1+1)*(2*j1+1)))
            else if(j3==j1-2)then
               res=sqrt(((j1-m3-1d0)*(j1-m3)*(j1-m3+1)*(j1-m3+2))/((2*j1-2)*(2*j1-1)*2*j1*(2*j1+1)))
            end if
         else if(m2==1)then
            if(j3==j1+2)then
               res=sqrt(((j1-m3+2d0)*(j1+m3+2)*(j1+m3+1)*(j1+m3))/((2*j1+1)*(j1+1)*(2*j1+3)*(j1+2)))
            else if(j3==j1+1)then
               res=-(j1-2d0*m3+2)*sqrt(((j1+m3+1d0)*(j1+m3))/(2*j1*(2*j1+1)*(j1+1)*(j1+2)))
            else if(j3==j1)then
               res=(1d0-2*m3)*sqrt(((3*(j1-m3+1d0)*(j1+m3)))/((2*j1-1)*j1*(2*j1+2)*(2*j1+3)))
            else if(j3==j1-1)then
               res=(j1+2d0*m3-1)*sqrt(((j1-m3+1d0)*(j1-m3))/((j1-1)*j1*(2*j1+1)*(2*j1+2)))
            else if(j3==j1-2)then
               res=-sqrt(((j1-m3+1d0)*(j1-m3)*(j1-m3-1)*(j1+m3-1))/((j1-1)*(2*j1-1)*j1*(2*j1+1)))
            end if
         else if(m2==0)then
            if(j3==j1+2)then
               res=sqrt(((3d0*(j1-m3+2)*(j1-m3+1)*(j1+m3+2)*(j1+m3+1)))/((2*j1+1)*(2*j1+2)*(2*j1+3)*(j1+2)))
            else if(j3==j1+1)then
               res=m3*sqrt(((3d0*(j1-m3+1)*(j1+m3+1)))/(j1*(2*j1+1)*(j1+1)*(j1+2)))
            else if(j3==j1)then
               res=(3d0*m3**2-j1*(j1+1))/sqrt((2d0*j1-1)*j1*(j1+1)*(2*j1+3))
            else if(j3==j1-1)then
               res=-m3*sqrt(((3d0*(j1-m3)*(j1+m3)))/((j1-1)*j1*(2*j1+1)*(j1+1)))
            else if(j3==j1-2)then
               res=sqrt((3d0*(j1-m3)*(j1-m3-1)*(j1+m3)*(j1+m3-1))/((2*j1-2)*(2*j1-1)*j1*(2*j1+1)))
            end if
         else if(m2==-1)then
            if(j3==j1+2)then
               res=sqrt(((j1-m3+2d0)*(j1-m3+1)*(j1-m3)*(j1+m3+2))/((2*j1+1)*(j1+1)*(2*j1+3)*(j1+2)))
            else if(j3==j1+1)then
               res=(j1+2d0*m3+2)*sqrt(((j1-m3+1d0)*(j1-m3))/(j1*(2*j1+1)*(2*j1+2)*(j1+2)))
            else if(j3==j1)then
               res=(2d0*m3+1)*sqrt((3d0*(j1-m3)*(j1+m3+1))/((2*j1-1)*j1*(2*j1+2)*(2*j1+3)))
            else if(j3==j1-1)then
               res=-(j1-2d0*m3-1)*sqrt(((j1+m3+1d0)*(j1+m3))/((j1-1)*j1*(2*j1+1)*(2*j1+2)))
            else if(j3==j1-2)then
               res=-sqrt(((j1-m3-1d0)*(j1+m3+1)*(j1+m1)*(j1+m3-1))/((j1-1)*(2*j1-1)*j1*(2*j1+1)))
            end if
         else if(m2==-2)then
            if(j3==j1+2)then
               res=sqrt(((j1-m3-1d0)*(j1-m3)*(j1-m3+1)*(j1-m3+2))/((2*j1+1)*(2*j1+2)*(2*j1+3)*(2*j1+4)))
            else if(j3==j1+1)then
               res=sqrt(((j1-m3-1d0)*(j1-m3)*(j1-m3+1)*(j1+m3+2))/(j1*(2*j1+1)*(j1+1)*(2*j1+4)))
            else if(j3==j1)then
               res=sqrt((3d0*(j1-m3-1)*(j1-m3)*(j1+m3+1)*(j1+m3+2))/((2*j1-1)*(2*j1+2)*(2*j1+3)))
            else if(j3==j1-1)then
               res=sqrt(((j1-m3-1d0)*(j1+m3)*(j1+m3+1)*(j1+m3+2))/((j1-1)*j1*(2*j1+1)*(2*j1+2)))
            else if(j3==j1-2)then
               res=sqrt(((j1+m3-1d0)*(j1+m3)*(j1+m3+1)*(j1+m3+2))/((2*j1-2)*(2*j1-1)*2*j1*(2*j1+1)))
            end if
         end if
         return
      end if
      
      A = (FACTN(j1+j2-j3)*FACTN(j3+j1-j2)*FACTN(j2+j3-j1))/(FACTN(j1+j2+j3+1))
      B = FACTN(j1+m1)*FACTN(j1-m1)*FACTN(j2+m2)*FACTN(j2-m2)*FACTN(j3+m3)*FACTN(j3-m3)
      
      if((A==0d0).OR.(B==0d0))then
         res=0d0
         return
      end if
      
      C = 0
      buf = 0
      do i = max(-j3+j2-m1,-j3+j1+m2,0),min(j1+j2-j3,j1-m1,j2+m2),1
         C = (FACTN(j1+j2-j3-i)*FACTN(j1-m1-i)*FACTN(j2+m2-i)*FACTN(j3-j2+m1+i)*FACTN(j3-j1-m2+i))
         buf = buf + (-1)**i/FACT(i)*C**(-1)
      end do
      
      res = (dble(2*j3+1)*A*B)**(0.5)*buf
      if(res>huge(1d0))then
         res=0d0
      end if
   end function CG
   
   !FACT(n) function checking for zero, returns double
   pure function FACTN(n) result (res)
      Integer,intent(in) :: n
      double precision :: res
      res = 0d0
      if(n>=0) then
         res = FACT(n)
      end if
   
   end function FACTN

   !Factorial of double type value, returns double
   pure function FACT(n) result (res)
      double precision :: res
      Integer,intent(in) :: n
      INTEGER :: i
      
      if((n==0).OR.(n==1))then
         res = 1d0
         return
      end if
      
      if(factFlag) then
         res = factorials(n)
         return
      end if
      
      res = 1d0
            
      do i= 1,n,1
         res = res*i
      end do
            
   end function FACT
   
   pure function FACT_10(n) result(res)
      real(10) :: res
      Integer,intent(in) :: n
      INTEGER :: i
      
      if((n==0).OR.(n==1))then
         res = 1d0
         return
      end if
      
      if(factFlag) then
         res = factorials(n)
         return
      end if
      
      res = 1.0_10
            
      do i= 1,n,1
         res = res*i
      end do
            
   end function FACT_10
end module wigner

module rotor
   use constants
   use wigner
   implicit none
   
      
   contains

   !J. Quant. Spectrosc. Radiat. Transfer Vol. 21, pp. 309-313
   !FOR VOIGT BAND PROFILE, UNTESTED, UNUSED, COPIED HERE IF YOU WANT IT
   SUBROUTINE CPFI2(X,Y,WR,WI)
      ! ROUTINE COMPUTES THE REAL (WR) AND IMAGINARY (WI) PARTS OF THE COMPLEX
      ! PROBABILITY FUNCTION W(Z)=EXP(-Z**2)*ERFC(-I*Z) IN THE UPPER HALF-PLANE Z=X+i*Y (I.E. FOR Y>=O)
      ! MAXIXUM RELATIVE ERROR OF WR IS <2*10**(-6), THAT OF WI <5*10**(-6)
      double precision :: T(6),C(6),S(6),X,Y,WR,WI,Y1,Y2,D1,D2,D3,D4,R1,R2,D,R,Y3
      INTEGER :: I
      !DIMENSION T(6),C(6),S(6)
      DATA(T(I),I=1,6)/.314240376,.947788391,1.59768264,2.27950708,3.02063703, &
    & 3.8897249/,(C(I),I=1,6)/1.01172805,-.75197147,1.2557727E-2, &
    & 1.00220082E-2,-2.42068135E-4,5.00848061E-7/,(S(I),I=1,6)/1.393237, &
    & .231152406,-.155351466,6.21836524E-3,9.19082986E-5,-6.27525958E-7/
      WR=0.
      WI=9.
      Y1=Y+1.5
      Y2=Y1*Y1
      ! IF HIGH ~ELATIVE ACCURACY IN REGION II IS NOT REQUIRED° THE FOLLOWING
      ! 19 LI~ES ~AY BE DELETED
      IF((Y.GT.0.85).OR.(ABS(X).LT.18.1*y+1.65))GOTO 2
      ! ***REGION II
      IF(ABS(X).LT.12.)WR=EXP(-X*X)
      Y3=Y+3.
      DO I=1,6
         R=X-T(I)
         R2=R*R
         D=1./(R2+Y2)
         D1=Y1*D
         D2=R*D
         WR=WR+Y*(C(I)*(R*D2-1.5*D1)+S(I)*Y3*D2)/(R2+2.25)
         R=X+T(I)
         R2=R*R
         D=1./(R2+Y2)
         D3=Y1*D
         D4=R*D
         WR=WR+Y*(C(I)*(R*D4-1.5*D3)-S(I)*Y3*D4)/(R2+2.25)
         WI=WI+C(I)*(D2+D4)+S(I)*(D1-D3)
      end do
      RETURN
      ! ***REGION I
2     DO I=1,6
         R=X-T(I)
         D=1./(R*R+Y2)
         D1=Y1*D
         D2=R*D
         R=X+T(I)
         D=1./(R*R+Y2)
         D3=Y1*D
         D4=R*D
         WR=WR+C(I)*(D1+D3)-S(I)*(D2-D4)
         WI=WI+C(I)*(D2+D4)+S(I)*(D1-D3)
      end do
      RETURN
   END SUBROUTINE CPFI2
   
   !matrix determinant
   !Shamelessly stolen from:
   !http://www.rosettacode.org/wiki/Determinant_and_permanent#Fortran
   recursive function DET(a,n,permanent) result(accumulation)
       ! setting permanent to 1 computes the permanent.
       ! setting permanent to -1 computes the determinant.
       Integer, intent(in) :: n, permanent
       double precision, dimension(n,n), intent(in) :: a
       double precision, dimension(n-1, n-1) :: b
       double precision :: accumulation
       INTEGER :: i, sgn
       if (n == 1) then
         accumulation = a(1,1)
       else
         accumulation = 0
         sgn = 1
         do i=1, n
           b(:, :(i-1)) = a(2:, :i-1)
           b(:, i:) = a(2:, i+1:)
           accumulation = accumulation + sgn * a(1, i) * det(b, n-1, permanent)
           sgn = sgn * permanent
         enddo
       endif
   end function DET   
   
!====================================
!WAVE FUNCTIONS, OPERATORS AND DIAGONALIZATIONS
!====================================

   !J+ (ladder + operator) eigenvalue for molecule fixed axis
   pure function Jplus_EVAL(J,K) result (res)
      Integer,intent(in) :: J,K
      double precision res
      
      if(abs(k)>j .or. K==-J)then
         res=0d0
      else
         res=sqrt(dble(j*(j+1)-k*(k-1)))
         !res=sqrt(dble(j*(j+1)-k*(k+1)))
      end if
   end function Jplus_EVAL
   
   function Jplus_EVAL_shift(J,K) result (res)
      Integer,intent(inout) :: J,K
      double precision res
      
      if(abs(k)>j .or. K==-J)then
         res=0d0
      else
         res=sqrt(dble(j*(j+1)-k*(k-1)))
         K=K-1
         !res=sqrt(dble(j*(j+1)-k*(k+1)))
      end if
   end function Jplus_EVAL_shift
   
   !J-
   pure function Jminus_EVAL(J,K) result (res)
      Integer,intent(in) :: J,K
      double precision res
      
      if(abs(k)>j .or. K==J)then
         res=0d0
      else
         res=sqrt(dble(j*(j+1)-k*(k+1)))
         !res=sqrt(dble(j*(j+1)-k*(k-1)))
      end if
   end function Jminus_EVAL
   
   function Jminus_EVAL_shift(J,K) result (res)
      Integer,intent(inout) :: J,K
      double precision res
      
      if(abs(k)>j .or. K==J)then
         res=0d0
      else
         res=sqrt(dble(j*(j+1)-k*(k+1)))
         K=K+1
         !res=sqrt(dble(j*(j+1)-k*(k-1)))
      end if
   end function Jminus_EVAL_shift
   
   !DIRECTION COSINE < N2 |COSE| N1 >
   pure function DIR_COSE_N_FROM_SPH(N2,N1)result(res)
      Integer,intent(in) :: N2,N1
      double precision :: res
      res = sqrt((2d0*N1+1d0)/(2d0*N2+1d0))
   end function DIR_COSE_N_FROM_SPH
   
   pure function DIR_COSE_N(N2,N1)result(res)
      Integer,intent(in) :: N2,N1
      double precision :: res
      select case(N2-N1)
         case(1)
            res=1d0/(4d0*N2*sqrt(4d0*N2**2-1))
         case(0)
            res=1d0/(4d0*N2*(N2+1))
         case(-1)
            res=1d0/(4d0*N1*sqrt(4d0*N1**2-1))
         case default
            res=0d0
      end select
   end function DIR_COSE_N
   
   !DIRECTION COSINE IN MOLECULAR AXIS < N2 K2 |COSE(molaxis)| N1 K1 >
   pure function DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1,molaxis)result(res)
      Integer,intent(in) :: N1,N2,k1,k2
      character(1),intent(in) :: molaxis
      double complex :: res
   
      select case(molaxis)
         case('X')
            res=1d0/dsqrt(2d0)*(CG(N1,1,N2,K1,-1,K2)-CG(N1,1,N2,K1,1,K2))
         case('Y')
            res=-iu/dsqrt(2d0)*(CG(N1,1,N2,K1,-1,K2)+CG(N1,1,N2,K1,1,K2))
         case('Z')
            res=CG(N1,1,N2,K1,0,K2)
         case default
            res=0d0
      end select
   end function DIR_COSE_NK_FROM_SPH
   
   
   !DIRECTION COSINE IN SPACE-FIXED AXIS < N2 K2 |COSE(labaxis)| N1 K1 >
   pure function DIR_COSE_NM_FROM_SPH(N2,M2,N1,M1,labaxis)result(res)
      Integer,intent(in) :: N1,N2,m1,m2
      character(1),intent(in) :: labaxis
      double complex :: res
      
      select case(labaxis)
         case('X')
            res=1d0/dsqrt(2d0)*(CG(N1,1,N2,M1,-1,M2)-CG(N1,1,N2,M1,1,M2))
         case('Y')
            res=iu/dsqrt(2d0)*(CG(N1,1,N2,M1,-1,M2)+CG(N1,1,N2,M1,1,M2))
         case('Z')
            res=CG(N1,1,N2,M1,0,M2)
         case default
            res=0d0
      end select
   end function DIR_COSE_NM_FROM_SPH
 
   pure function DIR_COSE_JM_FROM_SPH(J2,M2,J1,M1,labaxis)result(res)
      double precision,intent(in) :: J1,J2,m1,m2
      character(1),intent(in) :: labaxis
      double complex :: res
      
      select case(labaxis)
         case('X')
            res=1d0/dsqrt(2d0)*(CG_R(J1,1d0,J2,M1,-1d0,M2)-CG_R(J1,1d0,J2,M1,1d0,M2))
         case('Y')
            res=iu/dsqrt(2d0)*(CG_r(J1,1d0,J2,M1,-1d0,M2)+CG_R(J1,1d0,J2,M1,1d0,M2))
         case('Z')
            res=CG_R(J1,1d0,J2,M1,0d0,M2)
         case default
            res=0d0
      end select
   end function DIR_COSE_JM_FROM_SPH
 
   pure function DIR_COSE_NKM_FROM_SPH(N2,K2,M2,N1,K1,M1,F,alpha)result(res)
      integer,intent(in) :: N2,K2,M2,N1,K1,M1
      character(1),intent(in) :: F,alpha
      double complex res
      
      res=DIR_COSE_NM_FROM_SPH(N2,M2,N1,M1,F)*DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1,alpha)*DIR_COSE_N_FROM_SPH(N2,N1)
   end function DIR_COSE_NKM_FROM_SPH
 
   pure function DIR_COSE_NSJ(N2,S2,J2,N1,S1,J1)result(res)
      integer, intent(in) :: N2,N1
      double precision, intent(in) :: S2,J2,S1,J1
      double complex :: phase
      double complex :: res
      
      res=0d0
      if(S2/=S1)return
      phase = (-1d0)**(N2-N1)
      res = phase * (-1d0)**(S1+1-J2-N2) * dsqrt((2*J1+1)*(2*N1+1)) * RacahR(dble(N2),dble(N1),J2,J1,1d0,S2)
      
   end function DIR_COSE_NSJ
 
   !alpha is the direction cosine axis
   !beta is AM operator axis
   pure function DIR_COSE_NK_N_matel(N2,K2,N1,K1,alpha,beta)result(res)
      integer,intent(in) :: N2,K2,N1,K1
      double complex res,eig
      character(1),intent(in) :: alpha,beta
      
      
      select case(beta)
         case('X')
            res=0.5d0*(Jplus_EVAL(N1,K1)*DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1-1,alpha)+Jminus_EVAL(N1,K1)*DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1+1,alpha))
         case('Y')
            res=-iu*0.5d0*(Jplus_EVAL(N1,K1)*DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1-1,alpha)-Jminus_EVAL(N1,K1)*DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1+1,alpha))
         case('Z')
            res=DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1,alpha)*K1
         case default
            res=0d0
      end select
   end function DIR_COSE_NK_N_matel
 
   function Salpha_matel_NKSJM(N2,K2,S2,J2,M2,N1,K1,S1,J1,M1,axis,space)result(res)
      integer :: N2,N1,K2,K1,F,Mn1,Mn2,ms1_idx,ms2_idx
      logical space
      double precision :: S2,J2,M2,S1,J1,M1,ms2,ms1
      character(1) :: axis
      double complex :: res
      character(1),parameter :: axes(3)=['X','Y','Z']
      integer :: ks2_idx,Ks1_idx,S_mult
      double precision :: ks2_rev,ks1_rev,kj2,kj1,phase
      
      res=0d0
	   if(S2/=S1)return
      S_mult=NINT(2*S2)
      if(space)then
         do ms2_idx=0,S_mult
            ms2=ms2_idx-S2
            Mn2=NINT(M2-Ms2)
            do ms1_idx=0,S_mult
               ms1=ms1_idx-S1
               Mn1=NINT(M1-ms1)
               do F = 1,3
                  res=res+CG_SR(N2,S2,J2,Mn2,Ms2,M2)*CG_SR(N1,S1,J1,Mn1,Ms1,M1)*SAlpha_matel_space(S2,Ms2,S1,Ms1,axes(F))*DIR_COSE_NM_FROM_SPH(N2,Mn2,N1,Mn1,axes(F))
               end do
            end do
         end do
         res=res*DIR_COSE_N_FROM_SPH(N2,N1)*DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1,axis)
      else
         do ks2_idx=0,S_mult
            ks2_rev=ks2_idx-S2
            kj2=K2-ks2_rev
            do ks1_idx=0,S_mult
               ks1_rev=ks1_idx-S1
               kj1=kj2
               phase=(-1d0)**(N1-J1-Ks1_rev+N2-J2-Ks2_rev)
               !phase=1d0
               res=res-phase*CG_SR_REV(J2,S2,N2,Kj2,Ks2_rev,K2)*CG_SR_REV(J1,S1,N1,Kj1,Ks1_rev,K1)*Salpha_matel_SKs_rev(S2,Ks2_rev,s1,ks1_rev,axis)
            end do
         end do
      end if
   end function Salpha_matel_NKSJM
   
   function Salpha_matel_SKs_rev(S2,Ks2_rev,s1,ks1_rev,axis)result(res)
      double precision :: S2,S1,Ks2_rev,ks1_rev
      character(1) axis
      double complex :: res
      
      res=0d0
      if(abs(Ks2_rev)>s2 .or. abs(ks1_rev)>s1)return
      select case(axis)
         case('X')
            res=0.5d0*(Splus_matel_rev(S2,Ks2_rev,S1,ks1_rev)+Sminus_matel_rev(S2,Ks2_rev,S1,ks1_rev))
         case('Y')
            res=-iu*0.5d0*(Splus_matel_rev(S2,Ks2_rev,S1,ks1_rev)-Sminus_matel_rev(S2,Ks2_rev,S1,ks1_rev))
         case('Z')
            if(Ks2_rev==ks1_rev)res=ks1_rev
         case default
            res=0d0
      end select
   end function Salpha_matel_SKs_rev
   
   
   function Salpha_matel_SKs(S2,Ks2,s1,ks1,axis)result(res)
      double precision :: S2,S1,Ks2,ks1
      character(1) axis
      double complex :: res
      
      res=0d0
      select case(axis)
         case('X')
            res=0.5d0*(Splus_matel(S2,Ks2,S1,Ks1)+Sminus_matel(S2,Ks2,S1,Ks1))
         case('Y')
            res=-iu*0.5d0*(Splus_matel(S2,Ks2,S1,Ks1)-Sminus_matel(S2,Ks2,S1,Ks1))
         case('Z')
            if(ks2==ks1)res=ks1
         case default
            res=0d0
      end select
   end function Salpha_matel_SKs
   
   
   function Nalpha_matel_NKSJM(N2,K2,S2,J2,M2,N1,K1,S1,J1,M1,axis)result(res)
      integer :: N2,N1,K2,K1
      double precision :: S2,J2,M2,S1,J1,M1
      character(1) :: axis
      double complex :: res,other
      
      integer :: Ms2_idx,Ms1_idx,mult_S,mn2,mn1
      double precision :: Ms2,Ms1
      
      res=0d0
      if(N2/=N1.or.S2/=S1.or.m1/=m2)return
      
      ! mult_s=NINT(2*S2+1)
      ! do ms2_idx=0,mult_S-1
         ! ms2=ms2_idx-S2
         ! mn2=NINT(M2-Ms2)
         !! do ms1_idx=0,mult_S-1
            ! ms1=ms2
            ! Mn1=NINT(M1-Ms1)
            ! res=res+CG_SR(N2,S2,J2,Mn2,Ms2,M2)**2*Nalpha_matel_NKM(N2,K2,Mn2,N1,K1,Mn1,axis)
         !!end do
      ! end do
      
      res = Nalpha_matel_NKM(N2,K2,1,N1,K1,1,axis)
      !res=other
   end function Nalpha_matel_NKSJM
   
   function Nalpha_matel_NKM(N2,K2,M2,N1,K1,M1,axis)result(res)
      integer :: N2,K2,M2,N1,K1,M1
      character(1) :: axis
      double complex :: res
      
      res=0d0
      if(N2/=N1.or.m1/=m2)return
      select case(axis)
         case('X')
            res=0.5d0*(Jplus_matel(N2,K2,N1,K1)+Jminus_matel(N2,K2,N1,K1))
         case('Y')
            res=-iu*0.5d0*(Jplus_matel(N2,K2,N1,K1)-Jminus_matel(N2,K2,N1,K1))
         case('Z')
            if(K2==K1)res=K1
         case default
            res=0d0
      end select
   end function Nalpha_matel_NKM
   
   function Jplus_matel(N2,K2,n1,k1)result(res)
      integer :: n2,k2,n1,k1
      double precision :: res
      res=0d0
      if(abs(k1)>n1 .or. abs(k2)>n2)return
      if(N2==N1.and.K2==K1-1)then
         res=sqrt(dble(N1*(N1+1)-K1*(K1-1)))
      else
         res=0d0
      end if
   end function Jplus_matel
   
   function Jminus_matel(N2,K2,n1,k1)result(res)
      integer :: n2,k2,n1,k1
      double precision :: res
      res=0d0
      if(abs(k1)>n1 .or. abs(k2)>n2)return
      if(N2==N1.and.K2==K1+1)then
         res=sqrt(dble(N1*(N1+1)-K1*(K1+1)))
      else
         res=0d0
      end if
   end function Jminus_matel
   
   pure function Jplus_eval_space(J,M)result(res)
      Integer,intent(in) :: J,M
      double precision res
      if(abs(M)>j)then
         res=0d0
      else
         res=sqrt(dble(j*(j+1)-M*(M+1)))
         !res=sqrt(dble(j*(j+1)-k*(k+1)))
      end if
   end function Jplus_eval_space
   
   pure function Jminus_eval_space(J,M)result(res)
      Integer,intent(in) :: J,M
      double precision res
      if(abs(M)>j)then
         res=0d0
      else
         res=sqrt(dble(j*(j+1)-M*(M-1)))
         !res=sqrt(dble(j*(j+1)-k*(k+1)))
      end if
   end function Jminus_eval_space
   
   pure function Jq_sph_eval(J,K,q)result(res)
      Integer,intent(in) :: J,K,q
      double precision res
      
      select case(q)
         case(-1)
            res=Jminus_sph_EVAL(J,K)
         case(0)
            res=dble(K)
         case(1)
            res=Jplus_sph_eval(J,K)
      end select
   end function Jq_sph_eval
   
   pure function Jplus_sph_eval(J,K) result (res)
      Integer,intent(in) :: J,K
      double precision res
      if(abs(k)>j)then
         res=0d0
      else
         res=-1.0d0/sqrt(2.0d0)*Jplus_EVAL(J,K)
      end if
   end function Jplus_sph_EVAL
   
   pure function Jminus_sph_eval(J,K) result (res)
      Integer,intent(in) :: J,K
      double precision res
      if(abs(k)>j)then
         res=0d0
      else
         res=1.0d0/sqrt(2.0d0)*Jminus_EVAL(J,K)
      end if
   end function Jminus_sph_EVAL
   
   
   !At the point of writing these function below,
   !I was not sure whether all the operators have the exact same eigenvalues
   !in both the space and molecule-fixed axis systems
   !<S2 KS2| Sa |S1 KS1>
   pure function Splus_matel(S2,KS2,S1,KS1)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      double precision res
      if(S2==S1 .and. KS2==KS1+1)then
         res=sqrt(S1*(S1+1d0)-KS1*(KS1+1d0))
      else
         res = 0d0
      end if
   end function Splus_matel
   
   pure function Sminus_matel(S2,KS2,S1,KS1)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      double precision res
      if(S2==S1 .and. KS2==KS1-1)then
         res=sqrt(S1*(S1+1d0)-KS1*(KS1-1d0))
      else
         res = 0d0
      end if
   end function Sminus_matel
   
   pure function Splus_matel_space(S2,S1,MS2,MS1)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      double precision res

      if(S2==S1 .and. MS2==MS1+1)then
         res=sqrt(S1*(S1+1d0)-MS1*(MS1+1d0))
      else
         res = 0d0
      end if
   end function Splus_matel_space
   
   pure function Sminus_matel_space(S2,S1,MS2,MS1)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      double precision res
      if(S2==S1 .and. MS2==MS1-1)then
         res=sqrt(S1*(S1+1d0)-MS1*(MS1-1d0))
      else
         res = 0d0
      end if
   end function Sminus_matel_space
      
   pure function Sminus_matel_rev(S2,KSt2,S1,Kst1)result(res)   
      double precision,intent(in) :: S2,KSt2,S1,KSt1
      double precision res
      if(Kst2==kst1+1 .and. S2==S1)then
         res=-sqrt(S2*(S2+1d0)-KSt1*(Kst1+1d0))
      else
         res=0d0
      end if
   end function Sminus_matel_rev
      
   pure function Splus_matel_rev(S2,KSt2,S1,KSt1)result(res)  
      double precision,intent(in) :: S2,KSt2,S1,KSt1
      double precision res
      if(Kst2==kst1-1 .and. S2==S1)then
         res=-sqrt(S2*(S2+1d0)-KSt1*(Kst1-1d0))
      else
         res=0d0
      end if
   end function Splus_matel_rev   
      
   pure function Sx_matel_rev(S2,KSt2,S1,KSt1)result(res)
      double precision,intent(in) :: S2,KSt2,S1,KSt1
      double complex res
      res=0.5d0*(Splus_matel_rev(S2,KSt2,S1,KSt1)+Sminus_matel_rev(S2,KSt2,S1,KSt1))
   end function Sx_matel_rev
   
   pure function Sy_matel_rev(S2,KSt2,S1,KSt1)result(res) 
      double precision,intent(in) :: S2,KSt2,S1,KSt1
      double complex res
      res=-0.5d0*iu*(Splus_matel_rev(S2,KSt2,S1,KSt1)-Sminus_matel_rev(S2,KSt2,S1,KSt1))
   end function Sy_matel_rev 
      
   pure function Sz_matel_rev(S2,KSt2,S1,KSt1)result(res)
      double precision,intent(in) :: S2,KSt2,S1,KSt1
      double precision res
      if(KSt2==kst1 .and. s2==s1 .and. abs(kst2)<=s2)then
         res=kst2
      else
         res=0d0
      end if
   end function Sz_matel_rev   
      
   pure function Salpha_matel_rev(S2,KSt2,S1,KSt1,axis)result(res)
      double precision,intent(in) :: S2,KSt2,S1,KSt1
      character(1),intent(in) :: axis
      double complex res
      
      if(abs(Kst2)>S2 .or. abs(KSt1)>S1)then
         res=0d0
         return
      end if
      select case(axis)
         case('X')
            res=Sx_matel_rev(S2,KSt2,S1,KSt1)
         case('Y')
            res=Sy_matel_rev(S2,KSt2,S1,KSt1)
         case('Z')
            res=Sz_matel_rev(S2,KSt2,S1,KSt1)
         case default
            res=0d0
      end select
   end function Salpha_matel_rev
      
   function CommRel(J1,K1,M1,J2,K2,M2)result(code)
      Integer,intent(in) :: J1,K1,M1,J2,K2,M2
      INTEGER i,j,code
      character(1), parameter :: axes(3) = ['X','Y','Z']
      double precision left,right
      
      code=0
      
      
      if(abs(left-right)>1d-5)then
         code=1
         write(6,*)'WARNING: COMMUTATION RELATION NOT FULFILLED: '//axes(i)//axes(j)
         !goto 20
      end if
      
20    continue
   end function CommRel
   
   pure function SAlpha_matel(S2,S1,KS2,KS1,alpha)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      double complex res
      character(1),intent(in) :: alpha
      
      if(abs(Ks2)>S2 .or. abs(KS1)>S1)then
         res=0d0
         return
      end if
      select case(alpha)
         case('X')
            res=Sx_matel(S2,S1,KS2,KS1)
         case('Y')
            res=Sy_matel(S2,S1,KS2,KS1)
         case('Z')
            res=Sz_matel(S2,S1,KS2,KS1)
         case default
            res=0d0
      end select
   end function SAlpha_matel
      
   pure function Sp_matel_rev(S2,S1,KS2,KS1,p)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      Integer,intent(in) :: p
      double complex res
      
      if(abs(KS2)>S2 .or. abs(KS1)>S1)then
         res=0d0
         return
      end if
      
      select case(p)
         case(-1)
            if(S1==S2 .and. Ks2+1==Ks1)then
               res=-1d0/sqrt(2d0)*sqrt(S2*(S2+1)-Ks1*(Ks1+1))
            else
               res=0d0
            end if
         case(0)
            if(S1==S2 .and. Ks2==Ks1)then
               res=-Ks2
            else
               res=0d0
            end if
         case(1)
            if(S1==S2 .and. Ks2-1==Ks1)then
               res=1d0/sqrt(2d0)*sqrt(S2*(S2+1)-Ks1*(Ks1-1))
            else
               res=0d0
            end if
         case default
            res=0d0
      end select
   end function Sp_matel_rev
   
   pure function Sp_matel_space(S2,S1,MS2,MS1,p)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      Integer,intent(in) :: p
      double complex res
      
      if(abs(Ms2)>S2 .or. abs(MS1)>S1)then
         res=0d0
         return
      end if
      
      select case(p)
         case(-1)
            !res=1d0/sqrt(2d0)*(Sx_matel_space(S2,S1,Ms2,Ms1)-iu*Sy_matel_space(S2,S1,Ms2,Ms1))
            res=1d0/sqrt(2d0)*Sminus_matel_space(S2,S1,MS2,MS1)
         case(0)
            res=Sz_matel_space(S2,S1,Ms2,Ms1)
         case(1)
            !res=-1d0/sqrt(2d0)*(Sx_matel_space(S2,S1,Ms2,Ms1)+iu*Sy_matel_space(S2,S1,Ms2,Ms1))
            res=-1d0/sqrt(2d0)*Splus_matel_space(S2,S1,MS2,MS1)
         case default
            res=0d0
      end select
   end function Sp_matel_space
   
   pure function SAlpha_matel_space(S2,MS2,S1,MS1,alpha)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      double complex res
      character(1),intent(in) :: alpha
      
      if(abs(Ms2)>S2 .or. abs(MS1)>S1)then
         res=0d0
         return
      end if
      select case(alpha)
         case('X')
            res=Sx_matel_space(S2,S1,MS2,MS1)
         case('Y')
            res=Sy_matel_space(S2,S1,MS2,MS1)
         case('Z')
            res=Sz_matel_space(S2,S1,MS2,MS1)
         case default
            res=0d0
      end select
   end function SAlpha_matel_space
   
   pure function Sx_matel(S2,S1,KS2,KS1)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      double complex res
      res=iu*0.5d0*(Splus_matel(S2,S1,KS2,KS1)-Sminus_matel(S2,S1,KS2,KS1))
   end function Sx_matel
   
   pure function Sy_matel(S2,S1,KS2,KS1)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      double complex res
      res=0.5d0*(Sminus_matel(S2,S1,KS2,KS1)+Splus_matel(S2,S1,KS2,KS1))
   end function Sy_matel
   
   pure function Sz_matel(S2,S1,KS2,KS1)result(res)
      double precision,intent(in) :: S2,KS2,S1,KS1
      double complex res
      
      if(S1==S2 .and. KS2==KS1 .and. abs(KS2)<=S2)then
         res=KS2
      else
         res=0d0
      end if
   end function Sz_matel
   
   pure function Sy_matel_space(S2,S1,MS2,MS1)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      double complex res
      res=-iu*0.5d0*(Splus_matel_space(S2,S1,MS2,MS1)-Sminus_matel_space(S2,S1,MS2,MS1))
   end function Sy_matel_space
   
   pure function Sx_matel_space(S2,S1,MS2,MS1)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      double complex res
      res=0.5d0*(Sminus_matel_space(S2,S1,MS2,MS1)+Splus_matel_space(S2,S1,MS2,MS1))
   end function Sx_matel_space
   
   pure function Sz_matel_space(S2,S1,MS2,MS1)result(res)
      double precision,intent(in) :: S2,MS2,S1,MS1
      double complex res
      
      if(S1==S2 .and. MS2==MS1 .and. abs(MS2)<=S2)then
         res=MS2
      else
         res=0d0
      end if
   end function Sz_matel_space
   
   pure function S_matel(S2,S1)result(res)
      double precision,intent(in) :: S2,S1
      double precision res
      if(S2==S1)then
         res=sqrt(S1*(S1+1))
      else
         res=0d0
      end if
   end function S_matel
   
   
   recursive subroutine DIAG_SYM(Eval,Evec,dime,info,ordered)
      Integer, INTENT(IN) :: dime
      double precision, INTENT(OUT) :: Eval(dime)
      double precision, INTENT(INOUT) :: Evec(dime,dime)
      Integer, INTENT(OUT) :: info
      logical, intent(in) :: ordered
      double precision Work(7*dime,7*dime),evalc(dime),evecL(dime,dime),evecR(dime,dime)
      real(8) evec_4(dime,dime),eval_4(dime)
      
      info=0
      if(kind(eval(1))==8)then
      !LAPACK function, WORK array is set as large as possible, these eigenvalues are in ascending order
         if(ordered)then
            call dsyev('V','L',dime,evec,dime,eval,work,7*dime*7*dime,info)
         else
         !Another LAPACK function, these eigenvalues are not ordered (or perhaps they are shuffled ¯\_(ツ)_/¯, who knows)
            call dgeev('N','V',dime,evec,dime,eval,evalc,evecL,dime,evecR,dime,work,10*dime,info)
            eval=eval
            evec=evecR
         end if
      else if(kind(eval(1))==4)then
         if(ordered)then
            call ssyev('V','L',dime,evec,dime,eval,work,7*dime*7*dime,info)
         else
            call sgeev('N','V',dime,evec,dime,eval,evalc,evecL,dime,evecR,dime,work,10*dime,info)
            eval=eval
            evec=evecR
         end if
      else if(kind(eval(1))==16)then
         evec_4=real(Evec,kind=8)
         eval_4=real(Eval,kind=8)
         call DIAG_SYM(eval_4,evec_4,dime,info,ordered)
         evec=real(evec_4,kind=16)
         eval=real(eval_4,kind=16)
      end if
   end subroutine Diag_SYM
   
   subroutine DIAG_C(Eval,Evec,dime,info)
      Integer, INTENT(IN) :: dime
      double complex, INTENT(OUT) :: Eval(dime)
      double complex, INTENT(INOUT) :: Evec(dime,dime)
      Integer, INTENT(OUT) :: info
      double complex Work((7*dime)**2),evalc(dime),evecL(dime,dime),evecR(dime,dime)
      double precision rwork(2*dime)
      
      if(kind(eval(1))==8)then
         call zgeev('N','V',dime,evec,dime,eval,evecL,dime,evecR,dime,work,7*dime,rwork,info)
      else if(kind(eval(1))==4)then
         call cgeev('N','V',dime,evec,dime,eval,evecL,dime,evecR,dime,work,7*dime,rwork,info)
      end if
      eval=eval
      evec=evecR
   end subroutine DIAG_C
   
   subroutine DIAG_4(Eval,Evec,dime,info,CFCE,tol)
      logical, intent(in) :: CFCE !check for complex eigenvalues?
      Integer, INTENT(IN) :: dime
      real(4), intent(in) :: tol
      real(4), INTENT(OUT) :: Eval(dime)
      real(4), INTENT(INOUT) :: Evec(dime,dime)
      Integer, INTENT(OUT) :: info
      integer :: i
      real(4) Work((7*dime)**2),evalc(dime),evecL(dime,dime),evecR(dime,dime)
      real(4) evec_4(dime,dime),eval_4(dime)
      
      info=0
      if(kind(eval(1))==8)then
         call dgeev('N','V',dime,evec,dime,eval,evalc,evecL,dime,evecR,dime,work,7*dime,info)
      else if(kind(eval(1))==4)then
         call sgeev('N','V',dime,evec,dime,eval,evalc,evecL,dime,evecR,dime,work,7*dime,info)
      end if
      if(CFCE)then
         do i = 1,dime
            if(abs(evalc(i))>tol)then
               stop 'DIAG: Complex eigenvalue detected.'
            end if
         end do
      end if
      eval=eval
      evec=evecR
   end subroutine Diag_4
   
   subroutine DIAG(Eval,Evec,dime,info,CFCE,tol)
      logical, intent(in) :: CFCE !check for complex eigenvalues?
      Integer, INTENT(IN) :: dime
      double precision, intent(in) :: tol
      double precision, INTENT(OUT) :: Eval(dime)
      double precision, INTENT(INOUT) :: Evec(dime,dime)
      Integer, INTENT(OUT) :: info
      integer :: i
      double precision Work((7*dime)**2),evalc(dime),evecL(dime,dime),evecR(dime,dime)
      real(4) evec_4(dime,dime),eval_4(dime)
      
      info=0
      if(kind(eval(1))==8)then
         call dgeev('N','V',dime,evec,dime,eval,evalc,evecL,dime,evecR,dime,work,7*dime,info)
      else if(kind(eval(1))==4)then
         call sgeev('N','V',dime,evec,dime,eval,evalc,evecL,dime,evecR,dime,work,7*dime,info)
      else if(kind(eval(1))==16)then
         evec_4=real(Evec,kind=8)
         eval_4=real(Eval,kind=8)
         call DIAG_4(eval_4,evec_4,dime,info,cfce,real(tol,kind=kind(evec_4)))
         evec=real(evec_4,kind=16)
         eval=real(eval_4,kind=16)
         return
      end if
      if(CFCE)then
         do i = 1,dime
            if(abs(evalc(i))>tol)then
               stop 'DIAG: Complex eigenvalue detected.'
            end if
         end do
      end if
      eval=eval
      evec=evecR
   end subroutine Diag
   
   
   subroutine EVAL_TO_ENERGY(eval,dime,A,B,C,J)
      INTEGER dime,i,J
      double precision,intent(inout) :: eval(dime)
      double precision A,B,C
      eval=0.5d0*(A+C)*dble(J*(J+1d0))+0.5d0*(A-C)*eval
   end subroutine EVAL_TO_ENERGY
      
      
   function Reverse_IArr(arr,n)result(res1)
      INTEGER n,i,siz,buf
      Integer, INTENT(IN)::arr(n)
      INTEGER res1(n)
      siz = n/2
      
      do i = 1,n
         res1(i)=arr(i)
      end do
      do i = 1,siz
         buf=res1(i)
         res1(i)=res1(n-i+1)
         res1(n-i+1)=buf
      end do
      
      
   end function Reverse_IArr
      
   !Real To Half-Integer
   function R2HI(num) result(res)   
      double precision,intent(in) :: num
      INTEGER dpIDX
      character(:),allocatable :: res
      character(10) str
      
      if(mod(NINT(num*2),2)==0)then !is not a half-integer
         write(str,'(F5.0)')num
         str=adjustl(str)
         dpIdx=index(str,'.')
         res=trim(str(1:dpidx-1))
      else !is a half-integer
         write(str,'(F5.0)')num*2
         str=adjustl(str)
         dpIdx=index(str,'.')
         res=trim(str(1:dpidx-1))//'/2'
      end if
   end function R2HI
   
   !sigh
   function OldToNewNotation(J,T) result(res)
      INTEGER J,T,rep,i,ii,K1,Km1,Km1Arr(2*J+1),K1Arr(2*J+1)
      logical left
      INTEGER res(3)
            
      res(1)=J
      Km1=J+1
      K1=0
      K1Arr(1)=0
      
      do i=1,2*J
         K1Arr(i+1)=i/2+modulo(i,2)
      end do
      Km1Arr=Reverse_IArr(K1Arr,2*J+1)
      
      do i=1,2*J+1
         if(Km1ARr(i)-K1Arr(i) == T)exit
      end do
      
      res(2)=Km1ARr(i)
      res(3)=K1Arr(i)
      
   end function OldToNewNotation
      
   
   subroutine HAM_DIM_SEQUENCE(Jmin,Jmax,arr)
      INTEGER Jmin,Jmax,j,arr(Jmax-Jmin+1),car
      
      car=1
      do J = Jmin,Jmax
         arr(car)=2*J+1
         car=car+1
      end do
   end subroutine HAM_DIM_SEQUENCE
   
   function RepAxesSwitch_get_undo(rep,left)result(res)
      INTEGER rep,res(3)
      logical left
      if(.not.left)then!right handed representation
         select case(rep)
            case(1)
               !res=[2,3,1]
               res=[3,1,2]
            case(2)
               !res=[3,1,2]
               res=[2,3,1]
            case(3)
               !res=[1,2,3]
               res=[1,2,3]
            case default
               stop "RepAxesSwitch_get wrong rep"
         end select
      else
         select case(rep)
            case(1)
               !res=[3,2,1]
               res=[3,2,1]
            case(2)
               !res=[1,3,2]
               res=[1,3,2]
            case(3)
               !res=[2,1,3]
               res=[2,1,3]
            case default
               stop "RepAxesSwitch_get wrong rep"
         end select
      end if
      
   end function RepAxesSwitch_get_undo
   
   function RepAxesSwitch_get(rep,left) result(res)
      INTEGER rep,res(3)
      logical left
      
      if(.not.left)then!right handed representation
         select case(rep)
            case(1)
               res=[2,3,1]
            case(2)
               res=[3,1,2]
            case(3)
               res=[1,2,3]
            case default
               stop "RepAxesSwitch_get wrong rep"
         end select
      else
         select case(rep)
            case(1)
               res=[3,2,1]
            case(2)
               res=[1,3,2]
            case(3)
               res=[2,1,3]
            case default
               stop "RepAxesSwitch_get wrong rep"
         end select
      end if
   end function RepAxesSwitch_get
   
   function RepAxesSwitch_V(arr,rep,left,undo) result(res)
      double precision arr(3),res(3)
      logical left
      logical,optional :: undo
      INTEGER i,rep,order(3)
      
      
      res=0d0
      if(present(undo))then
         if(undo)then
            order=RepAxesSwitch_get_undo(rep,left)
         else
            order=RepAxesSwitch_get(rep,left)
         end if
      else
         order=RepAxesSwitch_get(rep,left)
      end if
      do i = 1,3
         res(i)=arr(order(i))
      end do   
   end function RepAxesSwitch_V
   
   pure function Permute_matrix(n,m,mat,perm,columns)result(res)
      Integer,intent(in) :: n,m,perm(:)
      INTEGER :: i
      double precision,intent(in) :: mat(n,m)
      double precision res(n,m)
      logical,intent(in) :: columns
   
      res=0d0
      
      if(columns)then
         do i = 1,m
            res(:,i)=mat(:,perm(i))
         end do
      else
         do i = 1,n
            res(i,:)=mat(perm(i),:)
         end do
      end if
   end function Permute_matrix
   
   function RepAxesSwitch_T(mat,rep,left,undo) result(res)
      double precision mat(3,3)
      double precision res(3,3),matbuf(3,3)
      logical left
      logical, optional :: undo
      INTEGER i,j,rep,order(3)
      
      res=0d0
      if(present(undo))then
         if(undo)then
            order=RepAxesSwitch_get_undo(rep,left)
         else
            order=RepAxesSwitch_get(rep,left)
         end if
      else
         order=RepAxesSwitch_get(rep,left)
      end if
      res=Permute_matrix(3,3,mat,order,.true.)
      res=Permute_matrix(3,3,res,order,.false.)
   end function RepAxesSwitch_T
   
   function AxesSwitch_V(arr,order) result(res)
      double precision,intent(in) :: arr(3)
      double precision :: res(3)
      INTEGER i
      Integer,intent(in) :: order(3)
      
      
      res=0d0
      do i = 1,3
         res(i)=arr(order(i))
      end do   
   end function AxesSwitch_V

   function AxesSwitch_T(mat,order)result(res)
      double precision,intent(in) :: mat(3,3)
      Integer,intent(in) :: order(3)
      double precision :: res(3,3)
      
      res=0d0
      res=Permute_matrix(3,3,mat,order,.true.)
      res=Permute_matrix(3,3,res,order,.false.)
   end function AxesSwitch_T
   
   function REPRESENTATION(RAY) result (res)
      double precision RAY,d,e,f
      INTEGER res
      
      if((RAY>1.0d0).OR.(RAY<-1.0d0)) STOP 1
      
      d = DABS(RAY-1.0)
      e = DABS(RAY-0.0)
      f = DABS(RAY+1.0)
      res=0
      if((d<e).AND.(d<f)) then
         res = 3
      else if((e<d).AND.(e<f)) then
         res = 2
      else if((f<e).AND.(f<d)) then
         res = 1
      end if
   end function REPRESENTATION
   
   !returns the dimensions of the full Hamiltonian for multiple Js in the range of [Jmin]-[Jmax]
   function DIM_HAMILTONIAN(Jmin,Jmax) result (res)
      implicit none
      INTEGER Jmin,Jmax,res,i
      
      res = 0
      do i=Jmin,Jmax
         res = res + (i*2+1)
      end do
   end function
   
   !imprints (copies) 1D array [arr2] with size [m] onto 1D array [arr1] with size [n], [arr1] has to be larger or the same size as [arr2],
   ![start] marks the starting index on the larger array [arr1]
   !eg. (0 0 0 0 0 0) <- (4 5 7), starting index 2 => (0 4 5 7 0 0) 
   subroutine COPY_ARRAY1D(arr1,n,arr2,m,start)
      implicit none
      Integer,intent(in) :: n,m,start
      double precision,intent(inout) :: arr1(n)
      double precision,intent(in) :: arr2(m)
      
      if(n<m) stop "Copied-to array smaller than copied array"
      
      arr1(start:(start+m-1))=arr2
   end subroutine COPY_ARRAY1D
   
   !imprints (copies) 2D array [arr2] with sizes [m1,m2] onto 2D array [arr1] with sizes [n1,n2], [arr1] has to be larger or the same size as [arr2] in any dimension,
   ![start1,start2] mark the starting indices on the larger array [arr1] 
   subroutine COPY_ARRAY2D(arr1,n1,n2,arr2,m1,m2,start1,start2)
      implicit none
      Integer,intent(in) :: n1,n2,m1,m2,start1,start2
      double precision, INTENT(INOUT) :: arr1(n1,n2)
      double precision, INTENT(in) :: arr2(m1,m2)
      arr1(start1:(start1+m1-1),start2:(start2+m2-1))=arr2
   end subroutine COPY_ARRAY2D
   
   function ThirdNum(a,b)result(c)
      integer a,b,c
      if(a==2 .and. b==2)then
         c=0
      else
         c=6-a-b
      end if
   end function ThirdNum
   
   
   function HAM_SR_OWN_MATEL_NKSJM(eps,N2,K2,S2,J2,M2,N1,K1,S1,J1,M1,space)result(res)
      integer :: N2,K2,N1,K1,N3,K3,a,b,c,j3_idx,m3_idx,j_mult,m_mult,Mn1,Mn2,Ms2_idx,Ms1_idx,F,mult_s
      logical space
      double precision :: S2,S1,J2,J1,M2,M1,J3,M3,S3,Ms2,Ms1
      double precision :: eps(3,3),LeviC(0:4,0:4,0:4)
      double complex :: buf,res,alt,bufAlt,bufalt2
      character(1),parameter :: axes(0:4) = ['N','X','Y','Z','N']
      
      res=0d0
      alt=0d0
      LeviC=Levi()
      mult_s=NINT(2*S2)
      if(S2/=S1)return
      if(space)then
         j_mult=NINT(2*max(J1,J2))
         s3=s2
         do a = 1,3
            do b = 1,3
               c=ThirdNum(a,b)
               buf=0d0
               ! goto 0451
               ! do k3 = k2-1,k2+1
                  ! j3=j1
                  ! m3=m1
                  ! do j3_idx = 0,j_mult
                     ! j3=j3_idx-j_mult
                     ! do m3_idx = 0,NINT(2*J3)
                        ! m3=m3_idx-J3
                        ! n3=n2
                        ! buf=buf+Nalpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(a))*Salpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(b),.true.)
                        
                        ! n3=n1
                        ! buf=buf+Salpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(b),.true.)*Nalpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(a))
                     ! end do
                  ! end do
               ! end do
               
! 0451           continue
               bufalt=0d0
               do Ms2_idx=0,mult_S
                  ms2=ms2_idx-S2
                  mn2=NINT(M2-Ms2)
                  do ms1_idx=0,mult_S
                     ms1=ms1_idx-S1
                     Mn1=NINT(M1-Ms1)
                     bufalt2=0d0
                     do F = 1,3
                        bufalt2=bufalt2+(2*DIR_COSE_N_FROM_SPH(N2,N1)*DIR_COSE_NM_FROM_SPH(N2,Mn2,N1,Mn1,axes(F))*DIR_COSE_NK_N_matel(N2,K2,N1,K1,axes(b),axes(a))-iu*LeviC(a,b,c)*DIR_COSE_NKM_FROM_SPH(N2,K2,Mn2,N1,K1,Mn1,axes(F),axes(c)))*SAlpha_matel_space(S2,Ms2,S1,Ms1,axes(F))
                     end do
                     bufalt=bufalt+bufalt2*CG_SR(N2,S2,J2,Mn2,Ms2,M2)*CG_SR(N1,S1,J1,Mn1,Ms1,M1)
                  end do
               end do
               
               !res=res+eps(a,b)*buf
               alt=alt+eps(a,b)*bufAlt
            end do 
         end do
         
      else
         if(K2==K1)then
            do a = 1,3
               do b = 1,3
                  buf=0d0
                  do K3 = max(K2-1,-N2),min(K2+1,N2)
                     N3=N2
                     J3=J1
                     M3=M1
                     S3=S2
                     buf=buf + Nalpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(a))*Salpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(b),.false.)
                     
                     N3=N1
                     J3=J2
                     M3=M2
                     S3=S1
                     buf=buf + Salpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(b),.false.)*Nalpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(a))
                  end do
                  res = res + eps(a,b)*buf
               end do
            end do
         else if(abs(K2-K1)==1)then
            do a = 1,3
               do b = 1,3
                  buf=0d0
                  do K3 = min(K2,K1),max(K2,K1)
                     N3=N2
                     J3=J1
                     M3=M1
                     S3=S2
                     buf=buf + Nalpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(a))*Salpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(b),.false.)
                     
                     N3=N1
                     J3=J2
                     M3=M2
                     S3=S1
                     buf=buf + Salpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(b),.false.)*Nalpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(a))
                  end do
                  res = res + eps(a,b)*buf
               end do
            end do
         else if(abs(K2-K1)==2)then
            K3=max(k2,k1)-1
            do a = 1,3
               do b = 1,3
                  buf=0d0
                  N3=N2
                  J3=J1
                  M3=M1
                  S3=S2
                  buf=buf + Nalpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(a))*Salpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(b),.false.)
                  
                  N3=N1
                  J3=J2
                  M3=M2
                  S3=S1
                  buf=buf + Salpha_matel_NKSJM(N2,K2,S2,J2,M2,N3,K3,S3,J3,M3,axes(b),.false.)*Nalpha_matel_NKSJM(N3,K3,S3,J3,M3,N1,K1,S1,J1,M1,axes(a))
                  res = res + eps(a,b)*buf
               end do
            end do
         end if
      end if
      
      !alt=alt*0.5d0
      res=alt*0.5d0
   end function HAM_SR_OWN_MATEL_NKSJM

   function HAM_RR_OWN_MATEL_NKM(abc,N2,K2,Mn2,N1,K1,Mn1)result(res)
      integer :: N2,K2,Mn2,N1,K1,Mn1
      double precision :: res,abc(3)
      
      res=0d0
	   if(N2/=N1 .or. (K2/=K1 .and. K2-2/=K1 .and. K2+2 /= K1) .or. Mn2 /= Mn1)return
	   
      if(K2==K1-2)then
         res=(abc(1)-abc(2))*(0.25d0*sqrt(dble(N1*(N1+1)-(K1-1)*(K1-2)))*sqrt(dble(N1*(N1+1)-K1*(K1-1))))
      else if(K2==K1+2)then
         res=(Abc(1)-abc(2))*(0.25d0*sqrt(dble(n1*(n1+1)-(K1+1)*(K1+2)))*sqrt(dble(N1*(N1+1)-K1*(K1+1))))
      else if(K2==K1)then
         res=(abc(1)+abc(2))*0.5d0*(N1*(N1+1)-K1**2)+abc(3)*K1**2
      end if
   end function HAM_RR_OWN_MATEL_NKM

	function HAM_RR_OWN_MATEL_NKSJM(abc,N2,K2,S2,J2,M2,N1,K1,S1,J1,M1)result(res)
      integer :: N2,K2,N1,K1,mn2,mn1,ms2_idx,ms1_idx
      double precision :: S2,S1,J2,J1,M2,M1,ms2,ms1
      double precision :: abc(3)
      double complex :: res
      
      res=0d0
      if(S2/=S1 .or. m2/=m1 .or. J2/=J1)return
      ! do ms2_idx=0,NINT(S2*2)
         ! ms2=ms2_idx-S2
         ! mn2=NINT(M2-ms2)
         ! do ms1_idx=0,NINT(S1*2)
            ! ms1=ms2
            ! mn1=NINT(M1-Ms1)
            ! res=res + CG_SR(N2,S2,J2,mn2,ms2,m2)*CG_SR(N1,S1,J1,mn1,ms1,m1)*HAM_RR_OWN_MATEL_NKM(abc,N2,K2,MN2,N1,K1,Mn1)
         ! end do
      ! end do
      res=HAM_RR_OWN_MATEL_NKM(abc,N2,K2,0,N1,K1,0)
   end function HAM_RR_OWN_MATEL_NKSJM

   function HAM_RR_OWN(N2,N1,A,B,C)result(HAM)
      integer N2,N1,k2,k1
      double precision HAM(2*N2+1,2*N1+1),A,B,C
      
      if(N2/=N1)then
         Ham=0d0
         return
      end if
      !left bra
      do k2 = -N2,N2
         !right ket
         do k1 = -N1,N1
            if(K2==K1-2)then
               HAM(k2+n2+1,k1+n1+1)=(A-B)*(0.25d0*sqrt(dble(N1*(N1+1)-(K1-1)*(K1-2)))*sqrt(dble(N1*(N1+1)-K1*(K1-1))))
            else if(K2==K1+2)then
               HAM(k2+n2+1,k1+n1+1)=(A-B)*(0.25d0*sqrt(dble(n1*(n1+1)-(K1+1)*(K1+2)))*sqrt(dble(N1*(N1+1)-K1*(K1+1))))
            else if(K2==K1)then
               HAM(k2+n2+1,k1+n1+1)=(A+B)*0.5d0*(N1*(N1+1)-K1**2)+C*K1**2
            else
               HAM(k2+n2+1,k1+n1+1)=0d0
            end if
           ! HAM(k2+n2+1,k1+n1+1)=(A-B)*( 0.25d0*KD(K2,K1-2)*sqrt(dble(N1*(N1+1)-(K1-1)*(K1-2)))*sqrt(dble(N1*(N1+1)-K1*(K1-1)))&
           ! + 0.25d0*KD(K2,K1+2)*sqrt(dble(n1*(n1+1)-(K1+1)*(K1+2)))*sqrt(dble(N1*(N1+1)-K1*(K1+1)))&
           ! + 0.5d0*KD(K2,K1)*(N1*(N1+1)-K1*(K1+1)))+C*K1**2*KD(K2,K1)
         end do
      end do
   end function HAM_RR_OWN

   function KD(a,b)result(res)
      integer a,b,res
      if(a==b)then
         res=1
      else
         res=0
      end if
   end function KD
   
   !gets the full block-diagonal assymetric rotor Hamiltonian [HamFull] with sizes [N,N] for Js in range [Jmin]-[Jmax]
   subroutine FULL_HAMILTONIAN_ASS(A,B,C,Jmin,Jmax,rep,HamFull,N,left)
      implicit none
      double precision, intent(in) :: A,B,C
      Integer, intent(in) :: Jmin,Jmax,N,rep
      double precision, intent(out) :: HamFull(N,N)
      double precision, allocatable :: Ham(:,:)
      INTEGER i,j,k,car,len1(2)
      INTEGER len2(2),curdim
      logical left
      
      len1 = shape(HamFull)
      car = 1
      do i=Jmin,Jmax,1
         allocate(Ham(i*2+1,i*2+1))
         curdim = 2*i+1
         
         len2=shape(Ham)
         do j = 1,len2(1),1
            do k = 1,len2(1),1
               Ham(j,k)=0.0d0
            end do
         end do
         call ASS_ROTOR([A,B,C],i,rep,Ham,curdim,left)
         !call W2dArr(Ham,curdim,curdim)
         call COPY_ARRAY2D(HamFull,len1(1),len1(2),Ham,len2(1),len2(2),car,car)
         car = car + curdim
         deallocate(Ham)
      end do
      
   end subroutine FULL_HAMILTONIAN_ASS
   
   !gets the assymetric rotor Hamiltonian [Ham] with sizes [N,N] for one specific [J], representation [rep] in range -1 to 1
   subroutine ASS_ROTOR(abc,J,rep,Ham,N,left)
      implicit none      
      double precision, intent(in) :: abc(3)
      double precision buf,kappa,Fv,Gv,Hv,A,B,C
      Integer,intent(in) :: J,rep,N
      INTEGER k1,k2,idx,idx2,idx3
      double precision,intent(out) :: Ham(N,N)
      logical left,mask1(3)
      
      A=abc(1)
      B=abc(2)
      C=abc(3)
      !sort abc if not given like A>B>C
      if(.not.A>=B .or. .not.B>=C)then
         idx=maxloc([A,B,C],dim=1)
         mask1=.true.
         mask1(idx)=.false.
         idx2=maxloc([A,B,C],dim=1,mask=mask1)
         idx3=6-idx-idx2
         A=abc(idx)
         B=abc(idx2)
         C=abc(idx3)
      end if
      
      kappa = (2.0d0*B-A-C)/(A-C)
      
      !dependent on the choice of [rep]
      Fv = Fval(kappa,rep)
      Gv = Gval(kappa,rep) 
      Hv = Hval(kappa,rep) 
      
      if(left)then
         Hv=-Hv
      end if
      
      ! Ham=0.0d0
      
      do k1=-J,J
         do k2=-J,J
            Ham(k1+J+1,k2+J+1)=MATEL_H_RIG(J,k1,k2,Fv,Gv,Hv)
            ! if(k1==k2)then
               ! Ham(k1,k2)=Fv*J*(J+1.0D0)+(Gv-Fv)*(k1-J-1.0D0)**2.0D0
            ! end if
            ! continue
            ! if((k1+2)==k2)then
               ! buf=Hv*sqrt(fun(J,k1-J-1))
               ! Ham(k1,k2)=buf
               ! Ham(k2,k1)=buf
            ! end if
         end do
      end do
   end subroutine ASS_ROTOR
   
   function MATEL_H_RIG(J,K1,K2,F,G,H)result(res)
      INTEGER J,K1,K2
      double precision res,F,G,H
      res=0d0
      if(K1==K2)then
         res=F*dble(J*(J+1)-K1**2)+G*dble(K1**2)
      else if(K1==K2+2 .or. K1==K2-2)then
         res=H*sqrt(Fun_new(J,min(K1,K2)+1))
      end if
   end function MATEL_H_RIG
   
   subroutine LIN_ROTOR()
      
   end subroutine LIN_ROTOR
   
   subroutine SYM_ROTOR()
   
   end subroutine SYM_ROTOR
   
   subroutine SPH_ROTOR()
      
   end subroutine SPH_ROTOR
   
   function Ray(A,B,C) result (res)
      double precision A,B,C,res
      res = (2.0d0*B-A-C)/(A-C)
   end function Ray
   
   function Fun(J,K) result (res)
      INTEGER J,K
      double precision res
      res = 0.25D0*((J*(J+1)-K*(K+1))*(J*(J+1)-(K+2)*(K+1)))
   end function Fun
   
   function Fun_new(J,K) result (res)
      INTEGER J,K
      double precision res
      ! res = 0.25d0*dble((J*(J+1)-K*(K+minus))*(J*(J+1)-(K+minus)*(K+minus*2)))
      !res = 0.25d0*dble((J*(J+1)-K*(K+minus))*(J*(J+1)-(K+minus)*(K+minus*2)))
      K=abs(K)
      res = 0.25d0*dble((J*(J+1)-K*(K+1))*(J*(J+1)-K*(K-1)))
   end function Fun_new
      
      
   !https://doi.org/10.1103/RevModPhys.23.213   
   function Fun_vanVleck(J,K) result(res)
      INTEGER J,K
      double precision res
      
      res=sqrt(dble(J*(J+1)-K*K+1))
   end function Fun_vanVleck
      
   function Fval(kappa,n) result (res)
      implicit none
      double precision kappa,res
      INTEGER n
      if(n==1) then
         res = 0.5d0*(kappa-1.0d0)
         return
      else if(n==3) then
         res = 0.5d0*(kappa+1.0d0)
         return
      else if(n==2) then
         res = 0
         return
      else
         stop "Incorrect ass. top representation"
      end if
   end function Fval
   
   function Gval(kappa,n) result(res)
      implicit none
      double precision kappa,res
      INTEGER n
      if(n==1) then
         res = 1d0
         return
      else if(n==3) then
         res = -1d0
         return
      else if(n==2) then
         res = kappa
         return
      else
         stop "Incorrect ass. top representation"
      end if
   end function Gval
   
   function Hval(kappa,n) result(res)
      implicit none
      double precision kappa,res
      INTEGER n
      if(n==1) then
         res = -0.5d0*(kappa+1.0d0)
         return
      else if(n==3) then
         res = 0.5d0*(kappa-1.0d0)
         return
      else if(n==2) then
         res = 1
         return
      else
         stop "Incorrect ass. top representation"
      end if
   end function Hval
   
!  write a very fancy 2D array, [n,m] dimensions from SHAPE function
!  useful, if u dont know GDB
   subroutine W2dArr(arr,n,m)
      INTEGER n,m
      double precision, dimension(n,m) :: arr
      Integer, dimension(2) :: leng
      INTEGER i,j
      leng = shape(arr)

      write(*,"(1A8)", advance = "no") "        "
      do i = 1,leng(2)
         write(6,"(1I4)", advance = "no") i
         write(6,"(1A4)", advance = "no") "     "
      end do

      write (*,"(1A1)") " "
      i = 1

      do i = 1,leng(1)
         write(6,"(1I4)", advance = "no") i
         do j = 1,leng(2)
            write (6,"(1F8.3)", advance = "no") arr(j,i)
         end do
         write(6,"(1A1)") " "
      end do
   end subroutine W2darr

   subroutine W2dArr_log(arr,n,m)
      INTEGER n,m
      logical, dimension(n,m) :: arr
      Integer, dimension(2) :: leng
      INTEGER i,j
      leng = shape(arr)

      write(*,"(1A6)", advance = "no") "      "
      do i = 1,leng(2)
         write(6,"(1I2)", advance = "no") i
         !write(6,"(1A4)", advance = "no") "     "
      end do

      write (*,"(1A1)") " "
      i = 1

      do i = 1,leng(1)
         write(6,"(1I4)", advance = "no") i
         do j = 1,leng(2)
            write (6,"(1L2)", advance = "no") arr(j,i)
         end do
         write(6,"(1A1)") " "
      end do
   end subroutine W2darr_log

!  write a very fancy 2D array, [n,m] dimensions from SHAPE function, can set starting indices when printed
   subroutine W2dArr_index(arr,n,m,start1,start2)
      INTEGER n,m
      double precision, dimension(n,m) :: arr
      Integer, dimension(2) :: leng
      INTEGER start1,start2
      INTEGER i,j
      leng = shape(arr)

      write(*,"(1A8)", advance = "no") "        "
      do i = 1,leng(2)
         write(6,"(1I4)", advance = "no") (start1 + i - 1)
         write(6,"(1A4)", advance = "no") "     "
      end do

      write (*,"(1A1)") " "
      i = 1

      do i = 1,leng(1)
         write(6,"(1I4)", advance = "no") (start2 + i - 1)
         do j = 1,leng(2)
            write (6,"(1F8.3)", advance = "no") arr(j,i)
         end do
         write(6,"(1A1)") " "
      end do
   end subroutine W2darr_index
   
   
end module rotor

module input
   use constants   
   use rotor
   use arrayOperations
   use molecule
   
   implicit none
   INTEGER n0,MENDELEV
   parameter (n0=100,MENDELEV=89)
   
   double precision quad(3,3),quadtr(3,3),qtot(3,3),g(3,3),gn(3,3),ge(3,3),gt(3,3),gnew(3,3),g_el(3,3)
   double precision xp(3,3),ama(n0),eps(3,3)
   double precision dip_el(3)
   double precision abc(3),mass,com_shift(3),charge
   double precision,target :: dip(3),dip_nuc(3)
   double precision,allocatable :: a_masses(:),r(:,:),r0(:,:),r1(:,:)
   logical dip_found, g_found, q_found, ABC_found,r_Found,eps_found
   
   INTEGER :: unitt,nAtoms,line
   character(150) :: FN
   character(6),parameter :: FNForm='(A150)'
   character(2),allocatable :: asy(:)
   
   public :: quad,quadtr,qtot,g,gn,ge,gt,g_el,gnew,xp,r,ama,dip,abc,mass,dip_found, g_found, q_found, ABC_found,r_Found,eps_found,eps,dip_nuc,dip_el
   private :: FN,unitt,line
   public :: nAtoms,com_shift,r0,r1,charge
   
   contains
   
   
   
   !What I've learned when writing fortran:
   !Any attempt at modularity will add double the time required on a functional solution
   !There are no docs
   !Strings are a pain
   !Have fun hunting SEGFAULTS which appear at -O3
   
   
   
      
   !write double 2D array to sequential file
   subroutine ARR2D_2_FILE(arr,dime,unitt,line)
      INTEGER dime,unitt,line,i,j
      double precision arr(dime,dime)
      character(3)str
      character(15)str2
      
      rewind(unitt)
      write(str,*)dime
      !FORMAT
      str2 = (str//'(1F15.8,A1))')   
      do i = 1,line-1
         read(unitt,*)
      end do
      do i = 1,dime
         do j = 1,dime
            write(unitt,str2)arr(i,j)
         end do
      end do
   end subroutine ARR2D_2_FILE
 
   subroutine ARR2D_2_FILE_spec(arr,n,m,maxidx,unitt)
      INTEGER unitt,i,j,maxidx,n,m
      double precision arr(n,m)
      character(3)str
      character(15)str2
      
      !rewind(unitt)
      str='   '
      write(str,'(I3)')m
      !FORMAT
      str2 = ("("//trim(adjustl(str))//"(E25.15E4,' '))")   
      !do i = 1,line-1
      !   read(unitt,*)
      !end do
      do i = 1,maxidx-1
         write(unitt,str2)(arr(i,j),j=1,m)
      end do
   end subroutine ARR2D_2_FILE_spec
 
   subroutine RotateVector(a,b,c,v)
      double precision,intent(in) :: a,b,c
      double precision,intent(inout) :: v(3)
      
      v = MATMUL(RotationMatrix(a,b,c),v)
   end subroutine RotateVector
   
   subroutine RotateMatrix(a,b,c,v)
      double precision,intent(in) :: a,b,c
      double precision,intent(inout) :: v(3,3)
      
      v = MATMUL(RotationMatrix(a,b,c),v)
   end subroutine RotateMatrix
   
   function RotationMatrix(a,b,c) result(mat)
      double precision,intent(in) :: a,b,c
      double precision mat(3,3)
      
      mat(1,1)=cos(a)*cos(b)
      mat(2,1)=sin(a)*cos(b)
      mat(3,1)=-sin(b)
      mat(1,2)=cos(a)*sin(b)*sin(c)-sin(a)*cos(c)
      mat(2,2)=sin(a)*sin(b)*sin(c)+cos(a)*cos(c)
      mat(3,2)=cos(b)*sin(c)
      mat(1,3)=cos(a)*sin(b)*cos(c)+sin(a)*sin(c)
      mat(2,3)=sin(a)*sin(b)*cos(c)-cos(a)*sin(c)
      mat(3,3)=cos(b)*cos(c)
   end function RotationMatrix
   
   function Shift_Coordinate_dip(dip,transform)result(res)
      double precision dip(3),transform(3),res(3)
      INTEGER i
      
      do i = 1,3
         res(i)=dip(i)*transform(i)
      end do
   end function Shift_Coordinate_dip
   
   function Shift_Coordinate_quad(quadr,dipo,transform) result(res)
      double precision quadr(3,3),dipo(3),transform(3),res(3,3)
      INTEGER i,j
      
      do i = 1,3
         do j = 1,3
            res(i,j)=quadr(i,j)+0.5*(transform(i)*dipo(j)+dipo(i)*transform(j))
         end do
      end do
   end function Shift_Coordinate_quad
   
   
   function Str2Num(str) result(num)
      character(*) str
      double precision num(3)
      INTEGER i
      
      read(str,*)num
   end function Str2Num
   
   
   
   !gets a more accurate mass than the one from gaussian
   function getMass(n,m)result(mass)
      INTEGER :: n,atomic,nucleon,length,atomic2,i
      double precision :: m,mass
      double precision,allocatable :: masses(:),masses_copy(:)
      character(4) :: symbol
      open(99,file='isotopes.par')
!19       continue
19       read(99,1996)atomic,symbol,nucleon,mass
         length=0
         if(atomic==n)then
            length=1
            read(99,1996)atomic2,symbol,nucleon,mass
            do while(atomic2==atomic)
               length=length+1
               read(99,1996)atomic2,symbol,nucleon,mass
            end do
            
            do i = 1,length+1
               backspace(99)
            end do
            
            allocate(masses(length),masses_copy(length))
            do i=1,length
               read(99,1996)atomic2,symbol,nucleon,masses(i)
            end do
            masses_copy=masses
            masses=abs(masses-m)
            mass = masses_copy(minloc(masses,1))
            goto 20
         else
            goto 19
         end if
         
!21       if(atomic==n .and.abs(mass/m-1d0)<1d-1)goto 20
      !goto 19
20    close(99)
      deallocate(masses,masses_copy)
1996  format(BN,I4,A4,I5,F15.11)
      !BN = left justified
      !/ = dont read the next record
   end function getMass
   
   subroutine Reorder_Tensor(arr,order)
      double precision,intent(inout) :: arr(3,3)
      double precision bufarr(3),bufarr2(3,3)
      INTEGER first,second,third,i
      Integer,intent(out) :: order(3)
      
      do i = 1,3
         bufarr(i)=abs(arr(i,i))
      end do
      third = maxloc(bufarr,dim=1)
      first = minloc(bufarr,dim=1)
      
      do i = 1,3
         if(i/=first .and. i/=third)then
            second=i
            exit
         end if
      end do
      
      order=[first,second,third] !bruh
      bufarr2=arr
      arr = axesSwitch_T(arr,order)
   end subroutine Reorder_Tensor
   
   subroutine Backward(unit,lines)
      integer(4) unit,lines
      integer(4) i
      
      do i = 1,lines
         backspace(unit)
      end do
   end subroutine Backward
   
   subroutine Forward(unitt,howmany)
      INTEGER unitt,howmany,i
      character(80) str
      
      do i = 1,howmany
         read(unitt,*)str
      end do
      line=line+i-1
   end subroutine Forward
   
   function GetJobEnds()result(jobEnds)
      INTEGER :: i,jobCount,idx
      Integer,allocatable :: jobEnds(:),jobEndsHoldover(:)
      character(150) FN
      
      allocate(jobEnds(100))
      jobEnds=0
      rewind(unitt)
      line=0
      jobCount=0
      do while(.true.)
         line=line+1
         read(unitt,2000,end=2001)FN
         idx=INDEX(FN,'Normal termination of Gaussian')
         if(idx>0)then
            jobCount=jobCount+1
            jobEnds(jobcount)=line
         end if
      end do
      
2001  jobEndsHoldover=jobEnds
      
      deallocate(jobEnds)
      
      allocate(jobEnds(jobCount))
      
      do i = 1,jobCount
         jobEnds(i)=jobEndsHoldover(i)
      end do
      
      deallocate(jobEndsHoldover)
2000  format(A150)   
      
      rewind(unitt)
      line=0
   end function GetJobEnds
   
   subroutine GoToOpts(n,jobEnds,i)
      INTEGER :: n,jobEnds(n),i,idx,idx2
      logical :: switch
      
      rewind(unitt)
      line=0
      if(i>1)then
         call Forward(unitt,jobEnds(i-1))
      end if
      
      !line=0
      do while(.true.)
         line=line+1
         read(unitt,FNForm,end=2002)FN
         idx = Index(FN,'******************************************')
         if(idx>0)then
            read(unitt,FNForm,end=2002)FN
            idx2 = Index(FN,'Gaussian 16:')
            read(unitt,FNForm,end=2002)FN
            read(unitt,FNForm,end=2002)FN
            line=line+3
            idx = Index(FN,'******************************************')
            if(idx==0 .or. idx2==0)cycle
            
            switch=.false.
            do while(.true.)
               
               read(unitt,FNForm,end=2002)FN
               line=line+1
               idx = Index(FN,'------------')
               if(idx/=0 .and. switch)then
                  return
               else if(idx/=0 .and. .not.switch)then
                  switch=.true.
               end if
            end do
         end if
      end do
2002  return
   end subroutine GoToOpts
   
   !untested
   function StrContains(str,n,charArr)result(res)
      INTEGER n,i
      character(*) str,charArr(n)
      logical res
      
      res=.false.
      do i = 1,n
         if(index(str,charArr(i))>0)then
            res=.true.
            return
         end if
      end do
      
   end function StrContains
   
   function findOpt(left,right)result(valuee)
      INTEGER valuee,i,left,right,endIdx,valueIdx,eqIdx,lines
      character(20) leftstr,rightstr,holder
      character(:), allocatable :: leftstrA, rightstrA   
            
      
      write(leftstr,'(I10)')left
      leftstrA=trim(adjustL(leftstr))
      leftstrA=leftstrA//'/'
      write(rightstr,'(I10)')right
      rightstrA=trim(adjustL(rightstr))
      rightstrA=rightstrA//'='
      
      lines=1
      read(unitt,FNForm)FN
      do while(index(FN,'----')==0)
         if(index(FN,leftStrA)>0)then
            valueIdx=index(FN,rightStrA)
            if(valueIdx>0)then
               eqIdx =index(FN(valueIdx:150),'=')+valueIdx-1
               
               i=eqIdx+1
               do while(FN(i:i)/='/' .and. FN(i:i)/=';' .and. FN(i:i)/=',')
                  i=i+1
               end do
               endIdx=i
               
               read(FN(eqIdx+1:endIdx),*)valuee
               return
            end if
         end if
         lines=lines+1
         read(unitt,FNForm)FN
      end do
      
      do i = 1,lines
         backspace(unitt)
      end do
      
   end function findOpt
   
   function getNumOfAtoms()result(res)
      INTEGER res
      
      rewind(unitt)
      do while(.true.)
         read(unitt,FNForm)FN
         if(index(FN,'NAtoms=')>0)then
            read(FN(9:13),'(I5)')res
            exit
         end if
      end do
   end function getNumOfAtoms
   
   
   
   subroutine Select_sort(arr,n,order)
      Integer,intent(in) :: n
      double precision,intent(inout) :: arr(n)
      Integer,intent(out),optional :: order(n)
      
      double precision :: bufi
      INTEGER i,j,idx
      
      if(present(order))then
         do i = 1,n
            order(i)=i
         end do
         
         do i = 1,n
            bufI=arr(i)
            idx=i
            do j = i+1,n
               if(arr(j)<bufI)then
                  bufI=arr(j)
                  idx=j
               end if
            end do
            call SwapEl(arr,n,i,idx)
            call SwapEl_int(order,n,i,idx)
         end do
      else
         do i = 1,n
            bufI=arr(i)
            idx=i
            do j = i+1,n
               if(arr(j)<bufI)then
                  bufI=arr(j)
                  idx=j
               end if
            end do
            call SwapEl(arr,n,i,idx)
         end do
      end if
   end subroutine Select_sort
   
   !iz a mess, dont judge
   subroutine readFile(filename,outputfilen,flags,sta,separateElNuc,rep,lh)
      INTEGER nat,ia,ierr,sta
      logical lex,outp,flags(10),rot_found,dipinpax,separateElNuc
      CHARACTER(2) atsy(MENDELEV)
      INTEGER lines,i,j,k,m,ka,iz(n0),l,arr(2),ii,n_atoms,ng
      double precision bohr,mp,c,u,amu,curmass
      double precision q(3,3),gq(3,3),tvec(3),dipbuf(3),quadbuf(3,3),qt(3,3),gdalton(3,3)
      double precision mt,cm(3),arrbuf(3,3),bufr1(3),bufr2(3),t0(3,3),t1(3,3),rotmat(3,3)
      double precision a(3),t(3,3),z(3,3),e(3,3)
      double precision,dimension(:),pointer :: dip_ptr
      character(80) filename
      character(:),allocatable :: outputfilen
      INTEGER idx,idx2,order(3)
      INTEGER rep
      logical lh
      
      !get number of gaussian jobs and their endings
      INTEGER jobCount,jobIndex,jobEnd,jobstart,lineStart
      !INTEGER line
      Integer,allocatable :: jobEnds(:)
      
      !dipole moment
      logical :: isNuclear
      
      
      
      data atsy/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
      'Na','Mg','Al','Si',' P',' S','Cl','Ar', &
      ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
               'Ga','Ge','As','Se','Br','Kr',&
      'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',&
               'In','Sn','Sb','Te',' I','Xe',&
      'Cs','Ba','La',&
                    'Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho',&
                    'Er','Tm','Yb','Lu',&
      'Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',&
                'Tl','Pb','Bi','Po','At','Rn',&
      'Fr','Ra','Ac'/
     
      double precision  amas(MENDELEV)
      data amas/1.008,4.003,&
      6.941,9.012,10.810,12.011,14.007,15.999,18.998,20.179,&
      22.990,24.305,26.981,28.086,30.974,32.060,35.453,39.948,&
      39.098,40.080,44.956,47.900,50.941,51.996,54.938,55.847,&
      58.933,58.700,63.546,65.380,&
      69.720,72.590,74.922,78.960,79.904,83.800,&
      85.468,87.620,88.906,91.220,92.906,95.940,98.906,101.070,&
      102.906,106.400,107.868,112.410,&
      114.82,118.69,121.75,127.600,126.905,131.300,&
      132.905,137.330,138.906,&
      140.120,140.908,144.240,145.000,150.400,&
      151.960,157.250,158.925,162.500,164.930,167.260,168.934,&
      173.040,174.970,&
      178.490,180.948,183.850,186.207,190.200,192.220,195.090,&
      196.967,207.590,204.370,207.200,208.980,210.000,210.001,&
      222.02,&
      223.000,226.025,227.028/  

!     conversions
      bohr=0.52917721092d0
      mp=1.00727646688d0
      c=2.99792458d10
      u=219474.63068d0
      amu=1822.88839d0
     
     
!     File Check
      unitt=2
      line=0
      inquire(file=filename,exist=lex)
      if(.not.lex)then
         sta=-1
         return
      end if
      sta=0
      OPEN(unitt,FILE=filename,ERR=1111)
            
      
      
2000  format(A150)
      
      jobEnds=GetJobEnds()
      jobCount=size(jobEnds,dim=1)
      nAtoms=getNumOfAtoms()
      allocate(r(3,nAtoms))
      allocate(a_masses(nAtoms))
      
2002  iz=0
      ama=0d0
      nat=0
      l=0
      !line=0
      rot_found=.false.
      dipinpax=.false.
!     Reads every line of the gaussian output file and checks every line for wanted parameters

      rewind(unitt)
      do i = 1,jobCount
         jobEnd=jobEnds(i)
         
         if(i==1)then
            jobStart = 1
         else
            jobStart = jobEnds(i-1)+1
         end if
         
         ! if(separateElNuc)then
            ! call GoToOpts(jobcount,jobEnds,i)
            ! isNuclear=findOpt(6,120)>0
         ! else
            ! isNuclear=.false.
         ! end if
         ! if(isNuclear)then
            ! dip_ptr=>dip_nuc
         ! else
            ! dip_ptr=>dip
         ! end if
         
         !lineStart=line
         !line=lineStart
         
         dip_ptr=>dip
         do while(line<=jobEnd)
            line=line+1
            read(unitt,2000,end=50)FN
            if(FN(1:12)==' Atom  1 has')then
               do j = 1,natoms
                  read(FN(36:44),*)a_masses(j)
                  read(FN(28:29),'(I2)')iz(j)
                  a_masses(j)=GetMass(iz(j),a_masses(j)) !???
                  read(unitt,2000)FN
                  line=line+1
               end do
            end if
            
            if(FN(27:44)=='Input orientation:')THEN
               read(unitt,2000)FN
               read(unitt,2000)FN
               read(unitt,2000)FN
               read(unitt,2000)FN
               line=line+4
               do j = 1,natoms
                  READ(UNITT,5001)IA,KA,IA,r(1,j),r(2,j),r(3,j)
               end do
5001           format(1x,I6,6x,I6,6x,I6,3x,3(1x,F11.6))
            ENDIF
            
            if(FN(1:44)==' electronic spin - molecular rotation tensor')then
               read(unitt,2000)FN
5002           format(3(6X,F11.4))   
               read(FN(1:51),5002)(eps(1,j),j=1,3)
               read(unitt,2000)FN
               read(FN(1:51),5002)(eps(2,j),j=1,3)
               read(unitt,2000)FN
               read(FN(1:51),5002)(eps(3,j),j=1,3)
               line=line+3
               eps=eps/1000.0d0 !To GHz
            end if

            if(FN(1:49)==' Dipole moment (field-independent basis, Debye):' .and. .not.dip_found)then !by default, at the end of output
               read(unitt,2000)FN
               read(FN(1:79),207)(dip_ptr(j),j=1,3)
207            format(3(13x,f15.4))
               line=line+1
               dip_ptr=dip_ptr*SI2AU_debye
            endif
            
            if(FN(2:23)=='Dipole moment (Debye):' .and. .not.dip_found)then !from Output(Pickett) option, this one is presumed to be the last dipole moment
               read(unitt,2000)FN
               read(FN(1:45),208)(dip_ptr(j),j=1,3)
208            format(3(1x,f14.7))
               dip_ptr=dip_ptr*SI2AU_debye
               line=line+1
               DipInPax=.true.
            endif
            
            IF(FN(2:40)=='Paramagnetic susceptibility tensor (au)')then
               do k=1,3
                  read(unitt,203)(xp(k,j),j=1,3)
                  line=line+1
               end do
   203         format(3(6x,f15.4))
            endif
                     
            IF(FN(2:43)=='Quadrupole moment (field-independent basis' .and. .not.q_found)then
               read(unitt,204)quad(1,1),quad(2,2),quad(3,3)
               read(unitt,204)quad(1,2),quad(1,3),quad(2,3)
   204         format(3(6x,f20.4))
               quad(2,1)=quad(1,2)
               quad(3,1)=quad(1,3)
               quad(3,2)=quad(2,3)
               do k=1,3
                  do j=1,3
                     quad(k,j)=quad(k,j)/bohr*SI2AU_debye
                  end do
               end do

               quadtr(1,1)=(2.0d0*quad(1,1)-quad(2,2)-quad(3,3))/3.0d0
               quadtr(2,2)=(2.0d0*quad(2,2)-quad(3,3)-quad(1,1))/3.0d0
               quadtr(3,3)=(2.0d0*quad(3,3)-quad(1,1)-quad(2,2))/3.0d0
               quadtr(1,2)=quad(2,1)
               quadtr(2,1)=quad(1,2)
               quadtr(1,3)=quad(3,1)
               quadtr(3,1)=quad(1,3)
               quadtr(2,3)=quad(3,2)
               quadtr(3,2)=quad(2,3)
               line=line+2
            endif

            IF(FN(2:28)=='Rotational constants (GHZ):' .and. .not.abc_found)then
               read(FN(29:48),2316)abc(1)
               read(FN(49:68),2316)abc(2)
               read(FN(69:88),2316)abc(3)
            endif    
            
            IF(FN(2:28)=='Rotational constants (MHZ):' .and. .not.abc_found)then
               read(2,2000)FN
               read(FN(2:15),2317)abc(1)
               read(FN(17:30),2317)abc(2)
               read(FN(32:45),2317)abc(3)
               line=line+1
               abc=abc*1d-3
            endif       
   2316     format(f19.7)  
   2317     format(f14.7)  

            if(FN(2:11)=='g tensor [')then
               read(unitt,2000)FN
               read(FN(7:21),*)g_el(1,1)
               read(FN(28:42),*)g_el(2,1)
               read(FN(49:63),*)g_el(3,1)
               read(unitt,2000)FN
               read(FN(7:21),*)g_el(1,2)
               read(FN(28:42),*)g_el(2,2)
               read(FN(49:63),*)g_el(3,2)
               read(unitt,2000)FN
               read(FN(7:21),*)g_el(1,3)
               read(FN(28:42),*)g_el(2,3)
               read(FN(49:63),*)g_el(3,3)
               line=line+3
            end if
            
            if(FN(2:29)=='Rotation matrix to Principal')then
               read(unitt,2000)FN
               read(unitt,2000)FN
               read(FN(10:50),*)(rotmat(1,k),k=1,3)
               read(unitt,2000)FN
               read(FN(10:50),*)(rotmat(2,k),k=1,3)
               read(unitt,2000)FN
               read(FN(10:50),*)(rotmat(3,k),k=1,3)
   45666       format(3(D13.6))
               line=line+4
               rot_found=.true.
            end if
            
            IF(FN(2:24)=='GIAO rotational tensor:' .and. .not.g_found)then
               read(unitt,2000)FN
               
               read(FN(7:17),*)g(1,1)
               read(FN(24:34),*)g(2,1)
               read(FN(41:51),*)g(3,1)
               read(unitt,2000)FN
               read(FN(7:17),*)g(1,2)
               read(FN(24:34),*)g(2,2)
               read(FN(41:51),*)g(3,2)
               read(unitt,2000)FN
               read(FN(7:17),*)g(1,3)
               read(FN(24:34),*)g(2,3)
               read(FN(41:51),*)g(3,3)
               line=line+3
205            format(3(6x,f11.4))
               g=transpose(g)
            endif
            
            if(FN(2:8).eq.'Charge=')then
               read(FN(14:28),*)charge
            end if
         end do
      end do
      
      mass=1d0
50    close(unitt)
      
      nat=natoms
!     molecular mass
      mass=0d0
      do i = 1,nat
         idx = iz(i)
         mass=mass+a_masses(i)
      end do

!     coordinates to atomic units:
      do ia=1,nat
         do i=1,3
            r(i,ia)=r(i,ia)/bohr
         end do
      end do


      CM=CALC_COM(nat,a_masses,r)
      
      
! !     subtract mass center:
      ! DO I=1,3
         ! CM(I)=0.0d0
         ! mt=0.0d0
         ! DO J=1,nat
            ! mt=mt+ama(J)
            ! CM(I)=CM(I)+R(I,J)*ama(J)
         ! end do   
         ! CM(I)=CM(I)/mt
      ! end do

      unitt=2
      if(flags(1))then
         open(unitt,file=outputfilen//'.molp')
      end if
      !dip=0d0
      
      if(.not.dipinpax)then
         if(.not.dip_found)then 
            do j = 1,3
               dip(j)=dip(j)+charge*CM(j)
               dip_nuc(j)=dip_nuc(j)+charge*CM(j)
            end do
         end if
      end if
      t0=CALC_INERT(nat,a_masses,r)
      
      do ia=1,nat
         do i=1,3
            r(i,ia)=r(i,ia)-CM(i)
         end do
      end do
      t1=CALC_INERT(nat,a_masses,r)
      
      !dip=dipbuf
      quadbuf=0d0
      
      !quadrupole translation
      if(.not.q_found)then 
         do j = 1,3
            do k = 1,3
               quadbuf(j,k)=quad(j,k)-CM(j)*dip(k)-CM(k)*dip(j)
            end do
         end do
         quad=quadbuf
      end if
      
      !g tensor translation
      if(.not.g_found)then 
         do i=1,3
            do j=1,3
               g(i,j)=g(i,j)*t0(j,j)-mpr_amu*(2d0*kdelta(i,j)*DOT_PRODUCT(dip,CM)-dip(i)*CM(j)-CM(i)*dip(j))
               g(i,j)=g(i,j)/t1(j,j)
            end do
         end do
      end if  
      
      
!     make inertia-like:
      qt=0.5d0*quad
      qtot=0d0
      do i=1,3
         do j=1,3
            if(i.eq.j)then
               k=i+1
               if(k.gt.3)k=1
               l=k+1
               if(l.gt.3)l=1
               qtot(i,j)=qt(k,k)+qt(l,l)
            else
               qtot(i,j)=-qt(i,j)
            endif
         end do
      end do
      quadbuf=0d0
      outp=.true.
      if(flags(1))then
         write(unitt,*)'mass center (A):'
         write(unitt,6001)(cm(i)*bohr,i=1,3)
      end if
   !  inertia and quadrupole moment:
      do 6 i=1,3
         t(i,i)=0.0d0
         q(i,i)=0.0d0
         j=i+1
         if(j.gt.3)j=1
         k=j+1
         if(k.gt.3)k=1
         do 61 ia=1,nat
            q(i,i)=q(i,i)+dble(iz(ia))*(r(j,ia)*r(j,ia)+r(k,ia)*r(k,ia))
            t(i,i)=t(i,i)+a_masses(ia)*(r(j,ia)*r(j,ia)+r(k,ia)*r(k,ia))
61       continue
         do 9 j=1,i-1
            t(i,j)=0.0d0
            q(i,j)=0.0d0
            do 62 ia=1,nat
               q(i,j)=q(i,j)-dble(iz(ia))*r(i,ia)*r(j,ia)
               t(i,j)=t(i,j)-a_masses(ia)*r(i,ia)*r(j,ia)
62          continue
            q(j,i)=q(i,j)
            t(j,i)=t(i,j)
9        continue
6     continue
         
      t=0d0
      q=0d0
      do i = 1,nat
         do j=1,3
            do k=1,3
               q(j,k)=q(j,k)+dble(iz(i))*(dble(kdelta(j,k))*DOT_PRODUCT(r(:,i),r(:,i))-r(j,i)*r(k,i))
            end do
         end do
         ! q(1,1)=q(1,1) + dble(iz(i))*(r(2,i)**2+r(3,i)**2)
         ! q(2,2)=q(2,2) + dble(iz(i))*(r(1,i)**2+r(3,i)**2)
         ! q(3,3)=q(3,3) + dble(iz(i))*(r(2,i)**2+r(1,i)**2)
         ! q(1,2)=q(1,2) - dble(iz(i))*r(1,i)*r(2,i)
         ! q(1,3)=q(1,3) - dble(iz(i))*r(1,i)*r(3,i)
         ! q(2,3)=q(2,3) - dble(iz(i))*r(2,i)*r(3,i)
      
         t(1,1)=t(1,1) + a_masses(i)*(r(2,i)**2+r(3,i)**2)
         t(2,2)=t(2,2) + a_masses(i)*(r(1,i)**2+r(3,i)**2)
         t(3,3)=t(3,3) + a_masses(i)*(r(2,i)**2+r(1,i)**2)
         t(1,2)=t(1,2) - a_masses(i)*r(1,i)*r(2,i)
         t(1,3)=t(1,3) - a_masses(i)*r(1,i)*r(3,i)
         t(2,3)=t(2,3) - a_masses(i)*r(2,i)*r(3,i)
      end do
      t(3,2)=t(2,3)
      t(3,1)=t(1,3)
      t(2,1)=t(1,2)
      q(3,2)=q(2,3)
      q(3,1)=q(1,3)
      q(2,1)=q(1,2)
      
      ! t=AxesSwitch_T(t,[1,3,2])
      !if molecular axes are not ordered from the lowest to the highest moment of inertia
      if(.not.(t(1,1)<t(2,2).and.t(2,2)<t(3,3)))then
         call Reorder_Tensor(t,order)
         do i = 1,nat
            r(:,i)=AxesSwitch_V(r(:,i),order)
         end do
         if(.not.dip_found.and..not.dipinpax)then
            dip = AxesSwitch_V(dip,order)
            dip_nuc = AxesSwitch_V(dip_nuc,order)
         end if
         if(.not.g_found)then
            g = AxesSwitch_T(g,order)
         end if
         if(.not.q_found)then
            quad = AxesSwitch_T(quad,order)
         end if
      end if
         !write(6,*)'t matrix:'
         !call WriteMat(t,3,3)
         !write(6,*)'q matrix:'
         !call WriteMat(q,3,3)
         !write(6,*)'r matrix:'
         !call WriteMat(r,3,n0)
         
         e=t
         !CALL WriteMat(T,3,3)
         ierr=0_4
         !CALL TRED12(3,T,A,Z,1,IERR,E)
         call Diag(a,t,3,ierr,.false.,real(0.1,kind(a(1))))
         !a_help=-a_help
         call Select_sort(a,3,order)
         !a_help=-a_help
         call Reorder_arr_2d(t,3,order,.true.)
         !CALL WriteMat(T,3,3)
         !CALL WriteMat(T,arr(1),arr(2))
         !arr=shape(A)
         !write(6,*)(a(i),i=1,3)
         !write(6,*)ierr
!     check that we have a right-handed system
         ! if(t(3,3)*(t(1,1)*t(2,2)-t(1,2)*t(2,1)) < 0d0)then
            ! t(1,3)=-t(1,3)
            ! t(2,3)=-t(2,3)
            ! t(3,3)=-t(3,3)
            ! if(flags(1)) write(unitt,*)'chirality fixed'           
         ! endif
         
         z=0d0
         do i=1,3
            z(i,i)=a(i)
         end do
         a=u*c/1.0d9/a/2.0d0
         !z=e-z
         !z=matmul(z,t(:,1))
         !call writemat(MATMUL(e,t),3,3)
         !call writemat(MATMUL(t,z),3,3)
         !Write(6,*)'Eigenvalues:'
         !Write(6,*)(Z(i),i=1,3)
         !Write(6,*)'Eigenvector:'
         !Write(6,*)(E(i),i=1,3)
         
         write(6,*)'Centre of mass shift (au):'
         write(6,6000)cm
         write(6,*)'Custom rotation angles to principal system:'
         write(6,6000)rotmat_2_angles(t)
         write(6,*)'Gaussian rotation angles to principal system:'
         write(6,6000)rotmat_2_angles(rotmat)
         
         if(rot_found)then
            t=rotmat
         end if
         
         arrbuf=0d0
         !t=transpose(t)
         do i = 1,3
            do j = 1,3
               do k = 1,3
                  do ii = 1,3
                     arrbuf(i,j)=arrbuf(i,j)+t(i,k)*t(j,ii)*e(k,ii)
                  end do
               end do
            end do
         end do
         
         !arrbuf=RotateT(e,t)

         if(flags(1))then
            write(unitt,*)'Gaussian (GHz):'
            write(unitt,6000)(abc(i),i=1,3)
            write(unitt,*)'Recalc:'
            write(unitt,6000)(a(i)/amu,i=1,3)
         end if
6000     format(3f18.6)

         if(.not.g_found)then
            do i=1,3
               do j=1,3
                  ge(i,j)=-4.0d0*xp(i,j)/z(i,i)*mp
                  gn(i,j)=q(i,j)/z(i,i)*mp
                  gt(i,j)=ge(i,j)+gn(i,j)
               end do
            end do

            do i=1,3
               do j=1,3
                  gq(i,j)=qtot(i,j)/z(i,i)*mp
               end do
            end do
         else if(g_found .and. separateElNuc)then
            do i=1,3
               do j=1,3
                  gn(i,j)=q(i,j)/e(i,i)*mp
               end do
            end do
            ge=gt-gn
         end if

         !dip_el=dip-dip_nuc
         if(flags(1))then
         write(unitt,*)'Laboratory system:'
         ! write(6,*)'Positions:'
         ! do i = 1,nat
            ! write(6,6382)i,(r(j,i),j=1,3)
! 6382     format(I4,': ',3F18.6)       
         ! end do
         write(unitt,*)
         write(unitt,*)'Trafo:'
         write(unitt,6001)((t(i,j),i=1,3),j=1,3)
         write(unitt,*)'xp:'
         write(unitt,6001)((xp(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor,electronic:'
         write(unitt,6001)((ge(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor,nuclear:'
         write(unitt,6001)((gn(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor,electronic+nuclear:'
         write(unitt,6001)((gt(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor in output:'
         write(unitt,6001)((g(i,j),i=1,3),j=1,3)
         write(unitt,*)'electronic g-tensor in output:'
         write(unitt,6001)((g_el(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment (1/2) qrr:'
         write(unitt,6001)((qt(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment nuclear:'
         write(unitt,6001)((q(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment total, inertia-like:'
         write(unitt,6001)((qtot(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment total, traceless:'
         write(unitt,6001)((quadtr(i,j),i=1,3),j=1,3)
         write(unitt,*)'Inertia moment:'
         write(unitt,6001)((e(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor from quadrupole:'
         write(unitt,6001)((gq(i,j),i=1,3),j=1,3)
         write(unitt,*)'dipole (au):'
         write(unitt,6001)(dip(i),i=1,3)
         if(separateElNuc)then
            write(unitt,*)'electronic dipole (au):'
            write(unitt,6001)(dip_el(i),i=1,3)
            write(unitt,*)'nuclear dipole (au):'
            write(unitt,6001)(dip_nuc(i),i=1,3)
         end if
         end if
6001     format(/,3(3f18.6,/),/)

         !write(6,*)'dipole (au):'
         !write(6,6001)(dip(i),i=1,3)
         do i = 1,nat
            r(:,i)=rotatev(r(:,i),t)
         end do
         
         ! arrbuf=0d0
         ! do i = 1,3
            ! do j = 1,3
               ! bufr=0d0
               ! do k = 1,3
                  ! do m = 1,3
                     ! bufr=bufr+t(i,k)*t(j,m)*xp(k,m)
                  ! end do 
               ! end do
               ! arrbuf(i,j)=bufr
            ! end do
         ! end do
         
         if(.not.dip_found.and..not.dipinpax)then
            dip=RotateV(dip,t)
            dip_nuc=RotateV(dip_nuc,t)
         end if
         xp=RotateT(xp,t)
         e=RotateT(e,t)
         q=CALC_Q(nat,iz,r)
         qtot=Rotatet(qtot,t)
         quadtr=Rotatet(quadtr,t)
         g=rotatet(g,t)
         gq=Rotatet(gq,t)
         qt=rotatet(qt,t)
         g_el=rotatet(g_el,t)
         
         !quad=rotatet(quad,t)
         !arrbuf=rotatet(t,t)
         if(.not.g_found)then
            do i=1,3
               do j=1,3
                  ge(i,j)=-4.0d0*xp(i,j)/z(i,i)*mp
                  gn(i,j)=q(i,j)/z(i,i)*mp
                  gt(i,j)=ge(i,j)+gn(i,j)
                  gq(i,j)=qtot(i,j)/z(i,i)*mp
               end do
            end do
         else if(g_found .and. separateElNuc)then
            do i=1,3
               do j=1,3
                  gn(i,j)=q(i,j)/z(i,i)*mp
               end do
            end do
            ge=gt-gn
         end if
         
         !dip_el=dip-dip_nuc
         if(flags(1))then
         write(unitt,*)'Molecular system:'
         ! write(6,*)'Positions:'
         ! do i = 1,nat
            ! write(6,6382)i,(r(j,i),j=1,3)
         ! end do
         write(unitt,*)
         write(unitt,*)'xp:'
         write(unitt,6001)((xp(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor,electronic:'
         write(unitt,6001)((ge(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor,nuclear:'
         write(unitt,6001)((gn(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor,electronic+nuclear:'
         write(unitt,6001)((gt(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor in output:'
         write(unitt,6001)((g(i,j),i=1,3),j=1,3)
         write(unitt,*)'electronic g-tensor in output:'
         write(unitt,6001)((g_el(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment (1/2) qrr:'
         write(unitt,6001)((qt(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment nuclear:'
         write(unitt,6001)((q(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment total, inertia-like:'
         write(unitt,6001)((qtot(i,j),i=1,3),j=1,3)
         write(unitt,*)'Quadrupole moment total, traceless:'
         write(unitt,6001)((quadtr(i,j),i=1,3),j=1,3)
         write(unitt,*)'Inertia moment:'
         write(unitt,6001)((e(i,j),i=1,3),j=1,3)
         write(unitt,*)'g-tensor from quadrupole:'
         write(unitt,6001)((gq(i,j),i=1,3),j=1,3)
         write(unitt,*)'dipole (au):'
         write(unitt,6001)(dip(i),i=1,3)
         if(separateElNuc)then
            write(unitt,*)'electronic dipole (au):'
            write(unitt,6001)(dip_el(i),i=1,3)
            write(unitt,*)'nuclear dipole (au):'
            write(unitt,6001)(dip_nuc(i),i=1,3)
         end if
         close(unitt)
         write(6,*)'WRITTEN '//outputfilen//'.molp'
         end if
      
      t=calc_inert(nat,a_masses,r)
      dip_nuc=0d0
      do i = 1,nat
         dip_nuc=dip_nuc+r(:,i)*iz(i)
      end do
      
      
      if(.not.q_found)then
0451     quad=qt
      end if
      if(.not.g_found)then
         gt=g
      else
         gt=repAxesSwitch_T(gt,rep,lh,.true.)
      end if
      if(dip_found)then
         dip=repAxesSwitch_V(dip,rep,lh,.true.)
      end if
      dip_el=dip-dip_nuc
      ge=gt-gn
      if(dip_found)then
         dip=repAxesSwitch_V(dip,rep,lh)
      end if
      if(g_found)then
         gt=repAxesSwitch_T(gt,rep,lh)
      end if
      
      return
      
1111  write(6,*)'GAUSSIAN OUTPUT NOT FOUND: '//filename
      stop
      end subroutine readFile
      
      function RotateV(v,T) result(res)
         double precision v(3),T(3,3),res(3)
         INTEGER i,j
         
         !res=0.0d0
         res=matmul(t,v)
         ! do i=1,3
            ! do j=1,3
               ! res(i)=res(i)+T(j,i)*v(j)
            ! end do
         ! end do
      end function RotateV
      
      function RotateT(v,t) result(res)
         double precision v(3,3),T(3,3),res(3,3)
         INTEGER i,j,k,l
         
         ! res=0.0d0
         res=matmul(t,v)
         ! do i=1,3
            ! do j=1,3
               ! do k=1,3
                  ! do l=1,3
                     ! res(i,j)=res(i,j)+t(k,i)*t(l,j)*v(k,l)
                  ! end do
               ! end do
            ! end do
         ! end do
      end function RotateT
      
!     write 2d matrix to unit 6 (the console)
      subroutine WriteMat(mat,x,y)
         INTEGER :: x,y,i,j
         double precision :: mat(x,y)
         write(6,*)
         do i = 1,x
            do j = 1,y
               write(6,2000,advance='no')mat(i,j)
            end do
            write(6,*)
         end do
2000     format((f11.5))
         write(6,*)
      end subroutine WriteMat
      
      subroutine WriteMat_c(mat,x,y)
         INTEGER :: x,y,i,j
         double complex :: mat(x,y)
         write(6,*)
         do i = 1,x
            do j = 1,y
               write(6,2000,advance='no')dble(mat(i,j))
            end do
            write(6,*)
         end do
         
         write(6,*)
         do i = 1,x
            do j = 1,y
               write(6,2001,advance='no')aimag(mat(i,j))
            end do
            write(6,*)
         end do
         
2000     format((f11.5))
2001     format((f11.5),'i')
         write(6,*)
      end subroutine WriteMat_c
      
      subroutine WriteMat_L(mat,x,y)
         INTEGER :: x,y,i,j
         logical :: mat(x,y)
         write(6,*)
         do i = 1,x
            do j = 1,y
               write(6,9468,advance='no')mat(i,j)
            end do
            write(6,*)
         end do
9468     format(1X,1L)
         write(6,*)
      end subroutine WriteMat_L
      
      
!     generic error handling, writes argument
      subroutine report(s)
      character(*) s
      write(6,*)s
      stop 1
      end subroutine report

      
      
end module input

module spectra
   use constants
   use arrayOperations
   implicit none
   !molecular and measurement parameters
   double precision :: temp_,mass_,fmin_,fmax_,resol_,calwidth_,pressMes_,pressCal_,tempCal_,molDens_
   
   private temp_,mass_,fmin_,fmax_,resol_,calwidth_,pressMes_,pressCal_,tempCal_,molDens_
   
 
CONTAINS
 
   
   !Shamelessly borrowed cause i'm too dumb too write qsort
   !https://gist.github.com/t-nissie/479f0f16966925fa29ea
   recursive subroutine quicksort(a, first, last,order)
      double precision  a(*), x, t
      INTEGER first, last, order(*),tt
      INTEGER i, j

      x = a( (first+last) / 2 )
      i = first
      j = last
      do
         do while (a(i) < x)
           i=i+1
         end do
         do while (x < a(j))
           j=j-1
         end do
         if (i >= j) exit
         t = a(i);  a(i) = a(j);  a(j) = t
         tt=order(i); order(i)=order(j); order(j)=tt
         i=i+1
         j=j-1
      end do
      if (first < i-1) call quicksort(a, first, i-1,order)
      if (j+1 < last)  call quicksort(a, j+1, last,order)
   end subroutine quicksort   
   
   pure function Pressure2MolDens(press,temp)result(res)
      double precision,intent(in) :: press,temp
      double precision res
      res=press*Na/(Rgas*temp)
   end function Pressure2MolDens
   
   pure function CONST2MOM(A) result(res)
      double precision,intent(in) :: A
      double precision res
      res=h/(8*pi**2*cc*A)
   end function CONST2MOM
   
   !Expected to start from J=0; the higher end J, the better "accuracy"
   !ens in Hz
   pure function ExactPARTFUN(ens,n,temp,sym,endJ) result(res)
      Integer,intent(in) :: n,sym,endJ
      double precision,intent(in) :: ens(n),temp
      double precision res,buf
      INTEGER i,j,k,car
      
      res=0d0
      car=0
      do j=0,endJ
         buf=0d0
         do k = -j,j
            car=car+1
            buf=buf+dble(2*J+1)*exp(-h*ens(car)/(kb*temp))
         end do
         res=res+buf!*dble(2*j+1)
      end do
      res=res/sym
   end function ExactPARTFUN
   
   pure function ExactPARTFUN_SR(ens,n,temp,sym,endJ) result(res)
      Integer,intent(in) :: n,sym,endJ
      double precision,intent(in) :: ens(n),temp
      double precision res,buf
      INTEGER i,j,k,car
      
      res=0d0
      car=0
      do j=0,endJ
         buf=0d0
         do k = -j,j
            car=car+1
            buf=buf+dble(2*J+1)*exp(-h*ens(car)/(kb*temp))
         end do
         res=res+buf!*dble(2*j+1)
      end do
      res=res/sym
   end function ExactPARTFUN_SR
   
   
   
   !Approximate partition function from rotational constants
   !A,B,C in Hz
   pure function PARTFUN(A,B,C,temp,symnum) result(res)
      Integer,intent(in) :: symnum 
      double precision,intent(in) :: a,b,c,temp
      double precision res,radconst
      !res=(8d0*pi**2*kb*temp/h**2)**(1.5d0)*(pi*CONST2MOM(A)*CONST2MOM(B)*CONST2MOM(C))**(0.5d0)
      !res=5.34d6/dble(symnum)*sqrt(temp**3/(A*B*C)*(1d3**3))
      parameter(radconst=h*cc/kb)
      res=sqrt(pi*temp**3/((A*B*C*1d9**3*Hz_2_cm**3*1d6)*radconst**3))
   end function PARTFUN
    
   !state population in percent
   !nu in Hz
   elemental function Calc_Prob(nu,J,Q,temp) result(res)
      Integer,intent(in) :: J
      double precision,intent(in) :: nu,Q,temp
      double precision res
      res=dble(2*J+1)*exp(-(h*nu)/(kb*temp))/Q
   end function Calc_Prob
    
   !Calculates percentual populations of given J,Tau energy states, Ms are degenerate
   !The sum of these is expected to give unity if accurate partition function is used
   !freqs in Hz
   function Calc_Frac_Populations(freqs,Js,N,Q,temp) result(res)
      Integer,intent(in) :: N,Js(N)
      double precision,intent(in) :: freqs(N),Q,temp
      double precision :: res(N),rbuf
      INTEGER i
      
      ! do i=1,N
         ! res(i)=Calc_Prob(freqs(i),Js(i),Q,temp)
      ! end do
      res=Calc_Prob(freqs,Js,Q,temp) !is this undefined behaviour I wonder, ie. passing 2 concurrent arrays to an elemental function
      
      !rbuf=sum(res)
      !rbuf=MAXVAL(res)
      !res=res!/rbuf
   end function Calc_Frac_Populations
   
   !get index of given J,T state in the state population array (given state pop. = (exp(-h*nu/(kb*t)))/(partition function))
   !J T idx
   !0 0   1
   !1-1   2
   !1 0   3
   !1 1   4
   !2-2   5
   !...
   pure function Calc_Frac_index(J,T)result(idx)
      Integer,intent(in) :: J,T
      INTEGER idx
      idx=j**2+j+t+1
   end function Calc_Frac_index
   
   
   
   !find the smallest difference in energy transitions
   function Find_Smallest_Gap(trans_ens,num) result(gap)
      INTEGER num,i,j
      double precision trans_ens(num),gap,buf
   
      gap=HUGE(gap)
      do i = 1,num
         if(trans_ens(i)<0d0)STOP 'Transition frequency cannot be negative'
      end do
      
      do i = 1,num-1
         do j = i+1,num
            buf=abs(trans_ens(i)-trans_ens(j))
            if(buf<gap)then
               gap=buf
            end if
         end do
      end do
   end function Find_Smallest_Gap
   
   !are doubles equal?
   function D_EQ(num1,num2,tol)result(res)
      double precision num1,num2,tol
      logical res
      
      if(dabs(num1-num2)<=tol)then
         res = .true.
      else
         res = .false.
      end if
   end function D_EQ
   
   !gets the indices where overlap of intervals x1,x2 starts
   subroutine Intervals_Overlap_Start(x1,n,endIdx1,x2,m,endIdx2,idx1,idx2,tol)
      Integer,intent(in) :: n,m,endIdx1,endIdx2
      double precision,intent(in) :: x1(n),x2(m),tol
      Integer,intent(out) :: idx1,idx2
      INTEGER i,j,k
      do i = 1,endIdx1
         do j = 1,endIdx2
            if(D_EQ(x1(i),x2(j),tol))then
               idx1=i
               idx2=j
               return
            end if
         end do
      end do
   end subroutine Intervals_Overlap_Start
   
   !gets the indices where overlap of intervals x1,x2 ends
   subroutine Intervals_Overlap_End(x1,n,endIdx1,x2,m,endIdx2,idx1,idx2,tol)
      Integer,intent(in) :: n,m,endIdx1,endIdx2
      double precision,intent(in) :: x1(n),x2(m),tol
      Integer,intent(out) :: idx1,idx2
      INTEGER i,j,k
      do i = endIdx1,1,-1
         do j = endIdx2,1,-1
            if(D_EQ(x1(i),x2(j),tol))then
               idx1=i
               idx2=j
               return
            end if
         end do
      end do
   end subroutine Intervals_Overlap_End
   
   !gets where overlap of two bands starts and ends for both bands
   function Intervals_Overlap_Size(x1,n,endIdx1,x2,m,endIdx2,tol) result(res)
      Integer,intent(in) :: n,m,endIdx1,endIdx2
      double precision,intent(in) :: x1(n),x2(m),tol
      INTEGER res
      INTEGER ovs1,ovs2,ove1,ove2
      
      call Intervals_Overlap_Start(x1,n,endIdx1,x2,m,endIdx2,ovs1,ovs2,tol)
      call Intervals_Overlap_End(x1,n,endIdx1,x2,m,endIdx1,ove1,ove2,tol)
      
      if((ovs2==1).AND.(ove2==endIdx2).AND.(endIdx1>endIdx2))then
         res=endIdx1
         return
      else if((ovs1==1).and.(ove1==endIdx1).AND.(endIdx1<endIdx2))then
         res=endIdx2
         return
      end if
      
      
      if((x1(1)>x2(1)).AND.(x1(endIdx1)>x2(endIdx2)))then
         res=ovs2+(ove2-ovs2)+(endIdx1-ove1)
      else if((x1(1)<x2(1)).AND.(x1(endIdx1)<x2(endIdx2)))then
         res=ovs1+(ove1-ovs1)+(endIdx2-ove2)
      else if((x1(1)>x2(1)).AND.(x1(endIdx1)<x2(endIdx2)))then
         res=endIdx2
      else if((x1(1)<x2(1)).AND.(x1(endIdx1)>x2(endIdx2)))then
         res=endidx1
      else if((D_EQ(x1(1),x2(1),tol)).AND.(D_EQ(x1(endIdx1),x2(endIdx2),tol)))then
         res=endidx1
      else
         STOP 'iz wronk overlap size'
      end if
   end function Intervals_Overlap_Size
   
   !assumes ordered values
   !checks if two bands overlap
   function Intervals_Overlap(x1,n,endIdx1,x2,m,endIdx2)result(res)
      INTEGER n,m,endIdx1,endIdx2
      double precision x1(n),x2(m)
      logical res
      
      res=.false.
      if((x1(1)<=x2(endIdx2)).and.(x1(1)>=x2(1)))then 
         res=.true.
      else if((x2(endIdx2)>=x1(1)).and.(x2(endIdx2)<=x1(endIdx1)))then
         res=.true.
      else if((x1(endIdx1)>=x2(1)).and.(x1(endIdx1)<=x2(endIdx2)))then
         res=.true.
      else if((x2(1)<=x1(endIdx1)).and.(x2(1)>=x1(1)))then
         res=.true.
      end if
   end function Intervals_Overlap
   
   !nu - absorbing frequency in Hz (or GHz, MHz, KHz...)
   !temp - temperature in Kelvin
   !mass - mass of a single molecule in kg
   !result - Doppler width in units of nu
   function Doppler_Width(nu,temp,mass) result(res)
      double precision nu,temp,mass,res
      res=2d0*sqrt(2*kb*temp/(mass))*nu/cc/2.3548d0
   end function Doppler_Width
   
   
   !Pressure broadening, argument breadth at 1mm Hg and 300K, press in mm Hg, temp in K 
   function Pressure_width(breadth,t1,p1,t2,p2) result(res)
      double precision,intent(in) :: breadth,t1,p1,t2,p2
      double precision res
      res = t1*p2*breadth/(p1*t2)
   end function Pressure_width
   
   
   !returns a 2d array of logicals whether given bands overlap, diagonals are a band with itself, (2,1)==(1,2)
   function Find_band_overlaps(x,n,endIdx,siz,tol)result(res)
      INTEGER n,siz,endIdx(n)
      double precision x(n,siz),tol
      logical res(n,n)
      INTEGER i,j
      
      res=.false.
      do i = 1,n
         do j = 1,n
            if(i==j)then
               res(i,j)=.true.
               cycle
            end if
            res(i,j) = Intervals_Overlap(x(i,1:endIdx(i)),size(x(i,1:endIdx(i)),dim=1),endIdx(i),x(j,1:endidx(j)),size(x(j,1:endIdx(j)),dim=1),endIdx(j))
         end do
      end do
   end function find_band_overlaps
   
   !overlaps two bands x1,y1 and x2,y2 into x3,y3
   !endIdxs are the ending indices of arrays which have to be larger than the amount of information in them (workaround cause no Lists)
   !resol is the distance at which two points should overlap,
   !
   !Author's notes:
   !uses simple Lambert-Beer law, ie. only adds two y points
   !it works but there is a great loss of accuracy the more it is used, you would need to supply it the exact same points in every x array
   !the initial idea was to save HDD space while preserving super narrow bands in a SINGLE file for a SUPER WIDE spectrum
   subroutine Band_overlap(x1,y1,n,endIdx1,x2,y2,m,endIdx2,x3,y3,r,startIdx3,resol,info)
      Integer, intent(in) :: endIdx1,endIdx2,startIdx3,n,m,r
      double precision, intent(in) :: x1(n),y1(n),x2(m),y2(m),resol
      double precision, intent(out) :: x3(r),y3(r)
      logical, intent(out) :: info
      double precision :: curRes
      INTEGER i,j,idx,ovs1,ovs2,ove1,ove2
           
           
      if(x1(endIdx1)==0)then
         stop 'x1 zero'
      end if
      if(x2(endidx2)==0)then
         stop 'x2 zero'
      end if
              
              
      !check for not overlapping bands
      if(Intervals_Overlap(x1,n,endIdx1,x2,m,endIdx2).eqv..false.)then
         j=startIdx3-1
         if(x1(1)<x2(1))then
            do i = 1,endIdx1
               j=j+1
               x3(j)=x1(i)
               y3(j)=y1(i)
            end do
            
            do i = 1,endIdx2
               j=j+1
               x3(j)=x2(i)
               y3(j)=y2(i)
            end do
         else if(x2(1)<x1(1))then
            do i = 1,endIdx2
               j=j+1
               x3(j)=x2(i)
               y3(j)=y2(i)
            end do
            
            do i = 1,endIdx1
               j=j+1
               x3(j)=x1(i)
               y3(j)=y1(i)
            end do
         end if
         do i = 1,j-1
            if(x3(i)>x3(i+1))then
            stop 'wrong overlap 1'
            end if
         end do
         info=.false.
         return
      end if
      
      call Intervals_Overlap_Start(x1,n,endIdx1,x2,m,endIdx2,ovs1,ovs2,resol)
      call Intervals_Overlap_End(x1,n,endIdx1,x2,m,endIdx2,ove1,ove2,resol)
           
      if((ovs1==1).AND.(ovs2==1).AND.(ove2==endidx2).AND.(ove1==endidx1))then
         if(endidx1>endidx2)then
            do i = 1,endidx2
               x3(i)=(x2(i)+x1(i))/2d0
               y3(i)=y2(i)+y1(i)
            end do
            curRes=x3(i-1)-x3(i-2)
            do j = i,endidx1
               x3(j)=x3(j-1)+curRes
               y3(j)=y1(j)
            end do
         else if(endIdx1==endIdx2) then
            do i = 1,endidx2
               x3(i)=(x2(i)+x1(i))/2d0
               y3(i)=y2(i)+y1(i)
            end do
            j=endidx1+1
         else if(endIdx2>endIdx1) then
            do i = 1,endidx1
               x3(i)=(x2(i)+x1(i))/2d0
               y3(i)=y2(i)+y1(i)
            end do
            curRes=x3(i-1)-x3(i-2)
            do j = i,endidx2
               x3(j)=x3(j-1)+curRes
               y3(j)=y2(j)
            end do
         end if
         do i = 1,j-2
            if(x3(i)>x3(i+1))then
            stop 'wrong overlap 5'
            end if
         end do
         return
      end if
           
           
      !x2 is inside x1
      if((ovs2==1).AND.(ove2==endIdx2).AND.(x1(1)<x2(1)).AND.(x1(endIdx1)>x2(endIdx2)))then
         idx=startIdx3
         do i = 1,ovs1-1
            x3(idx)=x1(i)
            y3(idx)=y1(i)
            idx=idx+1
         end do
         j = ovs1
         do i = ovs2,ove2
            x3(idx)=(x2(i)+x1(j))/2d0
            y3(idx)=y2(i)+y1(j)
            idx=idx+1
            j=j+1
         end do
         do i = ove1+1,endIdx1
            x3(idx)=x1(i)
            y3(idx)=y1(i)
            idx=idx+1
         end do
         do i = 1,idx-2
            if(x3(i)>x3(i+1))then
            stop 'wrong overlap 2'
            end if
         end do
         return
      !x1 is inside x2
      else if((ovs1==1).and.(ove1==endIdx1).AND.(x1(1)>x2(1)).AND.(x2(endIdx2)>x1(endIdx1)))then
         idx=startIdx3
         do i = 1,ovs2-1
            x3(idx)=x2(i)
            y3(idx)=y2(i)
            idx=idx+1
         end do
         j = ovs2
         do i = ovs1,ove1
            x3(idx)=(x1(i)+x2(j))/2d0
            y3(idx)=y1(i)+y2(j)
            idx=idx+1
            j=j+1
         end do
         do i = ove2+1,endIdx2
            x3(idx)=x2(i)
            y3(idx)=y2(i)
            idx=idx+1
         end do
         do i = 1,idx-2
            if(x3(i)>x3(i+1))then
            stop 'wrong overlap 3'
            end if
         end do
         return
      end if
      
      
      if((x1(1)<x2(1)).AND.(x1(endidx1)<x2(endidx2)))then
         idx=startIdx3
         do i = 1,ovs1-1
            x3(idx)=x1(i)
            y3(idx)=y1(i)
            idx=idx+1
         end do
         j = ovs2
         do i = ovs1,ove1
            if(j>endidx2)then
               x3(idx)=x1(i)
               y3(idx)=y1(i)
               exit
            end if
            x3(idx)=(x1(i)+x2(j))/2d0
            y3(idx)=y1(i)+y2(j)
            idx=idx+1
            j=j+1
         end do
         do i = ove2+1,endIdx2
            x3(idx)=x2(i)
            y3(idx)=y2(i)
            idx=idx+1
         end do
         
      else if((x1(1)>x2(1)).AND.(x1(endidx1)>x2(endidx2)))then
         idx=startIdx3
         do i = 1,ovs2-1
            x3(idx)=x2(i)
            y3(idx)=y2(i)
            idx=idx+1
         end do
         j = ovs1
         do i = ovs2,ove2
            if(j>endidx1)then
               x3(idx)=x2(i)
               y3(idx)=y2(i)
               exit
            end if

            x3(idx)=(x2(i)+x1(j))/2d0
            y3(idx)=y2(i)+y1(j)
            idx=idx+1
            j=j+1
         end do
         do i = ove1+1,endIdx1
            x3(idx)=x1(i)
            y3(idx)=y1(i)
            idx=idx+1
         end do
      end if
      !check for errors
      do i = 1,idx-2
         if(x3(i)>x3(i+1))then
         stop 'wrong overlap 4'
         end if
      end do
      info=.true.
   end subroutine Band_overlap
   
   recursive subroutine Make_Band_Overlaps(x_col,y_col,endIdxs,n,colSize,tole,xout,yout,m)
      INTEGER n,colSize,counter,endIdxs(colSize),curIndex,m,i,j
      double precision,allocatable :: x_col(:,:),y_col(:,:)
      double precision,allocatable :: x_col2(:,:),y_col2(:,:)
      double precision :: xout(m),yout(m),tole
      Integer,allocatable :: endIdxs2(:)
      logical infor
      
      curIndex=1
      counter=colsize/2+mod(colSize,2)
      allocate(x_col2(counter,n),y_col2(counter,n),endIdxs2(counter))
      do i = 1,colSize,2
         if((i+1)>colSize)then
            exit
         end if
         if(Intervals_overlap(x_col(i,:),n,endIdxs(i),x_col(i+1,:),n,endIdxs(i+1)).eqv..false.)then
         endIdxs2(curIndex) = endIdxs(i)+endIdxs(i+1)
         else
         endIdxs2(curIndex) = Intervals_Overlap_Size(x_col(i,:),n,endIdxs(i),x_col(i+1,:),n,endIdxs(i+1),tole)
         end if
         call Band_overlap(x_col(i,:),y_col(i,:),n,endIdxs(i),x_col(i+1,:),y_col(i+1,:),n,endIdxs(i+1),x_col2(curIndex,:),y_col2(curIndex,:),n,1,tole,infor)
         do while(x_col2(curIndex,endIdxs2(curIndex))==0d0)
            if(x_col2(curIndex,endIdxs2(curIndex))==0)then
               endIdxs2(curIndex)=endIdxs2(curIndex)-1
            end if
         end do
         do j = 1,endIdxs2(curIndex)-1
            if(x_col2(curIndex,j)>x_col2(curIndex,j+1))then
               stop 'wake me up inside' 
            end if
         end do
         !bufI=Intervals_Overlap_Size(x_col(i,:),,x_col(i+1,:))
         curIndex=curIndex+1
      end do
      if(mod(colSize,2)/=0)then
         x_col2(curIndex,:)=x_col(colSize,:)
         y_col2(curIndex,:)=y_col(colSize,:)
         endIdxs2(curIndex)=endIdxs(colSize)
      end if
      deallocate(x_col,y_col)
      if(counter==1)then
         xout=x_col2(1,1:endIdxs(1))
         yout=y_col2(1,1:endIdxs(1))
         return
      end if
      call Make_Band_Overlaps(x_col2,y_col2,endIdxs2,n,counter,tole,xout,yout,m)
   end subroutine Make_Band_Overlaps
   
   
   function CutFilter_Darr(arr,n,m) result(res)
      INTEGER n,m,i,j
      double precision arr(n),res(m)
      
      if(n==m)then
         res = arr
         return
      end if
      
      !always watchout for unassigned memory, kids
      ! if((FINDLOC(isNan(res),.true.,1))>1)then
         ! !allocate(bufarr(cou))
         ! stop (FINDLOC(isNan(arr),.true.,1))
      ! end if
      
      res=0d0
      j=1
      do i = 1,n
         if(arr(i)/=0d0)then
            res(j)=arr(i)
            j = j+1
         end if
      end do
      j=j+1
   end function CutFilter_Darr
   
   function CutFilter_Iarr(arr,n,m) result(res)
      INTEGER n,m,i,j,arr(n),res(m)
      
      j=1
      do i = 1,n
         if(arr(i)/=0)then
            res(j)=arr(i)
            j = j+1
         end if
      end do
   end function CutFilter_Iarr
   
   !SoS
   !https://rosettacode.org/wiki/Sorting_algorithms/Insertion_sort#Fortran
   pure subroutine Insertion_sort(a,n,comp)
      implicit none
      Integer,intent(in) :: n
      Integer,intent(out) :: comp(n)
      INTEGER :: i, j,y
      double precision,intent(inout) :: a(n)
      double precision :: x
    
      do i = 1,n
         comp(i)=i
      end do
      
      do i = 2, n
         x = a(i)
         y = comp(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            comp(j+1)=comp(j)
            j = j - 1
         end do
         a(j + 1) = x
         comp(j+1)=y
      end do
      
   end subroutine Insertion_sort 
   
   pure subroutine Insertion_sort_part(a,n,startidx,endidx,comp)
      implicit none
      Integer,intent(in) :: n,startidx,endidx
      Integer,dimension(startidx:endidx),intent(out) :: comp
      INTEGER :: i, j, k,y
      double precision,intent(inout) :: a(n)
      double precision :: x
      INTEGER :: dime
      if(startidx==endidx)return
      dime=endIdx-startidx+1
      j=1
      do i = startIdx,endIdx
         comp(i)=j
         j=j+1
      end do
      do i = startIdx+1,endidx
         x = a(i)
         y = comp(i)
         j = i - 1
         do while (j >= startidx)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            comp(j+1)=comp(j)
            j = j - 1
         end do
         a(j + 1) = x
         comp(j+1)=y
      end do
   end subroutine Insertion_sort_part
   
   pure subroutine Insertion_sort_part_C(a,n,startidx,endidx,comp)
      implicit none
      Integer,intent(in) :: n,startidx,endidx
      Integer,dimension(startidx:endidx),intent(out) :: comp
      INTEGER :: i, j, k,y
      double complex,intent(inout) :: a(n)
      double precision :: x
      INTEGER :: dime
      if(startidx==endidx)return
      dime=endIdx-startidx+1
      j=1
      do i = startIdx,endIdx
         comp(i)=j
         j=j+1
      end do
      do i = startIdx+1,endidx
         x = dble(a(i))
         y = comp(i)
         j = i - 1
         do while (j >= startidx)
            if (abs(a(j)) <= x) exit
            a(j + 1) = a(j)
            comp(j+1)=comp(j)
            j = j - 1
         end do
         a(j + 1) = x
         comp(j+1)=y
      end do
   end subroutine Insertion_sort_part_C
   
   subroutine Insertion_sort_WORK(a,n,endIdx,comp)
      implicit none
      INTEGER :: n, i, j,comp(n),endIdx,y
      double precision :: a(n), x
    
      do i = 1,endIdx
         comp(i)=i
      end do
      
      do i = 2,endIdx
         x = a(i)
         y = comp(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            comp(j+1)=comp(j)
            j = j - 1
         end do
         a(j + 1) = x
         comp(j+1)=y
      end do
      
   end subroutine Insertion_sort_WORK
   
   subroutine Insertion_sort_abs(a,n,comp)
      implicit none
      INTEGER :: n, i, j,comp(n),y
      double precision :: a(n), x
    
      do i = 1,n
         comp(i)=i
      end do
      
      do i = 2, n
         x = a(i)
         y = comp(i)
         j = i - 1
         do while (j >= 1)
            if (dabs(a(j)) <= dabs(x)) exit
            a(j + 1) = a(j)
            comp(j+1)=comp(j)
            j = j - 1
         end do
         a(j + 1) = x
         comp(j+1)=y
      end do
      
   end subroutine Insertion_sort_abs
   
   subroutine Insertion_sort_xcol(xcol,colSize,m,comp)
      implicit none
      INTEGER :: ColSize, i, j,comp(ColSize),m,y
      double precision :: xcol(colSize,m),a(ColSize), x
    
      do i = 1,ColSize
         comp(i)=i
      end do
      a = xcol(:,1)
      do i = 2, ColSize
         x = a(i)
         y = comp(i)
         j = i - 1
         do while (j >= 1)
            if (a(j) <= x) exit
            a(j + 1) = a(j)
            comp(j+1)=comp(j)
            j = j - 1
         end do
         a(j + 1) = x
         comp(j+1)=y
      end do
      call Reorder_arr_col(xcol,colSize,m,comp)
   end subroutine Insertion_sort_xcol
  
   pure subroutine Select_sort_int(arr,n,order)
      Integer,intent(in) :: n
      Integer,intent(inout) :: arr(n)
      Integer,intent(out),optional :: order(n)
      
      INTEGER i,j,bufI,idx
      if(present(order))then
         do i = 1,n
            order(i)=i
         end do
         
         do i = 1,n
            bufI=arr(i)
            idx=i
            do j = i+1,n
               if(arr(j)<bufI)then
                  bufI=arr(j)
                  idx=j
               end if
            end do
            call SwapEl_int(arr,n,i,idx)
            call SwapEl_int(order,n,i,idx)
         end do
      else
         do i = 1,n
            bufI=arr(i)
            idx=i
            do j = i+1,n
               if(arr(j)<bufI)then
                  bufI=arr(j)
                  idx=j
               end if
            end do
            call SwapEl_int(arr,n,i,idx)
         end do
      end if
   end subroutine Select_sort_int
   
   
   !trans_str - <jkm|mu|j'k'm'><j'k'm'|mu|jkm> , in Debye^2
   !nu - GHz
   function Abs_coeff(nu,N1,g1,N2,g2,trans_str,MW,absem,units) result(res)
      double precision nu,N1,N2,trans_str
      INTEGER g1,g2
      double precision res,bufr,nu_new
      logical MW,absem(2)
      character(*) units
      
      if(.not.absem(1) .and. N1/dble(g1)>N2/dble(g2))then
         res = 0d0
         return
      else if(.not.absem(2) .and. N1/dble(g1)<N2/dble(g2))then
         res = 0d0
         return
      end if
      
      nu_new=nu
      select case(units)
         case('cm-1')
            !1d-36 is conversion from Debye^2 to esu^2
            !1d-2 is conversion from m^-1 to cm^-1
            bufr=8*pi**3/(3d0*h*cc)*1d-36*1d-2
         case('l.mol-1.cm-1')
            !this is essentially a conversion from absorption cross section
            !1d-9 I have no idea (conversion from nm to m ???)
            bufr=8*pi**3*NA/(3d0*h*cc*1d3*log(10d0))*1d-36*1d-9
         case('nm2.mhz')
            !bufr=4.16231d-5
            bufr=8*pi**3/(3d0*h*cc)*1d-31*1d3
         case default
            bufr=8*pi**3/(3d0*h*cc)*1d-31*1d3
      end select
      
      if(.not.MW)bufr=bufr*4d0
      !nu seems to be anything you wish in scientific literature
      !so here in an academia sense, it will be in frequency units (GHz,Hz). Deal with it.
      res=bufr*nu_new*dabs(N1/dble(g1)-N2/dble(g2))*trans_str
   end function Abs_coeff

   function Reduce_IArr(arr,n,endidx,valu) result(res)
      INTEGER endidx,n,j,i,arr(n),res(endidx),valu
      
      j=1
      do i = 1,n
         if(arr(i)/=valu)then
            res(j)=arr(i)
            j=j+1
         end if
      end do
   end function Reduce_IArr
      
   function Reduce_filter_DArr_col(arr_col,n,colsize,filter,filteredSize) result(rescol)
      INTEGER i,j,n,colSize,filter(colSize),filteredSize
      double precision arr_col(colSize,n),rescol(filteredSize,n)
      
      j=1
      do i = 1,colsize
         if(filter(i)/=0)then
            rescol(j,:)=arr_col(i,:)
            j=j+1
         end if
      end do
   end function Reduce_filter_DArr_col
   
   !Polavarapu 2017
   function Width2FWHM_Gauss(sigma) result(W)
      double precision sigma,W
      W=2d0*sigma*sqrt(2d0*log(2d0))
   end function Width2FWHM_Gauss
   
   function FWHM2Width_Gauss(W) result(sigma)
      double precision sigma,W
      sigma=W/2d0/sqrt(2d0*log(2d0))
   end function FWHM2Width_Gauss
   
   pure function Intensity_Gauss(area,W)result(y0)
      double precision,intent(in) :: area,W
      double precision y0
      y0=area/(W*sqrt(2d0*pi))
   end function Intensity_Gauss
   
   pure function Intensity_Lorentz(area,W)result(y0)
      double precision,intent(in) :: area,W
      double precision y0
      y0=area/(W*pi)
   end function Intensity_Lorentz
   
   pure function Intensity_Voigt(area,W)result(y0)
      double precision,intent(in) :: area,W
      double precision y0
      y0=area/(W*pi)
   end function Intensity_Voigt
   
   !W = width (written as sigma, is not FWHM or HWHM or any other)
   !y0 = height
   !x0 = location of maximum
   elemental function Profile_Gauss(x,x0,W) result(y)
      double precision,intent(in) :: x,x0,W
      double precision y
      y=1d0/(W*sqrt(2d0*pi))*exp(-0.5*((x-x0)/W)**2)
   end function Profile_Gauss
   
   elemental function Profile_Lorentz(x,x0,W) result(y)
      double precision,intent(in) :: x,x0,W
      double precision y
      y=1d0/(W*pi)*((W**2)/((x-x0)**2+W**2))
   end function Profile_Lorentz
   
   pure function Profile_Voigt_Wrap(x,x0,AL,AD) result(res)
      double precision,intent(in) :: x,x0,AL,AD
      double precision res
      res = Profile_Voigt(x,x0,AL,AD,-10000,1d0,10000)
   end function Profile_Voigt_Wrap
   
   pure function Profile_Voigt(nu,nu0,AL,AD,start,step,endd) result(res)
      double precision,intent(in) :: nu,nu0,AL,AD,step
      Integer,intent(in) :: start,endd
      double precision x,y,u,res,first,second,uSecond
      INTEGER i
      x = sqrt(log(2d0))*(nu-nu0)/AD
      y = sqrt(log(2d0))*AL/AD
      res=0d0
      
      !numerical integration :>>>
      do i = start,endd
         u=start+step*(i-start)
         uSecond=u+step
         first=exp(-u**2)/(y**2+(x-u)**2)
         second=exp(-uSecond**2)/(y**2+(x-uSecond)**2)
         res=res+(first+second)/2*step
      end do
      res=res*y/pi
   end function Profile_Voigt
   
end module spectra

module dir_cose_storage
   implicit none
   
   private :: DIR_COSE_INNER
   
   !defines k2, m2
   type :: DIR_COSE_INNER
      INTEGER :: k2,m2
      character (1) :: labaxis,molaxis
      double complex,allocatable :: cube(:,:,:)
   end type
   
   !defines n2
   !k2 and m2 are in arrays
   type :: DIR_COSE_ROOT
      INTEGER :: n2
      type(DIR_COSE_INNER),allocatable :: dir_coses(:,:,:,:)
   end type
   
   interface
      pure double complex function DIR_COSE_FUNCTION(N2,K2,M2,N1,K1,M1,labaxis,molaxis)
         Integer,intent(in) :: n2,k2,m2,n1,k1,m1
         character(1),intent(in) :: molaxis,labaxis
      end function DIR_COSE_FUNCTION
   end interface   
   
end module dir_cose_storage

module strings
   implicit none
   !tried making my own String type, does not work with gfortran
   
   contains
   
   !shamelessly borrowed from: http://rosettacode.org/wiki/String_case#Fortran
   subroutine To_upper(str)
     character(*), intent(inout) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("a":"z")
           str(i:i) = achar(iachar(str(i:i))-32)
       end select
     end do 
   end subroutine To_upper
 
   subroutine To_lower(str)
     character(*), intent(inout) :: str
     integer :: i
 
     do i = 1, len(str)
       select case(str(i:i))
         case("A":"Z")
           str(i:i) = achar(iachar(str(i:i))+32)
       end select
     end do  
   end subroutine To_Lower   
   
   
   function StrAreEqual(str1,str2,caseSensitive)result(res)
      character(*),intent(in) :: str1,str2
      logical,optional,intent(in) :: caseSensitive
      
      logical res
      character(len(str1)) :: one,two
      
      ! if(len(str1)/=len(str2))then
         ! res=.false.
         ! return
      ! end if
      
      one=str1
      two=str2
      call to_upper(one)
      call to_upper(two)
      if(present(caseSensitive))then
         if(.not.caseSensitive)then
            one = (adjustl(trim(one)))
            two = (adjustl(trim(two)))
            res = one == two
         else
            one = adjustl(trim(one))
            two = adjustl(trim(two))
            res = one == two
         end if
      else
         one = (adjustl(trim(one)))
         two = (adjustl(trim(two)))
         res = one == two
      end if
      
   end function StrAreEqual
   
   function splitString(str,n,sep)result(arr)
      character(*),intent(in) :: str,sep
      integer,intent(in) :: n
      character(n),allocatable :: arr(:)
      
      integer :: m
      integer :: arr_help(50,2)
      integer :: idx1,idx2,i,idx3,leng,lengcur
      
      idx1=1
      idx2=index(str,sep)-1
      leng=len(str)
      if(idx2==-1)then
         allocate(arr(1))
         arr(1)=str
         !m=1
         return
      end if
      i=1
      do while(idx2>0)
         !lengcur=idx2-idx1+1
         arr_help(i,:)=[idx1,idx2]
         idx1=idx2+2
         idx3=idx2
         idx2=index(str(idx1:),sep)
         i=i+1
         if(idx2==0)then
            arr_help(i,:)=[idx1,leng]
            exit
         end if
         idx2=idx3+idx2
      end do
      m=i
      allocate(arr(m))
      do i = 1,m
         arr(i)=str(arr_help(i,1):arr_help(i,2))
      end do
   end function splitString
   
   function findString(arr,n,str)result(res)
      integer(4) res,i,n
      character(2) arr(n),str
      
      res=0
      do i = 1,n
         if(str==arr(i))then
            res=i
            exit
         end if
      end do
   end function findString
   
end module strings

program main
   use input
   use wigner
   use constants
   use rotor
   use spectra
   use dir_cose_storage
   use arrayOperations
   use molecule
   use strings
   
   !$ use omp_lib_kinds
   !$ use omp_lib
   
   implicit none
   
   !include "omp_lib.h"
   !===================================================
   !array(countX,countY...)
   !
   !SPH_dip(3) - spherical tensor dipole moment
   !SPH_g(3,3) - mixed spherical tensor g-factor !TODO check validity
   !SPH_quad1(3) - spherical tensor of rank 1 quadrupole moment
   !SPH_quad2(5) - spherical tensor of rank 2 quadrupole moment
   !
   !fmin - start of the spectrum, cannot be negative
   !fmax - end of the spectrum
   !resolution - constant distance between adjacent x-values of the spectrum
   !
   !press - pressure
   !temp - temperature
   !MolDens - molecular density in molecules/m^3, calculated from pressure and temperature
   !Q - rotational partition function
   !
   !startJ,endJ - Duh!
   !cou - number of possible transitions in the range(startJ,endJ,1)
   !dime_full - number of different J,Tau energy states
   !trans_arr(cou,4) - Saves the lower and upper states of transitions
   !ROT_STRS(cou) - array of rotational strengths for transitions
   !DIP_STRS(cou) - array of dipole strengths for transitions
   !TRANS_PROBS(cou) - array of CD transition probabilities !TODO
   !trans_en_arr(cou) - array of transition energies, normally in GHz
   !trans_abs(cou) - array of MW absorption coefficients, dependent on energy state populations
   !trans_probs_MW(cou) - array of MW absorption probabilities, independent of energy state populations
   !Fracpop(dime_full) - percentual populations of energy states, their sum should add up to one if exact partition function is used
   !===================================================
   
   
   double complex SPH_dip(3),SPH_g(3,3),SPH_quad1(3),SPH_quad2(5),sph_sr0(1),sph_sr1(3),sph_sr2(5),T_B0(1),T_B2(5),T_e0(1),T_e1(3),T_e2(5),T_dip1(3)
   double complex ga(3,3),gcor(3),g_qq11(3,3)
   double precision gt_own(3,3),Ra,fmin,fmax,resolution,k,time1,time2,rbuf,press,temp,Q,outTol,buf(10),intensity_cutoff,abc_tensor(3,3)
   double precision sr(3,3),S,beta
   double precision, allocatable :: TRANS_PROBS(:),trans_Fracpop(:),Fracpop(:),trans_en_arr_JPL(:)
   double precision, allocatable :: inten_JPL(:)
   INTEGER i,j,ii,jj,dime,info,rep,dime_full,m,idx,symnum
   INTEGER startJ,endJ,unitt,car,cou,spectraWorkArraySize,W3JArraySize,repaxesInt(3)
   Integer, allocatable :: Js(:),trans_arr_JPL(:,:)
   
   character(:),allocatable :: tempchar,outputfile,spr_units
   character(80) holder
   character(80) JPLFile,parfile,mainparfile,altfile
   character(80) inputfile,inputfileGau
   
   INTEGER inputLength,tempcharlen,tempcharlenchar,sta,pwidthh,tempp,idx_dot
   character(3)  profileFun
   character(1)  leftright,repaxes(3)
   character(2)  workType
   character(3)  OutUnitsStr,OutUnitsEn
   character(10) pwidthchar
   character(6) WhichPartFun
   character(7) :: catTag='      1'
   character(4) :: orderBy
   character(1) :: spinMultChar
   
   double precision pwidth
   
   

   logical :: FFPEFLAG=.true.,dopflag=.false.,presflag=.true.,approxPartFFlag=.true.
   logical :: lefthandedFlag=.false.,SkipTransFlag=.false.
   logical :: PreviousCalcFlag=.false.,MWFlag=.true.,ownPartFunFlag=.false.,switchUpAxes=.false.,dipstr_cutoff_percentual=.false.,RCDCalcGoing=.false.
   logical :: separateElNuc=.false.,gaufile=.false.,trnfile=.false.,reorder_principally=.false.,skipTransFlagLowPop=.false.
   logical :: flag_sprCM=.true.
   double precision :: LowPopTol=1d-14
   
   logical :: shift_Dipole=.false.,shift_quad=.false.,shift_gt=.false.,shift_eps=.false.
   logical :: rotate_Dipole=.false.,rotate_quad=.false.,rotate_gt=.false.,rotate_eps=.false.
   logical :: repsw_Dipole=.false.,repsw_quad=.false.,repsw_gt=.false.,repsw_eps=.false.
   double precision :: rotate_matrix(3,3)=0d0
   integer :: axes_order(3)=0

   !for error checking purposes TODO
   logical :: debugFlag=.false.,isOMP=.false.,doNothing=.false.,forceLAEvec=.false.
   INTEGER :: OMP_nthr=-1
   
   double precision dipstr_cutoff,rotstr_cutoff,kuhnpar_cutoff
   
   !tolerance for double equality comparison (usually equality to zero)
   parameter(outTol=1d-14)
   double precision, allocatable :: Ham_full(:,:),Evec_full(:,:)
   double precision, allocatable :: eval_full(:)
   double precision, allocatable :: Evec(:,:),Eval(:),Ham(:,:)
            
   !output and input
   logical :: file_molParams(10)=[.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.]
   logical :: file_transitions(10)=[.true.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.]
   logical :: file_lineSpectrum(10)=[.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.]
   logical :: file_spectrum(10)=[.true.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.]
   logical :: file_moments(10)
   logical :: file_prn(10)
   logical :: file_cat(10)
   logical :: flag_absem(2)=[.true.,.false.]
   logical :: magdip_contr(10)=[.true.,.true.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.]
   logical :: contrs(4)
   

   type :: List_2D_D
      double precision,allocatable :: arr(:,:) 
   end type
      
   type :: List_1D_D
      double precision,allocatable :: arr(:) 
   end type
   
   type :: List_1D_C
      double complex,allocatable :: arr(:) 
   end type
         
   type :: TUPLE_EVEC_SPINROT !one eigenvector for N,T,J,S state, separated into mixing states of N'
      type(List_1D_C),allocatable :: Ns(:)
      Integer,allocatable :: N_arr(:)
   end type
   
   type :: EVECS_SPINROT !all eigenvectors of given J state
      type(TUPLE_EVEC_SPINROT),allocatable :: Evecs(:)
      double precision :: J
      double complex,allocatable :: Evals(:) 
      Integer,allocatable :: N(:)
   end type
   
   type :: TRANS_JT
      INTEGER :: j1=0,t1=0,j2=0,t2=0
      double precision freq,dipstr,rotstr
      double precision en1,en2
      double precision,allocatable :: dips_contr(:),rots_contr(:)
   end type TRANS_JT
   ! type :: EVEC_LIST_SPINROT !all eigenvectors of all Js, neatly separate
      ! type(EVECS_SPINROT),allocatable :: evecs_mat(:)
   ! end type
   
   type :: TRANS_NTSJ
      INTEGER N2,N1,T2,T1
      double precision J2,J1,S2,S1
      double complex :: en2,en1
      double precision :: freq,dipstr,rotstr
   end type
   
   !type(EVEC_LIST_SPINROT) :: EVECS_MAT_SR !evec collection
   type(EVECS_SPINROT),allocatable :: EVECS_MAT_SR(:)
   type(TRANS_NTSJ),allocatable :: TRANS_NTSJ_col(:)
   type(TRANS_JT),allocatable,target :: TRANS_JT_col(:)
   
   type(List_2D_D),allocatable :: evec_full_list(:) !eigenvectors of the rigid rotor, separated for each |JT> state
   type(List_1D_D),allocatable :: eval_full_list(:) !eigenvalues --||--
   
   type(DIR_COSE_ROOT),allocatable :: DIR_COSES(:)
   
2003 format(A40,1G15.7)
2002 format(A40,1F15.7)
2456 format(A40,1F8.2)
2001 format(A40,1I10)
2469 format(A20,1F15.7,A20)
   write(6,*)'========================'
   write(6,*)'RCD program starting...'

   !initialize flags
   factFlag=.false. !dont touch >:( , not for user
   
   !allocate(evec_full_list(3))
   !allocate(evec_full_list(1)%arr(5,5))
   
   !program defaults
   !this is my spaghetti Bolognese, there are many like it but this one is mine!
   DopFlag=.false.  !include Doppler broadening of bands?
   PresFlag=.false. !include pressure broadening of bands?
   approxPartFFlag = .false. !use approximate partition function?
   lefthandedFlag=.false. !left handed asymmetric rotor representation
   FFPEFlag=.false. !write floating point exceptions?
   SkipTransFlag=.false. !not implemented, only do the transitions in the (fmin–fmax) range
   PreviousCalcFlag=.false.
   MWFlag=.false.
   ownPartFunFlag=.false.
   switchUpAxes=.false.
   
   dipstr_cutoff_percentual=.false.
   separateElNuc=.false.
   doNothing=.false.
   orderBy='QN'
   !---------------------
   !OMP VARIABLES
   !DECLARED IN SEQUENTIALLY COMPILED PROGRAM AS WELL
   isOMP=.false.
   !$ isOMP=.true.
   !$ write(6,*) 'OMP is enabled'
   OMP_nthr=1
   !$ OMP_nthr=OMP_GET_MAX_THREADS()
   !$ write(6,*) 'Count of used threads: ',OMP_nthr
   !---------------------
   
   startJ=0
   endJ=8
   unitt=8!output unit
   press=100!Pascal
   temp=20!Kelvin
   fmin=190.0d0
   fmax=220.0d0
   inputfile='G.OUT'!deprecated
   outUnitsStr='DEB'!DEB or AU
   outUnitsEn='GHz'
   
   resolution=1d-3
   

   !READ WORK FILE FROM COMMAND ARGUMENT #1
   sta=1
   call GET_COMMAND_ARGUMENT(1,mainParFile,status=sta)
   if(sta/=0)then
      mainParFile='RCD.par'
   end if
   sta=1
   call ReadFlags(mainParFile,sta)
   if(sta==0)write(6,*)'Found '//trim(mainParFile)
   
   
   !READ WORK FILE FROM COMMAND ARGUMENT #2
   sta=1
   call GET_COMMAND_ARGUMENT(2,inputfile,status=sta)
   
   if(sta/=0)then
      write(6,*)'Usage: rcd.exe [gaussian file] [work file] [[spectrum - temperature in K] [spectrum - pressure HWHM]]'
      write(6,*)'The program also looks for .par file in the same directory as the gaussian file. The .par file has the same name before the dot.'
      stop
   end if
   
   idx_dot=index(inputfile,'.',.true.)
   if(trim(inputfile(idx_dot:))=='.rcdmol')then
      call readCustomFile(inputfile)
      outputfile=inputfile(index(inputfile,'/',.true.)+1:index(inputfile,'.',.true.)-1)
   else if(trim(inputfile(idx_dot:))=='.out')then
      inputfileGau=inputfile
      call readFile(inputfileGau,outputfile,file_molParams,sta,separateElNuc,rep,lefthandedFlag)
      outputfile=inputfileGau(index(inputfileGau,'/',.true.)+1:index(inputfileGau,'.',.true.)-1)
      gaufile=.true.
   else if(trim(inputfile(idx_dot:))=='.trn')then
      call LOAD_TRN(inputfile)
      trnfile=.true.
   end if

   if(sum(abc)==0d0)write(6,*)'WARNING: Rotational constants are nil'
   if(sum(abs(dip))==0d0)write(6,*)'WARNING: Electric dipole moments are nil'
   if(sum(abs(gt))==0d0)write(6,*)'WARNING: G-Tensor components are nil'
   if(sum(abs(eps))==0d0 .and. worktype == 'SR')write(6,*)'WARNING: Spin-rotation tensor components are nil'
      
   !READ TEMPERATURE FROM COMMAND ARGUMENT #3
   sta=1
   tempchar='     '
   call GET_COMMAND_ARGUMENT(3,tempchar,status=sta)
   
   if(sta==0)then
      read(tempChar,*)temp
      write(6,*)'Setting temperature to '//trim(tempChar)//' K'
   else
      write(holder,'(I3)')NINT(temp)
      tempchar=trim(adjustL(holder))
      write(6,*)'Setting temperature to '//tempchar//' K'
   end if
   beta=1d0/(kb*temp)
   
   !READ PRESSURE HWHM FROM COMMAND ARGUMENT #4
   sta=1
   call GET_COMMAND_ARGUMENT(4,pwidthchar,status=sta)
   if(sta==0)then
      read(pwidthchar,*)pwidth
      pwidth=pwidth*1d-3
      write(6,*)'Setting Lorentzian HWHM to '//trim(pwidthchar)//' MHz'
      PresFlag=.true.
   end if
111 format(I6)  
112 format(F6.3) 
   call GET_COMMAND_ARGUMENT(5,cattag,status=sta)
   call Compute_factorials(abc(1))
   
230 format(A40)   
231 format(1F15.3,' seconds') 
  
   
   !APPEND SPIN-ROTATION CALCULATION OUTPUT WITH "_S{MULTIPLICITY}"
   if(worktype=='SR' .and. file_transitions(9))then
      write(spinMultChar,'(I1)')NINT(2*S+1)
   end if

   !JUST FOR FUN RAY PARAMETER
   Ra = Ray(abc(1),abc(2),abc(3))
   if(rep<1 .or. rep>3)then
      rep = REPRESENTATION(ra)
   end if

   !SWITCH UP THE COMPONENTS OF MOLECULAR PROPERTIES ACCORDING TO THE REPRESENTATION CHOSEN
   if(switchUpAxes .and. gaufile)then
      dip=RepAxesSwitch_V(dip,rep,lefthandedFlag)
      dip_Nuc=RepAxesSwitch_V(dip_nuc,rep,lefthandedFlag)
      dip_el=RepAxesSwitch_V(dip_el,rep,lefthandedFlag)


      gt=RepAxesSwitch_T(gt,rep,lefthandedFlag)
      ge=RepAxesSwitch_T(ge,rep,lefthandedFlag)
      gn=RepAxesSwitch_T(gn,rep,lefthandedFlag)
      quad=RepAxesSwitch_T(quad,rep,lefthandedFlag)
      quadtr=RepAxesSwitch_T(quadtr,rep,lefthandedFlag)
      eps=RepAxesSwitch_T(eps,rep,lefthandedFlag)
   end if
   
      
   
   if(shift_gt)gt = SHIFT_GT_do(gt,dip,com_shift/bohrR,a_masses*amu_2_au,natoms,r0/bohrR,r1/bohrR)
   if(shift_Dipole)dip = SHIFT_DIPOLE_do(dip,com_shift,charge)
   if(shift_quad)quad = SHIFT_QUAD_TR_do(quad,dip,com_shift,charge)
   if(shift_eps)eps = shift_eps_do(eps,com_shift) !not implemented
   
   if(rotate_dipole)dip=rotateV(dip,rotate_matrix)
   if(rotate_quad)quad=rotateT(quad,rotate_matrix)
   if(rotate_gt)gt=rotateT(gt,rotate_matrix)
   if(rotate_eps)eps=rotateT(eps,rotate_matrix)
   
   if(reorder_principally)then
      if(repsw_dipole)dip=RepAxesSwitch_V(dip,rep,lefthandedFlag)
      if(repsw_quad)quad=RepAxesSwitch_T(quad,rep,lefthandedFlag)
      if(repsw_gt)gt=RepAxesSwitch_T(gt,rep,lefthandedFlag) 
      if(repsw_eps)eps=RepAxesSwitch_T(eps,rep,lefthandedFlag) 
   end if
   
   dime_full=CALCULATE_HAMDIM(startJ,endJ)!only RR Hamiltonian dimensions
   
   !assymmetric rotor
   write(6,*)
   if(lefthandedFlag)then
      leftright='L'
   else
      leftright='R'
   end if
   repaxesInt = RepAxesSwitch_get(rep,lefthandedFlag)
   if(switchUpAxes)then
      do i = 1,3
         select case(repaxesInt(i))
            case(1)
               repaxes(i)='a'
            case(2)
               repaxes(i)='b'
            case(3)
               repaxes(i)='c'
         end select
      end do
   else
      repaxes=['a','b','c']
   end if
   write(6,2002)'ROTATIONAL CONSTANTS:               A = ',abc(1)
   write(6,2002)'                                    B = ',abc(2)
   write(6,2002)'                                    C = ',abc(3)
   write(6,2002)'RAY PARAMETER:                          ',Ra
   !rep = 3 !1
   write(6,2004)'ASSYMMETRIC TOP REPRESENTATION:         ',rep,leftright
   write(6,2112)'[x,y,z] --> [',repAxes(1),',',repaxes(2),',',repaxes(3),']'
   write(6,2001)'STARTING J:                             ',startJ
   write(6,2001)'ENDING J:                               ',endJ
   write(6,*)
   write(6,2007)'Properties in molecular axes: '
   if(separateElNuc)then
      write(6,2007)'Nuclear Dipole moment (Debye):        '
      write(6,2006)'X=              Y=              Z=           '
      write(6,2005)(dip_nuc(i)*AU2SI_debye,i=1,3)
      write(6,2007)'Electronic Dipole moment (Debye):        '
      write(6,2006)'X=              Y=              Z=           '
      write(6,2005)(dip_el(i)*AU2SI_debye,i=1,3)
   else
      write(6,2007)'Dipole moment (Debye):        '
      write(6,2006)'X=              Y=              Z=           '
      write(6,2005)(dip(i)*AU2SI_debye,i=1,3)
   end if

   write(6,*)
   write(6,2006)'X=              Y=              Z=           '
   if(separateElNuc)then
      write(6,2007)'Electronic G-Tensor (1):                 '
      write(6,2005)((ge(jj,i),i=1,3),jj=1,3)
      write(6,2007)'Nuclear G-Tensor (1):                 '
      write(6,2005)((gn(jj,i),i=1,3),jj=1,3)
   else
      write(6,2007)'G-Tensor (1):                 '
      write(6,2005)((gt(jj,i),i=1,3),jj=1,3)
   end if
   if(worktype == 'SR')then
      write(6,*)
      write(6,2007)'SR coupling tensor (GHz):                    '
      write(6,2005)((eps(jj,i),i=1,3),jj=1,3)
   end if
2004 format(A40,1I10,A1)
2007 format(A30)
2006 format(A45)
2005 format(3(F15.7,1X))
2112 format(A13,6(A1))
   if(doNothing)goto 11
   write(6,230,advance='no')'DIAGONALIZING HAMILTONIAN...            ' 
   call CPU_TIME(time1)
   abc=RepAxesSwitch_V(abc,rep,lefthandedFlag)
   if(worktype=='SR')then
      call CALCULATE_HAM_SR()
   else
      call CALCULATE_HAM()
   end if
   call CPU_TIME(time2)
   write (6,231)(time2-time1)
   write(6,*)

   if(worktype == 'RR')then !RIGID ROTOR PARTITION FUNCTION (PF)
      select case(whichPartfun)
         case('APPROX')
            Q=PARTFUN(abc(1),abc(2),abc(3),temp,symnum)
         case('EXACT')
            dime=0
            do i = startJ,endJ
               dime = dime + i*2+1
            end do

            allocate(eval_full(dime))
            car=1
            do i = startJ,endJ
               do ii = 1,2*i+1
                  eval_full(car)=eval_full_list(i-startJ+1)%arr(ii)
                  car = car + 1
               end do
            end do
            Q=ExactPARTFUN(Eval_full*1d9,dime_full,temp,symnum,endJ)
            deallocate(eval_full)
         case('SUPP')
      end select
   else if(worktype == 'SR')then !SPIN ROTATION RIGID ROTOR PARTITION FUNCTION
      select case(whichPartfun)
         case('APPROX')
            Q=PARTFUN(abc(1),abc(2),abc(3),temp,symnum)*(2*S+1)
         case('EXACT')
            Q=0d0
            do i = 1,SIZE(EVECS_MAT_SR,dim=1)
               do j = 1,SIZE(EVECS_MAT_SR(i)%Evals,dim=1)
                  Q = Q + ((EVECS_MAT_SR(i)%J)*2+1)*exp(-h*(dble(EVECS_MAT_SR(i)%Evals(j))*1d9)/(kb*temp))
               end do
            end do
            Q=Q/dble(symnum)
         case('SUPP')
      end select
   end if
   
   if(file_lineSpectrum(1).or.file_spectrum(1))then !IF SPECTRUM IS DESIRED
      write(6,2456)'TEMPERATURE (Kelvin):                   ',temp
      select case(whichPartfun)
         case('APPROX')
            write(6,2002)'PARTITION FUNCTION (Approximated):      ',Q
         case('EXACT')
            write(6,2002)'PARTITION FUNCTION (Summed):           ',Q
         case('SUPP')
            write(6,2002)'PARTITION FUNCTION (Provided):          ',Q
      end select
      if(presflag.and.file_spectrum(1))then
         write(6,2456)'PRESSURE HWHM (MHz):                    ',pwidth*1d3
      end if
      if(file_spectrum(1))then
         write(6,2456)'GRAIN (MHZ):                       ',resolution*1d3
      end if
      write(6,*)
      write(6,2456)'SPECTRUM START (GHz):                   ',fmin
      write(6,2456)'SPECTRUM END   (GHz):                   ',fmax
   end if
   
   !CALCULATE DIPOLE AND ROTATIONAL STRENGTHS
   write(6,230)'CALCULATING RCD...                           ' 
   call CPU_TIME(time1)
   call CALCULATE_RCD_ALL() 
   if(worktype=='RR')then
      deallocate(Evec_full_list,eval_full_list)
   end if
   call CPU_TIME(time2)
   write(6,231)(time2-time1)
   close(69)
      
   !CALCULATE SPECTRA
   if(.not.file_spectrum(1))goto 11
   write(6,*)
   write(6,230)'CALCULATING SPECTRA...                       ' 
   call CPU_TIME(time1)
   if(separateElNuc)then
      call CALCULATE_SPECTRUM(.false.,1,outputfile//'.spr_elel_'//trim(tempchar)//'K')
      call CALCULATE_SPECTRUM(.false.,2,outputfile//'.spr_elnuc_'//trim(tempchar)//'K')
      call CALCULATE_SPECTRUM(.false.,3,outputfile//'.spr_nucel_'//trim(tempchar)//'K')
      call CALCULATE_SPECTRUM(.false.,4,outputfile//'.spr_nucnuc_'//trim(tempchar)//'K')
      if(MWFlag)then
         call CALCULATE_SPECTRUM(.true.,1,outputfile//'.sprMW_elel_'//trim(tempchar)//'K')
         call CALCULATE_SPECTRUM(.true.,2,outputfile//'.sprMW_elnuc_'//trim(tempchar)//'K')
         call CALCULATE_SPECTRUM(.true.,3,outputfile//'.sprMW_nucel_'//trim(tempchar)//'K')
         call CALCULATE_SPECTRUM(.true.,4,outputfile//'.sprMW_nucnuc_'//trim(tempchar)//'K')
      end if
      contrs=[.true.,.true.,.true.,.true.]
      call CALCULATE_SPECTRUM_ELNUC(.false.,contrs,outputfile//'.spr_sum_'//trim(tempchar)//'K')
      call CALCULATE_SPECTRUM_ELNUC(.true.,contrs,outputfile//'.sprMW_sum_'//trim(tempchar)//'K')
   end if
   call CALCULATE_SPECTRUM(.false.)
   if(MWFlag)then
      call CALCULATE_SPECTRUM(.true.)
   end if
   call CPU_TIME(time2)
   !$ time2 = OMP_get_wtime()
   write (6,231)(time2-time1)
   
11 write(6,*)
   print *,'PROGRAM EXITED NORMALLY'
   
   if(FFPEFLAG)then
      write(6,*)'UNDERFLOW AND DENORMAL ARE OKAY, OVERFLOW AND ZERO ARE BAD'
      stop
   end if
!  -------------------------------------------------------------
!  FUNCTIONS
!  -------------------------------------------------------------
   
   contains
   
   subroutine LOAD_TRN(filee)
      character(*) filee
   end subroutine LOAD_TRN
   
   function unitMatrix(n)result(res)
      integer :: n,i
      double precision res(n,n)
      
      res=0d0
      do i = 1,n
         res(i,i)=1d0
      end do
   end function unitMatrix
   
   function vec2Mat(vec,n)result(mat)
      integer i,n
      double precision vec(n),mat(n,n)
      
      mat=0d0
      do i = 1,n
         mat(i,i)=vec(i)
      end do
   end function vec2Mat
   
   !END
   subroutine ReadCustomFile_prop(line)
      character(*),intent(inout) :: line
      character(:),allocatable :: newLine
      double precision r_help(3,1000),a_masses_help(1000),r_one(3),it(3,3),a,t,com(3),help(3),E(3,3),it_help(3,3),left(3,3),right(3,3),ang(3)
      double precision, parameter :: u=219474.63068d0,c=2.99792458d10
      double precision, allocatable :: r_other(:,:)
      character(2) asy_help(1000)
      integer,parameter :: maxIdentifierLength = 70
      character(maxIdentifierLength),allocatable :: arr(:)
      integer i,idx1,idx2,m,AN,NN,ierr,order(3),mL(2)
      character(6) sy
      
      call To_upper(line)
      if(index(line,'ABC')>0)then !rotational constants
         read(7,*)abc
         idx1=index(line,'(')+1
         idx2=index(line,')')-1
         if(idx1>0 .and. idx2>0 .and. idx2 > idx1)then
            arr = splitString(line(idx1:idx2),maxIdentifierLength,',')
            do i = 1,size(arr,dim=1)
               if(StrAreEqual(arr(i),'UNITS=GHZ'))then
                  abc=abc
               else if(StrAreEqual(arr(i),'UNITS=MHZ'))then
                  abc=abc*1d-3
               else if(StrAreEqual(arr(i),'UNITS=cm-1'))then
                  abc=abc*100d0*cc
               else if(StrAreEqual(arr(i),'UNITS=nm'))then
                  abc=cc/(abc*1d-9)*1d-9
               else if(StrAreEqual(arr(i),'UNITS=MHZ'))then
                  abc=abc*1d-3
               else if(StrAreEqual(arr(i),'UNITS=g.cm2'))then
                  abc=h/(abc/NA*1d-4*8d0*pi**2)*1d-6
               else if(StrAreEqual(arr(i),'UNITS=mc'))then
                  abc=abc*1d-3
               else if(StrAreEqual(arr(i),'UNITS=amu.a2'))then
                  abc=h/(abc/NA*1d-20*8d0*pi**2)*1d-6
               else
                  write(6,*)'ERROR: UNEXPECTED OPTION IN ROTATIONAL CONSTANTS DEFINITION'
                  stop
               end if
            end do
         else
            abc=abc
         end if
         
         abc_found=.true.
      else if(index(line,'DIP')>0)then !dipole moments
         read(7,*)dip
         idx1=index(line,'(')+1
         idx2=index(line,')')-1
         if(idx1>0 .and. idx2>0 .and. idx2 > idx1)then
            arr = splitString(line(idx1:idx2),maxIdentifierLength,',')
            do i = 1,size(arr,dim=1)
               if(StrAreEqual(arr(i),'REP_SWITCH'))then
                  !dip=RepAxesSwitch_V(dip,rep,lefthandedFlag)
                  repsw_dipole=.true.
               else if(StrAreEqual(arr(i),'CM_SHIFT'))then
                  shift_Dipole=.true.
               else if(StrAreEqual(arr(i),'ROT'))then
                  rotate_Dipole=.true.
               else if(StrAreEqual(arr(i),'UNITS=Debye'))then
                  dip=dip*si2au_debye
               else if(StrAreEqual(arr(i),'UNITS=C.M'))then
                  dip=dip*si2au_debye/cc*1d-21
               else if(StrAreEqual(arr(i),'UNITS=au'))then
                  dip=dip
               else
                  write(6,*)'ERROR: UNEXPECTED OPTION IN ELECTRIC DIPOLE MOMENT DEFINITION'
                  stop
               end if
            end do
         else 
            dip=dip*si2au_debye
         end if

         dip_found=.true.
      else if(index(line,'QUAD')>0)then !quadrupole moments
         read(7,*)quad
         read(7,*)quad
         read(7,*)quad
         idx1=index(line,'(')+1
         idx2=index(line,')')-1
         if(idx1>0 .and. idx2>0 .and. idx2 > idx1)then
            arr = splitString(line(idx1:idx2),maxIdentifierLength,',')
            do i = 1,size(arr,dim=1)
               if(StrAreEqual(arr(i),'REP_SWITCH'))then
                  !quad=RepAxesSwitch_t(quad,rep,lefthandedFlag)
                  repsw_quad=.true.
               else if(StrAreEqual(arr(i),'CM_SHIFT'))then
                  shift_quad=.true.
               else if(StrAreEqual(arr(i),'ROT'))then
                  rotate_quad=.true.
               else if(strAreEqual(arr(i),'UNITS=DEBYE.A'))then
                  quad=quad*SI2AU_debye/bohrR
               else if(strAreEqual(arr(i),'UNITS=DEBYE.AU'))then
                  quad=quad*SI2AU_debye
               else if(strAreEqual(arr(i),'UNITS=C.M2'))then
                  quad=quad*SI2AU_debye/cc*1d-21*1D10/bohrR
               else
                  write(6,*)'ERROR: UNEXPECTED OPTION IN ELECTRIC QUADRUPOLE MOMENT DEFINITION'
                  stop
               end if
            end do
         else 
            quad=quad*SI2AU_debye/bohrR
         end if
         q_found=.true.
      else if(index(line,'G_TEN')>0)then !g tensor
         read(7,*)gt(1,:)
         read(7,*)gt(2,:)
         read(7,*)gt(3,:)
         
         idx1=index(line,'(')+1
         idx2=index(line,')')-1
         if(idx1>0 .and. idx2>0 .and. idx2 > idx1)then
            arr = splitString(line(idx1:idx2),maxIdentifierLength,',')
            do i = 1,size(arr,dim=1)
               if(StrAreEqual(arr(i),'REP_SWITCH'))then
                  !gt=RepAxesSwitch_t(gt,rep,lefthandedFlag)
                  repsw_gt=.true.
               else if(StrAreEqual(arr(i),'CM_SHIFT'))then
                  shift_gt=.true.
               else if(StrAreEqual(arr(i),'ROT'))then
                  rotate_gt=.true.
               else
                  write(6,*)'ERROR: UNEXPECTED OPTION IN G-TENSOR DEFINITION'
                  stop
               end if
            end do
         else 
            gt=gt
         end if
         g_found=.true.
      else if(index(line,'SR_TEN')>0)then !spin rotation tensor
         read(7,*)eps
         read(7,*)eps
         read(7,*)eps
         
         idx1=index(line,'(')+1
         idx2=index(line,')')-1
         if(idx1>0 .and. idx2>0 .and. idx2 > idx1)then
            arr = splitString(line(idx1:idx2),maxIdentifierLength,',')
            do i = 1,size(arr,dim=1)
               if(StrAreEqual(arr(i),'REP_SWITCH'))then
                  !eps=RepAxesSwitch_t(eps,rep,lefthandedFlag)
                  repsw_eps=.true.
               else if(StrAreEqual(arr(i),'ROT'))then
                  !eps=RepAxesSwitch_t(eps,rep,lefthandedFlag)
                  rotate_eps=.true.
               else if(StrAreEqual(arr(i),'CM_SHIFT'))then
                  !eps=RepAxesSwitch_t(eps,rep,lefthandedFlag)
                  shift_eps=.true.
               else if(StrAreEqual(arr(i),'UNITS=MHZ'))then
                  eps=eps*1d-3
               else if(StrAreEqual(arr(i),'UNITS=GHZ'))then
                  eps=eps
               else
                  write(6,*)'ERROR: UNEXPECTED OPTION IN SPIN-ROTATION TENSOR DEFINITION'
                  stop
               end if
            end do
         else 
            eps=eps
         end if
         eps_found=.true.
      else if(index(line,'COORDS')>0)then !molecule coordination
         natoms=0
119      read(unit=7,fmt=*,err=121)sy,r_one
         if(sy == '      ')goto 120
         call Symbol2ANNN(trim(adjustl(sy)),AN,NN)
         if(AN==-1 .and. nn==-1)then
            backspace(7)
            backspace(7)
            goto 120
         end if
         natoms=natoms+1
         r_help(:,natoms)=r_one
         asy_help(natoms)=sy
         a_masses_help(natoms)=getIsotopeMass(AN,NN)
         goto 119
121      continue
         backspace(7)
         backspace(7)
120      continue

         allocate(r(3,natoms),a_masses(natoms),asy(natoms))
         do i = 1,natoms
            r(:,i)=r_help(:,i)
            asy(i)=asy_help(i)
            a_masses(i)=a_masses_help(i)
         end do
         r0=r
         mass=sum(a_masses)
         com=CALC_COM(natoms,a_masses,r)
         com_shift=com
         call SHIFT_COORD(natoms,r,com)
         r1=r
         com=CALC_COM(natoms,a_masses,r)
         it = CALC_INERT(natoms,a_masses,r)
         mL = maxloc(abs(it))
         axes_order=[1,2,3]
         if(mL(1)/=3 .and. ml(2)/=3)then
            axes_order=[3-ml(1),3,ml(1)]
            reorder_principally=.true.
         end if
         mL = minloc(abs(it))
         if(mL(1)/=1 .and. ml(2)/=1)then
            axes_order=[2,1,3]
            
            reorder_principally=.true.
         end if
         if(reorder_principally)then
            r_other=r
            do i = 1,natoms
               r(1,i)=r_other(axes_order(1),i)
               r(2,i)=r_other(axes_order(2),i)
               r(3,i)=r_other(axes_order(3),i)
            end do
            !TODO com_shift,dipole moments, etc.
         end if
         it_help=it
         
         call Diag(abc,it,3,ierr,.false.,real(0.1,kind(abc(1))))
         left=matmul(it_help,it)
         right=vec2Mat(abc,3)
         right=matmul(it,right)
         call Select_sort(abc,3,axes_order)
         call Reorder_arr_2d(it,3,axes_order,.true.)
         rotate_matrix=it
         !it(:,3)=-it(:,3)
         ang=ROTMAT_2_ANGLES(rotate_matrix)
         do i = 1,nAtoms
            r(:,i) = matmul(r(:,i),it)
         end do
         it_help = CALC_INERT(natoms,a_masses,r)
         ! E=unitMatrix(3)
         ! help=matmul(it_help-abc(1)*E,it(:,1))
         if(.not.abc_found)abc=h/(abc/NA*1d-20*8d0*pi**2)*1d-6
         
         r_found=.true.
      else if(strAreEqual(line,'MASS'))then
         read(7,*)mass
      else if(strAreEqual(line,'SYMM'))then
         read(7,*)symnum
      else if(strAreEqual(line,'CONST'))then
         !call setConstant()
      end if
   end subroutine ReadCustomFile_prop

   subroutine ReadCustomFile(filename)
      character(80) filename,check,units,exception
      character(120) :: curLine
      character(:),allocatable :: outfile
      INTEGER :: commentIndex
      logical FEx,useThis,switchUp,transp,altf
      double precision bufr(3,3)
      
      useThis=.true.
      abc_found=.false.
      g_found=.false.
      dip_found=.false.
      q_found=.false.
      r_Found=.false.
      eps_found=.false.
      switchUp=.false.
      transp=.false.
      
      open(7,file=filename)
      outfile=filename
      
1001  read(7,'(A120)',end=1002,err=1003)curLine
         commentIndex = index(curLine,'!')
         if(commentIndex>1)then
            call ReadCustomFile_prop(curline(1:commentIndex-1))
         else if(commentIndex==1)then
            continue
         else
            call ReadCustomFile_prop(curline)
         end if
      goto 1001
1003  stop -77
1002  continue
   end subroutine readCustomFile

   !sorta binary flags, implemented through strings :)
   subroutine ReadFlags_Binary(str,arr,n)
      Integer,intent(in) :: n
      character(n),intent(in) :: str
      logical,intent(out) :: arr(n)
      INTEGER i

      arr=.false.
      do i = 1,n
         if(str(i:i)=='1')then
            arr(i)=.true.
         else if(str(i:i)==' ')then
            return
         end if
      end do
   
   end subroutine ReadFlags_Binary
   
   !INDEX "intrinsic" BUT FOR ALLOCATABLE STRINGS
   function index_al(str,substr,back)result(idx)
      character(:),allocatable,intent(in) :: str
      character(*),intent(in) :: substr
      logical,optional :: back
      integer :: idx
      
      integer :: i,n,m
      
      n=len(str)
      m=len(substr)
     
      idx=0
      if(m>n)then
         return
      else if(m==n)then
         idx=1
         return
      end if
      if(present(back))then
         if(back)then
            do i = n-m,1,-1
               if(str(i:i+m-1) == substr)then
                  idx=i
                  return
               end if
            end do
         else
            do i = 1,n-m
               if(str(i:i+m-1) == substr)then
                  idx=i
                  return
               end if
            end do
         end if
      else
         do i = 1,n-m
            if(str(i:i+m-1) == substr)then
               idx=i
               return
            end if
         end do
      end if
   end function index_al
   
   !READ THE WORK .PAR FILE
   subroutine ReadFlags(fille,sta)
      logical lex,logbuf
      integer, parameter :: FN_length = 120
      character(FN_length) FN
      character(80) :: fille
      character(:),allocatable :: fn_help
      character(2) :: rep_char
      INTEGER line,i,j,idx,sta,idx_comment
      
      inquire(file=fille,exist=lex)
      if(.not.lex)then
         write(6,*)fille//' not found, using program defaults'
         sta=-1
         return
      end if
      sta=0
      open(unitt,file=fille)
      line=0
      flag_absem=.true.
2149  continue    
      read(unitt,2151,end=2150)FN
      
      idx_comment = index(FN,'!')
      if(idx_comment>1)then
         fn_help = FN(1:idx_comment-1)
      else if(idx_comment==1)then
         goto 2149
      else
         fn_help = FN
      end if

      if(StrAreEqual(fn_help(1:9),'CALC_TYPE'))then
         read(unitt,*)fn_help
         if(StrAreEqual(fn_help(1:2),'RR'))then
            worktype='RR'
         else if(StrAreEqual(fn_help(1:2),'SR'))then
            worktype='SR'
         end if
      end if
          
      if(StrAreEqual(fn_help,'MOLPAR_ONLY'))then
         read(unitt,*)doNothing
         line=line+1
      end if
      
      if(strAreEqual(fn_help(1:16),'GAUSSIAN_REPAXES'))then
         read(unitt,*)switchUpAxes
         line=line+1
      end if
      
      if(strAreEqual(fn_help(1:9),'SPIN_MULT'))then
         read(unitt,2151)FN
         read(FN,*)S
         S=(S-1)/2.0
         line=line+1
      end if
      
      if(strAreEqual(fn_help(1:8),'ORDER_BY'))then
         read(unitt,2151)fn_help
         if(strAreEqual(fn_help,'QN'))then
            orderBy='QN'
         else if(strAreEqual(fn_help,'FREQ'))then
            orderBy='FREQ'
         else
            write(6,*)'UNKNOWN ORDER_BY VALUE, WILL ORDER BY QUANTUM NUMBERS'
            orderBy='QN'
         end if
         line=line+1
      end if
      
      if(strAreEqual(fn_help,'MIN_DIPSTR'))then
         
         read(unitt,2151)FN
         idx=index(FN,'%')
         if(idx>0)then
            read(FN(1:idx-1),2154)dipstr_cutoff
            dipstr_cutoff=dipstr_cutoff*1d-2
            dipstr_cutoff_percentual=.true.
         else
            read(FN,2154)dipstr_cutoff
         end if
         
         line=line+1
      end if
      if(strAreEqual(fn_help,'MIN_ROTSTR'))then
         read(unitt,2154)rotstr_cutoff
         line=line+1
      end if
      if(strAreEqual(fn_help,'MIN_KUHN'))then
         read(unitt,2154)kuhnpar_cutoff
         line=line+1
      end if
      
      if(strAreEqual(fn_help,'SPR_UNITS'))then
         read(unitt,*)fn_help
         if(strAreEqual(fn_help,'l.mol-1.cm-1'))then
            spr_units=trim(adjustl(fn_help))
         else if(strAreEqual(fn_help,'cm-1'))then
            spr_units=trim(adjustl(fn_help))
         else if(strAreEqual(fn_help,'nm2.MHz'))then
            spr_units=trim(adjustl(fn_help))
         else
            write(6,*)'Unknown spectra units: ',trim(adjustl(fn_help))
            write(6,*)'Will use nm2.MHz'
            spr_units='nm2.MHz'
         end if
         line=line+1
      end if
      
      if(strAreEqual(fn_help(1:14),'ReP'))then
         read(unitt,*)rep_char
         read(rep_char(1:1),*)rep
         if(rep_char(2:2)=='r')then
            lefthandedFlag=.false.
         else if(rep_char(2:2)=='l')then
            lefthandedFlag=.true.
         end if
      end if

      
      
      
      if(fn_help(1:5)=='.molp')then
         read(unitt,2151)FN
         call ReadFlags_Binary(FN(1:10),file_molParams,10)
         line=line+1
      end if
      if(fn_help(1:4)=='.trn')then
         read(unitt,2151)FN
         call ReadFlags_Binary(FN(1:10),file_transitions,10)
         line=line+1
      end if
      if(fn_help(1:5)=='.lspr')then
         read(unitt,2151)FN
         call ReadFlags_Binary(FN(1:10),file_lineSpectrum,10)
         line=line+1
      end if
      if(fn_help(1:4)=='.spr')then
         read(unitt,2151)FN
         call ReadFlags_Binary(FN(1:10),file_spectrum,10)
         line=line+1
      end if

      if(fn_help(1:4)=='.cat')then
         read(unitt,2151)FN
         call ReadFlags_Binary(FN(1:10),file_cat,10)
         line=line+1
      end if

      if(fn_help(1:4)=='.PRN')then
         read(unitt,2151)FN
         call ReadFlags_Binary(FN(1:10),file_prn,10)
         line=line+1
      end if

      



      
      if(StrAreEqual(fn_help(1:7),'PARTFUN'))then
         read(unitt,*)fn_help
         call to_upper(fn_help)
         read(fn_help,*,err=50)q
         WhichPartFun='SUPP'
         goto 2149
50       continue
         if(fn_help == 'APPROX')WhichPartFun='APPROX'
         if(fn_help == 'EXACT')WhichPartFun='EXACT'
      end if
            
      if(StrAreEqual(fn_help,'TEMP'))then
         read(unitt,*)fn
         read(fn,*)temp
      end if
      
      if(StrAreEqual(fn_help,'GRAIN'))then
         read(unitt,*)fn
         read(fn,*)resolution
         resolution=resolution*1d-3
      end if
            
      if(StrAreEqual(fn_help,'MW_SPECTRUM'))then
         read(unitt,*)fn
         read(fn,*)MWFlag
      end if
            
      if(StrAreEqual(fn_help,'SPECTRAL_SHAPE'))then
         read(unitt,*)fn
         read(fn,*)profileFun
         if(StrAreEqual(profileFun,'LOR'))then
            presflag=.true.
            dopflag=.false.
         else if(StrAreEqual(profileFun,'GAU'))then
            presflag=.false.
            dopflag=.true.
         else if(StrAreEqual(profileFun,'VOI'))then
            presflag=.true.
            dopflag=.true.
         else
            write(6,*)'Unknown spectral profile in .par file'
            write(6,*)'Will use Lorentzian'
            presflag=.true.
            dopflag=.false.
         end if
      end if
            
      if(StrAreEqual(fn_help,'SKIP_TRANS_NOTINFREQ'))then
         read(unitt,*)fn
         read(fn,*)skipTransFlag
      end if
      
      if(StrAreEqual(fn_help,'SKIP_TRANS_LOWPOP'))then
         read(unitt,*)fn
         read(fn,*)skipTransFlagLowPop
      end if
      
      if(StrAreEqual(fn_help,'LOWPOP_TOLERANCE'))then
         read(unitt,*)fn
         read(fn,*)LowPopTol
      end if
            
      if(StrAreEqual(fn_help(1:7),'J_START'))then
         read(unitt,2153)startJ
         line=line+1
      end if
      
      if(StrAreEqual(fn_help(1:5),'J_END'))then
         read(unitt,2153)endj
         line=line+1
      end if
            
      if(StrAreEqual('HWHM',fn_help))then
         read(unitt,*)fn
         read(fn,*)pwidth
         pwidth=pwidth*1d-3
         presflag=.true.
      end if
      
      if(StrAreEqual('FREQ_RANGE',fn_help))then
         read(unitt,*)fn
         idx=index(fn,'-')
         read(fn(1:idx-1),*)fmin
         read(fn(idx+1:),*)fmax
      end if
            
      
      goto 2149
2150  close(unitt)
2151  format(A80)    
2152  format(1L)
2153  format(I8)
2154  format(F15.8)
2155  format(3X,1L)
2156  format(3X,F15.8)
   end subroutine ReadFlags
         
   !CALCULATE BAND SPECTRUM, SPIN-ROTATION RR
   subroutine CALCULATE_SPECTRUM_SR(MW)
      logical MW
      type(TRANS_NTSJ) tr
      INTEGER finalSize,g1,g2
      double precision transstr,abscoeff,dwidth,frpop1,frpop2
      double precision,allocatable :: x(:),y(:),bufr(:)
      character(:),allocatable :: filee
      
      finalSize=NINT((fmax-fmin)/resolution)
      allocate(x(finalSize),y(finalSize),bufr(finalSize))
      
      do i = 1,finalSize
         x(i)=fmin+resolution*(i-1)
         y(i)=0d0
      end do
      
      if(presflag.and.dopflag)then
         do i = 1,car-1
         end do
      else if(presflag)then
         do i = 1,car-1
            tr=TRANS_NTSJ_col(i)
            if(tr%freq==0)cycle
            if(MW)then
               transstr=tr%dipstr
               if(file_spectrum(2).and.transstr<dipstr_cutoff)cycle
            else
               transstr=tr%rotstr
               if(file_spectrum(2).and.abs(transstr)<rotstr_cutoff)cycle
            end if
            if(file_spectrum(3))then
               frpop2=(2*tr%J2+1)*exp(-h*dble(tr%en2)*1d9/(kb*temp))/Q
               frpop1=(2*tr%J1+1)*exp(-h*dble(tr%en1)*1d9/(kb*temp))/Q
               g2=NINT((tr%J2)*2+1)
               g1=NINT((tr%J1)*2+1)
            else
               frpop1=1
               frpop2=0
               g1=1
               g2=1
            end if            
            abscoeff=Abs_coeff(abs(tr%freq),frpop2,g2,frpop1,g1,transstr*(AU2SI_debye)**2,MW,flag_absem,spr_units)
            bufr=Profile_Lorentz(x,abs(tr%freq),pwidth)
            y=y+abscoeff*bufr
         end do
      else if(dopflag)then
         do i = 1,car-1
            tr=TRANS_NTSJ_col(i)
            if(tr%freq==0)cycle
            if(MW)then
               transstr=tr%dipstr
               if(file_spectrum(2).and.transstr<dipstr_cutoff)cycle
            else
               transstr=tr%rotstr
               if(file_spectrum(2).and.abs(transstr)<rotstr_cutoff)cycle
            end if
            if(file_spectrum(3))then
               frpop2=(2*tr%J2+1)*exp(-h*dble(tr%en2)*1d9/(kb*temp))/Q
               frpop1=(2*tr%J1+1)*exp(-h*dble(tr%en1)*1d9/(kb*temp))/Q
               g2=NINT((tr%J2)*2+1)
               g1=NINT((tr%J1)*2+1)
            else
               frpop1=1
               frpop2=0
               g1=1
               g2=1
            end if            
            dwidth=Doppler_width(abs(tr%freq),temp,mass*1D-3/NA)
            bufr=Profile_Gauss(x,abs(tr%freq),dwidth)
            abscoeff=Abs_coeff(abs(tr%freq),frpop2,g2,frpop1,g1,transstr*(AU2SI_debye)**2,MW,flag_absem,spr_units)
            y=y+abscoeff*bufr
         end do
      end if
      if(MW)then
         if(file_spectrum(4))then
            filee=outputfile//'_S'//spinMultChar//'.sprMW_'//trim(tempChar)//'K'
         else
            filee=outputfile//'.sprMW'
         end if
      else
         if(file_spectrum(4))then
            filee=outputfile//'_S'//spinMultChar//'.spr_'//trim(tempChar)//'K'
         else
            filee=outputfile//'.spr'
         end if
      end if
      open(unitt+1,file=filee)
      call Write_Spectrum(x,y,finalSize,unitt+1)
      if(.not.RCDCalcGoing)then
         write(6,*)'WRITTEN '//filee
      end if
      close(unitt+1)
   end subroutine CALCULATE_SPECTRUM_SR
   
   !CALCULATE BAND SPECTRUM, RIGID ROTOR
   subroutine CALCULATE_SPECTRUM(MW,contr,filename)
      logical,intent(in) :: MW
      character(*),optional,intent(in) :: filename
      Integer,optional,intent(in) :: contr
      
      double precision, parameter :: tol = 1d-6
      double precision,allocatable :: x_col(:,:),y_col(:,:),x_col2(:,:),y_col2(:,:)
      double precision,allocatable :: x(:),y(:),xy(:,:),y_thr(:,:)
      Integer,allocatable :: orderXY(:)
      double precision,allocatable :: bufR(:),stupidFortranBuffer(:)
      double precision dwidth,inten,curY,curX,abscoeff,stuff(cou),dmax,dips,rots,gk
      double precision frpop1,frpop2
      Integer,allocatable :: trans_arr_(:,:),order(:),orderR(:),endIdxs(:),colOrder(:),endIdxs2(:)
      INTEGER idx1,idx2,colIdx,filtered(cou),finalSize,g1,g2
      INTEGER i,j,numX,alArrSize,curIndex,bufI,counter,filteredCount,zeroIDX
      character(3) profile
      character(:),allocatable :: filee
      type(trans_jt) tr
      
            
      if(.not.(PresFlag.or.dopflag).or..not.file_spectrum(1))return
      if(worktype == 'SR')then
         call CALCULATE_SPECTRUM_SR(MW)
         return
      end if
      
      
      
      finalSize=NINT((fmax-fmin)/resolution)
      allocate(x(finalSize),y(finalSize),bufr(finalSize))
      y=0d0
      do i = 1,finalSize
         x(i)=fmin+resolution*i
      end do
      
      
      
      zeroIDX=cou+1
      do i = 1,cou
         tr=trans_jt_col(i)
         if(tr%j1 == 0 .and. tr%j2 == 0)then
            zeroIDX = i
            exit
         end if
      end do
      if(presflag.and.dopflag)then
         do i = 1,zeroIDX-1
            tr=trans_jt_col(i)
            rots=tr%rotstr
            dips=tr%dipstr
            if(present(contr))then
               if(MW)then
                  inten=tr%dips_contr(contr)
               else
                  inten=tr%rots_contr(contr)
               end if
            else
               if(MW)then
                  inten=dips
               else
                  inten=rots
               end if
            end if
            if(file_spectrum(2).and.((abs(rots)<rotstr_cutoff) .or. (dips<dipstr_cutoff)))cycle
            if(dips>tol)then
               gk=4d0*rots/dips
            else
               gk=0d0
            end if
            if(file_spectrum(2).and.(abs(gk)<kuhnpar_cutoff))cycle
            if(file_spectrum(3))then
               idx1=Calc_Frac_index(tr%j1,tr%t1)
               idx2=Calc_Frac_index(tr%j2,tr%t2)
               frpop1=Fracpop(idx1)
               frpop2=Fracpop(idx2)
               g1=tr%j1*2+1
               g2=tr%j2*2+1
            else
               frpop1=1
               frpop2=0
               g1=1
               g2=1
            end if
            dwidth=Doppler_width(tr%freq,temp,mass*1D-3/NA)
            abscoeff=Abs_coeff(tr%freq,frpop1,g1,frpop2,g2,inten*(AU2SI_debye)**2,MW,flag_absem,spr_units)
            do j=1,finalSize
               y(j)=y(j)+abscoeff*Profile_Voigt(x(j),tr%freq,pwidth,dwidth,-10,5d-2,10)
            end do
         end do
      else if(presflag)then
         do i = 1,zeroIDX-1
            tr=trans_jt_col(i)
            rots=tr%rotstr
            dips=tr%dipstr
            if(present(contr))then
               if(MW)then
                  inten=tr%dips_contr(contr)
               else
                  inten=tr%rots_contr(contr)
               end if
            else
               if(MW)then
                  inten=dips
               else
                  inten=rots
               end if
            end if
            if(file_spectrum(2).and.((abs(rots)<rotstr_cutoff) .or. (dips<dipstr_cutoff)))cycle
            if(file_spectrum(3))then
               idx1=Calc_Frac_index(tr%j1,tr%t1)
               idx2=Calc_Frac_index(tr%j2,tr%t2)
               frpop1=Fracpop(idx1)
               frpop2=Fracpop(idx2)
               g1=tr%j1*2+1
               g2=tr%j2*2+1
            else
               frpop1=1
               frpop2=0
               g1=1
               g2=1
            end if
            abscoeff=Abs_coeff(tr%freq,frpop1,g1,frpop2,g2,inten*(AU2SI_debye)**2,MW,flag_absem,spr_units)
            bufr=Profile_Lorentz(x,tr%freq,pwidth)
            y=y+abscoeff*bufr       
         end do
      else if(dopflag)then
         do i = 1,zeroIDX-1
            tr=trans_jt_col(i)
            rots=tr%rotstr
            dips=tr%dipstr
            if(present(contr))then
               if(MW)then
                  inten=tr%dips_contr(contr)
               else
                  inten=tr%rots_contr(contr)
               end if
            else
               if(MW)then
                  inten=dips
               else
                  inten=rots
               end if
            end if
            if(file_spectrum(2).and.((abs(rots)<rotstr_cutoff) .or. (dips<dipstr_cutoff)))cycle
            if(file_spectrum(3))then
               idx1=Calc_Frac_index(tr%j1,tr%t1)
               idx2=Calc_Frac_index(tr%j2,tr%t2)
               frpop1=Fracpop(idx1)
               frpop2=Fracpop(idx2)
               g1=tr%j1*2+1
               g2=tr%j2*2+1
            else
               frpop1=1
               frpop2=0
               g1=1
               g2=1
            end if
            
            dwidth=Doppler_width(tr%freq,temp,mass*1D-3/NA)
            abscoeff=Abs_coeff(tr%freq,frpop1,g1,frpop2,g2,inten*(AU2SI_debye)**2,MW,flag_absem,spr_units)
            bufr=Profile_Gauss(x,tr%freq,dwidth)
            y=y+abscoeff*bufr
         end do
      end if
      
      if(present(filename))then
         filee=filename
      else
         if(MW)then
            if(file_spectrum(4))then
               filee=outputfile//'_S'//spinMultChar//'.sprMW_'//trim(tempChar)//'K'
            else
               filee=outputfile//'.sprMW'
            end if
         else
            if(file_spectrum(4))then
               filee=outputfile//'_S'//spinMultChar//'.spr_'//trim(tempChar)//'K'
            else
               filee=outputfile//'.spr'
            end if
         end if
      end if
      open(unitt,file=filee)
      call Write_Spectrum(x,y,finalSize,unitt)
      close(unitt)
230   format(A45)   
      write(6,*)'WRITTEN '//filee
      deallocate(x,y)
      return
      
   end subroutine CALCULATE_SPECTRUM
      
   !CALCULATE ELECTRONIC AND NUCLEAR CONTRIBUTIONS BAND SPECTRA
   subroutine CALCULATE_SPECTRUM_ELNUC(MW,contr,filename)
      logical,intent(in) :: MW
      character(*),optional,intent(in) :: filename
      logical,intent(in) :: contr(4)
   
      double precision,allocatable :: x_col(:,:),y_col(:,:),x_col2(:,:),y_col2(:,:)
      double precision,allocatable :: x(:),y(:),xy(:,:),y_thr(:,:)
      Integer,allocatable :: orderXY(:)
      double precision,allocatable :: bufR(:)
      double precision dwidth,inten,curY,curX,abscoeff,stuff(cou),dmax,dips,rots
      double precision frpop1,frpop2
      Integer,allocatable :: trans_arr_(:,:),order(:),orderR(:),endIdxs(:),colOrder(:),endIdxs2(:)
      INTEGER idx1,idx2,colIdx,filtered(cou),finalSize
      INTEGER g1,g2
      INTEGER i,j,k,numX,alArrSize,curIndex,bufI,counter,filteredCount,zeroIDX
      character(3) profile
      character(:),allocatable :: filee
      type(trans_jt) tr
      
            
      if(.not.(PresFlag.or.dopflag).or..not.file_spectrum(1))return
      if(worktype == 'SR')then
         call CALCULATE_SPECTRUM_SR(MW)
         return
      end if
      
      
      
      finalSize=NINT((fmax-fmin)/resolution)
      allocate(x(finalSize),y(finalSize),bufr(finalSize))
      y=0d0
      do i = 1,finalSize
         x(i)=fmin+resolution*i
      end do
      
      
      
      zeroIDX=cou+1
      do i = 1,cou
         tr=trans_jt_col(i)
         if(tr%j1 == 0 .and. tr%j2 == 0)then
            zeroIDX = i
            exit
         end if
      end do
      do k=1,4
         if(.not.contr(k))cycle
         if(presflag.and.dopflag)then
            do i = 1,zeroIDX-1
               tr=trans_jt_col(i)
               rots=tr%rotstr
               dips=tr%dipstr
               if(MW)then
                  inten=tr%dips_contr(k)
               else
                  inten=tr%rots_contr(k)
               end if
               if(((abs(rots)<rotstr_cutoff) .or. ((dips<dipstr_cutoff))).and.file_spectrum(2))cycle
               if(file_spectrum(3))then
                  idx1=Calc_Frac_index(tr%j1,tr%t1)
                  idx2=Calc_Frac_index(tr%j2,tr%t2)
                  frpop1=Fracpop(idx1)
                  frpop2=Fracpop(idx2)
                  g1=tr%j1*2+1
                  g2=tr%j2*2+1
               else
                  frpop1=1
                  frpop2=0
                  g1=1
                  g2=1
               end if
               dwidth=Doppler_width(tr%freq,temp,mass*1D-3/NA)
               abscoeff=Abs_coeff(tr%freq,frpop1,g1,frpop2,g2,inten*(AU2SI_debye)**2,MW,flag_absem,spr_units)
               do j=1,finalSize
                  y(j)=y(j)+abscoeff*Profile_Voigt(x(j),tr%freq,pwidth,dwidth,-10,5d-2,10)
               end do
            end do
         else if(presflag)then
            do i = 1,zeroIDX-1
               tr=trans_jt_col(i)
               rots=tr%rotstr
               dips=tr%dipstr
               if(MW)then
                  inten=tr%dips_contr(k)
               else
                  inten=tr%rots_contr(k)
               end if
               if(((abs(rots)<rotstr_cutoff) .or. ((dips<dipstr_cutoff))).and.file_spectrum(2))cycle
               if(file_spectrum(3))then
                  idx1=Calc_Frac_index(tr%j1,tr%t1)
                  idx2=Calc_Frac_index(tr%j2,tr%t2)
                  frpop1=Fracpop(idx1)
                  frpop2=Fracpop(idx2)
                  g1=tr%j1*2+1
                  g2=tr%j2*2+1
               else
                  frpop1=1
                  frpop2=0
                  g1=1
                  g2=1
               end if
               abscoeff=Abs_coeff(tr%freq,frpop1,g1,frpop2,g2,inten*(AU2SI_debye)**2,MW,flag_absem,spr_units)
               bufr=Profile_Lorentz(x,tr%freq,pwidth)
               y=y+abscoeff*bufr       
            end do
         else if(dopflag)then
            do i = 1,zeroIDX-1
               tr=trans_jt_col(i)
               rots=tr%rotstr
               dips=tr%dipstr
               if(MW)then
                  inten=tr%dips_contr(k)
               else
                  inten=tr%rots_contr(k)
               end if
               if(((abs(rots)<rotstr_cutoff) .or. ((dips<dipstr_cutoff))).and.file_spectrum(2))cycle
               if(file_spectrum(3))then
                  idx1=Calc_Frac_index(tr%j1,tr%t1)
                  idx2=Calc_Frac_index(tr%j2,tr%t2)
                  frpop1=Fracpop(idx1)
                  frpop2=Fracpop(idx2)
                  g1=tr%j1*2+1
                  g2=tr%j2*2+1
               else
                  frpop1=1
                  frpop2=0
                  g1=1
                  g2=1
               end if
               dwidth=Doppler_width(tr%freq,temp,mass*1D-3/NA)
               abscoeff=Abs_coeff(tr%freq,frpop1,g1,frpop2,g2,inten*(AU2SI_debye)**2,MW,flag_absem,spr_units)
               bufr=Profile_Gauss(x,tr%freq,dwidth)
               y=y+abscoeff*bufr
            end do
         end if
      end do
      
      if(present(filename))then
         filee=filename
      else
         if(MW)then
            filee=outputfile//'.sprMW_'//trim(tempChar)//'K'
         else
            filee=outputfile//'.spr_'//trim(tempChar)//'K'
         end if
      end if
      open(unitt,file=filee)
      call Write_Spectrum(x,y,finalSize,unitt)
      close(unitt)
230   format(A45)   
      write(6,*)'WRITTEN '//filee
      deallocate(x,y)
      return
      
   end subroutine CALCULATE_SPECTRUM_ELNUC
      
   subroutine Write_Spectrum(x,y,n,unit)
      INTEGER n,i,unit
      double precision x(n),y(n)
      
      do i = 1,n
         write(unit,1325)x(i),y(i)
      end do
      
1325  format(2(E18.10E3,1x)) 

      
      !call ARR2D_2_FILE_spec(xy,n,2,idx,unitt)
   end subroutine Write_Spectrum
      
   !writes the lines to that god awful JPL .cat format, hooray for and intensity units
   subroutine ToJPLCAT()
      INTEGER i,Km1_1,Km1_2,K1_1,K1_2,j1,j2,t1,t2,bufI(3),idx1,idx2,orderR(cou)
      double precision jplInten,trans_en_arr_(cou),fp1,fp2
      type(trans_jt) tr

     
      open(100,file=outputfile//'.cat')
      
      trans_en_arr_=trans_jt_col(:)%freq
      call Insertion_sort(trans_en_arr_,cou,orderR)
      
      
      do i=1,car
         tr=trans_jt_col(i)
         j1=tr%j1
         j2=tr%j2
         if(j1==0 .and. j2==0)exit
         t1=tr%t1
         t2=tr%t2
         idx1=Calc_Frac_index(j1,t1)
         idx2=Calc_Frac_index(j2,t2)
         bufI=OldToNewNotation(j1,t1)
         Km1_1=bufI(2)
         K1_1=bufI(3)
         bufI=OldToNewNotation(j2,t2)
         Km1_2=bufI(2)
         K1_2=bufI(3)
         fp1=Fracpop(idx1)
         fp2=Fracpop(idx2)
         
         if(fp1>fp2)then
            jplInten=4.16231D-5*tr%freq*1D3*tr%dipstr*(2*j1+1)*(fp1-fp2)
         else
            jplInten=4.16231D-5*tr%freq*1D3*tr%dipstr*(2*j1+1)*(fp2-fp1)
         end if
         write(100,2000)tr%freq*1D3,0d0,log10(jplInten),3,0d0,j2*2+1,0,1404,j1,Km1_1,K1_1,0,0,0,j2,Km1_2,K1_2,0,0,0
      end do
      
2000  format(F13.4,F8.4,F8.4,I2,F10.4,I3,I7,I4,6I2,6I2)
      close(100)
   end subroutine ToJPLCAT
      
   !analytical functions derived by salzmann, kept here for reasons
   !http://dx.doi.org/10.1063/1.474598
   subroutine SALZMAN_COMPARISON()
      double precision Ds(18),Rs(18),buf
      INTEGER transs(18,4),i
      
      Ds=0d0
      Rs=0d0
      transs=0
      
      Rs(1)=1d0/2d0*dip(3)*(gt(1,2)-gt(2,1))
      Ds(1)=dip(3)**2
      transs(1,1)=0
      transs(1,2)=0
      transs(1,3)=1
      transs(1,4)=-1
      
      Rs(2)=1d0/2d0*dip(1)*(gt(2,3)-gt(3,2))
      Ds(2)=dip(1)**2
      transs(2,1)=0
      transs(2,2)=0
      transs(2,3)=1
      transs(2,4)=0
      
      Rs(3)=1d0/2d0*dip(2)*(gt(3,1)-gt(1,3))
      Ds(3)=dip(2)**2
      transs(3,1)=0
      transs(3,2)=0
      transs(3,3)=1
      transs(3,4)=1
      
      Rs(4)=-1d0/4d0*dip(2)*(gt(3,1)+gt(1,3))
      Ds(4)=dip(2)**2/2d0
      transs(4,1)=1
      transs(4,2)=-1
      transs(4,3)=1
      transs(4,4)=0
      
      Rs(5)=1d0/4d0*dip(1)*(gt(2,3)+gt(3,2))
      Ds(5)=dip(1)**2/2d0
      transs(5,1)=1
      transs(5,2)=-1
      transs(5,3)=1
      transs(5,4)=1
      
      Rs(6)=-1d0/4d0*dip(3)*(gt(1,2)+gt(2,1))
      Ds(6)=dip(3)**2/2d0
      transs(6,1)=1
      transs(6,2)=0
      transs(6,3)=1
      transs(6,4)=1
      
      Rs(7)=1d0/4d0*dip(1)*(gt(2,3)-3d0*gt(3,2))
      Ds(7)=dip(1)**2/2d0
      transs(7,1)=1
      transs(7,2)=-1
      transs(7,3)=2
      transs(7,4)=-1
      
      Rs(8)=1d0/4d0*dip(2)*(3d0*gt(3,1)-gt(1,3))
      Ds(8)=dip(2)**2/2d0
      transs(8,1)=1
      transs(8,2)=-1
      transs(8,3)=2
      transs(8,4)=0
      
      Rs(9)=1d0/4d0*dip(3)*(3d0*gt(1,2)-gt(2,1))
      Ds(9)=dip(3)**2/2d0
      transs(9,1)=1
      transs(9,2)=0
      transs(9,3)=2
      transs(9,4)=-1
      
      Rs(10)=1d0/4d0*dip(2)*(gt(3,1)-3d0*gt(1,3))
      Ds(10)=dip(2)**2/2d0
      transs(10,1)=1
      transs(10,2)=0
      transs(10,3)=2
      transs(10,4)=1
      
      Rs(11)=1d0/4d0*dip(3)*(gt(1,2)-3d0*gt(2,1))
      Ds(11)=dip(3)**2/2d0
      transs(11,1)=1
      transs(11,2)=1
      transs(11,3)=2
      transs(11,4)=0
      
      Rs(12)=1d0/4d0*dip(1)*(3d0*gt(2,3)-gt(3,2))
      Ds(12)=dip(1)**2/2d0
      transs(12,1)=1
      transs(12,2)=1
      transs(12,3)=2
      transs(12,4)=1
      
      Rs(13)=-1d0/4d0*dip(3)*(gt(1,2)+gt(2,1))
      Ds(13)=dip(3)**2/6d0
      transs(13,1)=2
      transs(13,2)=-1
      transs(13,3)=2
      transs(13,4)=0
      
      Rs(14)=1d0/4d0*dip(1)*(gt(2,3)+gt(3,2))
      Ds(14)=dip(1)**2/6d0
      transs(14,1)=2
      transs(14,2)=-1
      transs(14,3)=2
      transs(14,4)=1
      
      Rs(15)=-1d0/4d0*dip(2)*(gt(3,1)+gt(1,3))
      Ds(15)=dip(2)**2/6d0
      transs(15,1)=2
      transs(15,2)=0
      transs(15,3)=2
      transs(15,4)=1
      
      Rs(16)=1d0/2d0*dip(2)*(gt(3,1)-gt(1,3))
      Ds(16)=dip(2)**2/3d0
      transs(16,1)=2
      transs(16,2)=-1
      transs(16,3)=3
      transs(16,4)=0
      
      Rs(17)=1d0/2d0*dip(1)*(gt(2,3)-gt(3,2))
      Ds(17)=dip(1)**2/3d0
      transs(17,1)=2
      transs(17,2)=0
      transs(17,3)=3
      transs(17,4)=0
      
      Rs(18)=1d0/2d0*dip(3)*(gt(1,2)-gt(2,1))
      Ds(18)=dip(3)**2/3d0
      transs(18,1)=2
      transs(18,2)=1
      transs(18,3)=3
      transs(18,4)=0
      
      Rs=nucm_au*Rs
      
      open(unitt,file='COMPARISON.TAB')
      write(unitt,663)
      write(unitt,660)
660   format(60(1h-))
      
      do i = 1,18
      if(.not.(abs(Ds(i))<0.0d0))then
         buf = 4.0d0*Rs(i)/Ds(i)
      else
         buf = 0.0d0
      end if
      
      
      write(unitt,661)i,0d0,Ds(i),Rs(i),buf,transs(i,1),transs(i,2),transs(i,3),transs(i,4),0d0
      end do
661   format(i7,f15.5,3(' ',E15.7E3),4i3,'    ',E15.7E3)
663   format('                  freq              D              R           4R/D J1 T1 J2 T2   Probab',/,'transition         GHz           au^2           au^2')
      close(unitt)
   end subroutine SALZMAN_COMPARISON
   
   !shorten precision, originally made for comparison purposes with SPCAT, kept here for reasons
   elemental function ShortPres(num,pres)result(res)
      double precision,intent(in) :: num
      double precision res,help
      Integer,intent(in) :: pres !amount of accurate digits
      INTEGER newPres
      
      !if(abs(num<1d-1))
      help = num
      newPres=0
      do while(abs(help)>1d0)
         help=help/10
         newPres=newPres+1
      end do
      res=dble(NINT(num*(10**(pres-newPres)),kind=8))/(10**(pres-newPres))
      !res=help*10**newpres
   end function ShortPres
      
   !calculates dimensions of the full hamiltonian in J (or N if coupled) range
   function CALCULATE_HAMDIM(LJ,HJ) result (res)
      Integer, intent(in) :: LJ,HJ
      INTEGER res,arr(HJ-LJ+1),i,leng
   
      leng=HJ-LJ+1
      call HAM_DIM_SEQUENCE(LJ,HJ,arr)
      res = 0
      do i = 1,leng
         res = res + arr(i)
      end do
   end function CALCULATE_HAMDIM
   
   !CALCULATE AND DIAGONALIZE (read: Solve the Schrodinger equation) HAMILTONIAN FOR EACH J, RIGID ROTOR
   !THE WAVEFUNCTION IS |J T M >
   subroutine CALCULATE_HAM()
      INTEGER ii,idx,idx2,k1,k2,i
      double precision,allocatable :: E(:,:),rbuf(:)
      Integer,allocatable :: order(:)
      logical success

      allocate(evec_full_list(endJ-startJ+1),eval_full_list(endJ-startJ+1))
      car = 1
      
      do i = startJ,endJ
         dime = i*2+1
         allocate(Ham(dime,dime),Eval(dime),Evec(dime,dime),E(dime,dime),rbuf(dime),order(dime))
         E=0.0d0
         Ham = 0.0d0
         evec = 0.0d0
         eval = 0.0d0
         Ham=HAM_RR_OWN(i,i,abc(1),abc(2),abc(3))
         Evec=Ham
         call DIAG(Eval,Evec,dime,info,.false.,dblTol)
         call Insertion_sort(eval,dime,order)
         call Reorder_Evec(evec,dime,1,dime,order)

         if(info/=0)stop 'CALCULATE_HAM: ERROR AT DIAGONALIZATION'

         
         allocate(evec_full_list(car)%arr(dime,dime),eval_full_list(car)%arr(dime))
         evec_full_list(car)%arr = evec
         eval_full_list(car)%arr = eval
         deallocate(Ham,Eval,Evec,E,rbuf,order)
         car = car + 1
      end do
   end subroutine CALCULATE_HAM
         
   function WHICH_Ns(J,S)result(res)
      double precision J,S
      INTEGER i,Namount,N
      Integer,allocatable :: res(:)
      
      Namount=NUM_OF_Ns(J,S)
      allocate(res(Namount))
      
      N=NINT(J+S)
      i=1
      do while(dble(N)>=J-S .and. dble(N)<=J+S)
         res(i)=N
         N=N-1
         i=i+1
      end do
      !call Select_sort_int(res,Namount)
      call Reverse_arr_int(res,Namount)
   end function WHICH_Ns
   
   function NUM_OF_Ns(J,S)result(res)
      double precision J,S
      INTEGER i,res,N
      
      res=0
      N=NINT(J+S)
      do while(dble(N)>=J-S .and. dble(N)<=J+S)
         res=res+1
         N=N-1
      end do
   end function NUM_OF_Ns
   
   function I2str(i)result(str)
      integer :: i
      character(4) :: str
      write(str,'(1I4)')i
   end function I2str
   
   !CALCULATE AND DIAGONALIZE HAMILTONIAN FOR THE SPIN ROTATION CASE
   !THE WAVE FUNCTION IS: |N K S J M >
   subroutine CALCULATE_HAM_SR()
      INTEGER hamdim,N2,N1,K2,K1,startIDX,startIDXc,startIDXr,i,car,k,endjCustom,bufi,idx
      Integer,allocatable :: N_arr(:),K2_arr(:),K1_arr(:),N2_arr(:),N1_arr(:),N_arr_buf(:)
      INTEGER mult,Jamount,N,idx1,idx2
      double precision J,deter
      double complex,allocatable :: ham(:,:),eval(:),evec(:,:),eval_buf(:),evec_buf(:,:)
      double complex ,allocatable :: eval_help(:),evec_help(:,:)
      logical isOdd
      
      type :: N_Js
         INTEGER N
         double precision,allocatable :: Js(:)
      end type N_Js
      
      type(N_Js),allocatable :: N_Js_col(:)
      
      N2=0 !lower
      N1=1 !higher, pretty sure this is the convention :)))
      car=1
      mult=nint(2*s+1)
      !startIDX=1
      !allocate(evec_full_list(endJ+1),eval_full_list(endJ+1))
      if(WhichPartFun=='EXACT')then
         endjCustom=120 !any more than 120 and you are going to need REAL(10) for the factorials
      else
         endjCustom = endj
      end if
      allocate(EVECS_MAT_SR(endjCustom+1))
      info=1
      
      ! do i = 0,mult-1,1
         ! N_arr(i+1)=i
      ! end do
      
      
      if(MOD(NINT(2.0*S),2)==0)then !S is even (1, 2, 3, ...)
         isOdd=.false.
      else !S is odd (1/2, 3/2, 5/2, ...)
         isOdd=.true.
      end if
      
      car=1
      if(S==0d0)write(6,*)'WARNING: Spin is set to nil'
      do i = 0,endjCustom
         
         J=S+i
         !write(6,*)J
         N_arr=WHICH_Ns(J,S)
         hamdim = sum(N_arr*2+1)

         allocate(ham(hamdim,hamdim),eval(hamdim),evec(hamdim,hamdim))
         ham=0d0
         allocate(K2_arr(hamdim),K1_arr(hamdim),N2_arr(hamdim),N1_arr(hamdim))
         
         !prepare K2 and K1 indices
         idx=1
         do ii = 1,size(N_arr,dim=1)
            !if(N_arr(ii)==0)cycle
            do K2 = -N_arr(ii),N_arr(ii),1
               K2_arr(idx)=K2
               N2_arr(idx)=N_arr(ii)
               idx=idx+1
            end do
         end do
         K1_arr=K2_arr
         N1_arr=N2_arr
         

         do ii = 1,hamdim !row
            N2=N2_arr(ii)
            K2=K2_arr(ii)
            do jj = 1,hamdim !column
               K1=K1_arr(jj)
               N1=N1_arr(jj)
               
               ham(ii,jj)=HAM_RR_OWN_MATEL_NKSJM(abc,N2,K2,S,J,-J,N1,K1,S,J,-J)+HAM_SR_OWN_MATEL_NKSJM(eps,N2,K2,S,J,-J,N1,K1,S,J,-J,.true.)
            end do
         end do
         
         
         !ham(1,2)=0d0
         do ii = 1,hamdim
            do jj = 1,hamdim
               if(dble(ham(ii,jj)-conjg(ham(jj,ii)))>=1d-7)then
                  write(6,*) 'Warning: Hamiltonian not hermitian, J = '//I2STR(i)
                  goto 2000
               end if
            end do
         end do
2000     evec=ham
         call DIAG_C(Eval,Evec,hamdim,info)
         if(info/=0)then
            stop 'Error at diagonalization = SR'
         end if
         call EVEC_ASSIGN_TO_N_GEN_C(eval,evec,hamdim,N_arr,size(N_arr,dim=1))
         car=car+1

         allocate(EVECS_MAT_SR(i+1)%Evecs(hamdim))
         allocate(EVECS_MAT_SR(i+1)%Evals(hamdim))
         allocate(EVECS_MAT_SR(i+1)%N(size(N_arr,dim=1)))
         EVECS_MAT_SR(i+1)%N=N_arr
         EVECS_MAT_SR(i+1)%J=J
         EVECS_MAT_SR(i+1)%Evals=Eval
         do k = 1,hamdim
            allocate(EVECS_MAT_SR(i+1)%Evecs(k)%Ns(size(N_arr,dim=1)))
            allocate(EVECS_MAT_SR(i+1)%Evecs(k)%N_arr(size(N_arr,dim=1)))        
            do ii=1,size(N_arr,dim=1)
               N=N_arr(ii)
               allocate(EVECS_MAT_SR(i+1)%Evecs(k)%Ns(ii)%arr(2*N+1))
               EVECS_MAT_SR(i+1)%Evecs(k)%N_arr(ii)=N
               if(N>N_arr(1))then
                  idx=SumNs(N_arr(1),N-1)+1
               else
                  idx=1
               end if
               idx2=SumNs(N_arr(1),N)
               EVECS_MAT_SR(i+1)%Evecs(k)%Ns(ii)%arr=evec(idx:idx2,k)
            end do
         end do

         deallocate(ham,eval,evec,K2_arr,K1_arr,N2_arr,N1_arr)
      end do
   end subroutine CALCULATE_HAM_SR
   
   pure function SumNs(N2,N1)result(res)
      Integer,intent(in) :: n2,n1
      integer res
      INTEGER i
      
      res=0
      do i = min(N1,N2),max(N1,N2)
         res=res+2*i+1
      end do
   end function SumNs
      
   pure subroutine EVEC_ASSIGN_TO_N_GEN(eval,evec,n,n_arr,m)
      Integer,intent(in) :: n,m,n_arr(m)
      double precision,intent(inout) :: eval(n),evec(n,n)
      INTEGER i,ii,nend,nstart
      Integer,allocatable :: order(:)
      
      
      call EVEC_SORT_PICKETT(n,evec,eval)
      
      do i=1,m
         allocate(order(2*n_arr(i)+1))
         nstart=1
         do ii=1,i-1
            nstart=nstart+(N_arr(ii)*2+1)
         end do
         nend=nstart+(n_arr(i)*2)
         order=1
         call Insertion_sort_part(eval,n,nstart,nend,order)
         call Reorder_Evec(evec,n,nstart,nend,order)
         deallocate(order)
      end do
   end subroutine EVEC_ASSIGN_TO_N_GEN
   
   pure subroutine EVEC_ASSIGN_TO_N_GEN_C(eval,evec,n,n_arr,m)
      Integer,intent(in) :: n,m,n_arr(m)
      double complex,intent(inout) :: eval(n),evec(n,n)
      INTEGER i,ii,nend,nstart
      Integer,allocatable :: order(:)
      
      
      call EVEC_SORT_PICKETT_C(n,evec,eval)
      
      do i=1,m
         allocate(order(2*n_arr(i)+1))
         nstart=1
         do ii=1,i-1
            nstart=nstart+(N_arr(ii)*2+1)
         end do
         nend=nstart+(n_arr(i)*2)
         order=1
         call Insertion_sort_part_C(eval,n,nstart,nend,order)
         call Reorder_Evec_C(evec,n,nstart,nend,order)
         deallocate(order)
      end do
   end subroutine EVEC_ASSIGN_TO_N_GEN_C
   
   !eval - eigenvalue
   !evec - eigenvector
   pure subroutine EVEC_ASSIGN_TO_N(eval,evec,n,n2,n1)
      Integer,intent(in) :: n,n2,n1
      double precision,intent(inout) :: eval(n),evec(n,n)
      INTEGER i,n2end,n2start,curCol,idx,order2(2*n2+1),order1(2*n1+1)
      
      
      n2start=1 !lower
      n2end=2*n2+1
      curCol=1
      order2=1
      
      call EVEC_SORT_PICKETT(n,evec,eval)
      
      ! do i = 1,n
         ! idx=maxloc(dabs(evec(:,i)),dim=1)
         ! if(BetweenInc(n2start,idx,n2end))then
            ! call SwapCol(evec,n,i,curCol)
            ! call SwapEl(eval,n,i,curCol)
            ! curCol=curCol+1
         ! end if
      ! end do
      
      call Insertion_sort_part(eval,n,n2start,n2end,order2)
      call Reorder_Evec(evec,n,n2start,n2end,order2)
      call Insertion_sort_part(eval,n,n2end+1,n,order1)
      call Reorder_Evec(evec,n,n2end+1,n,order1)
   end subroutine EVEC_ASSIGN_TO_N
   
   pure subroutine Reorder_Evec(evec,n,startidx,endidx,order)
      Integer,intent(in) :: n,order(endIdx-startidx+1),startidx,endidx
      double precision,intent(inout) :: evec(n,n)
      INTEGER i,j
      double precision evecbuf(n,n)
      
      if(startidx==endidx)return
      
      j=1
      evecbuf=evec
      do i = startidx,endidx
         evec(:,i)=evecbuf(:,order(j)+startidx-1)
         j=j+1
      end do
   end subroutine Reorder_Evec
   
   pure subroutine Reorder_Evec_C(evec,n,startidx,endidx,order)
      Integer,intent(in) :: n,order(endIdx-startidx+1),startidx,endidx
      double complex,intent(inout) :: evec(n,n)
      INTEGER i,j
      double complex evecbuf(n,n)
      
      if(startidx==endidx)return
      
      j=1
      evecbuf=evec
      do i = startidx,endidx
         evec(:,i)=evecbuf(:,order(j)+startidx-1)
         j=j+1
      end do
   end subroutine Reorder_Evec_C
   
   pure subroutine EVEC_SORT_PICKETT_NOEVAL(n,evec)
      Integer,intent(in) :: n
      double precision,intent(inout) :: evec(n,n)
      INTEGER i,idx(2)
      double precision absevec(n,n)
      logical mask2d(n,n)
      
      mask2d=.true.
      do i = 1,n
         absevec=dabs(evec)
         idx=maxloc(absevec,mask2d)
         if(idx(1)==idx(2))goto 2222
         call ARRAY_SWAP_COLUMNS(n,evec,idx(1),idx(2),mask2d(:,idx(2)))
2222     mask2d(idx(1),:)=.false.
         mask2d(:,idx(1))=.false.
      end do 
   end subroutine EVEC_SORT_PICKETT_NOEVAL
   
   pure subroutine EVEC_SORT_PICKETT(n,evec,eval)
      Integer,intent(in) :: n
      double precision,intent(inout) :: evec(n,n),eval(n)
      INTEGER i,idx(2)
      double precision absevec(n,n)
      logical mask2d(n,n)
      
      mask2d=.true.
      do i = 1,n
         absevec=dabs(evec)
         idx=maxloc(absevec,mask2d)
         if(idx(1)==idx(2))goto 2222
         call ARRAY_SWAP_COLUMNS(n,evec,idx(1),idx(2),mask2d(:,idx(2)))
         call ARRAY1D_SWAP_COLUMNS(n,eval,idx(1),idx(2))
2222     mask2d(idx(1),:)=.false.
         mask2d(:,idx(1))=.false.
      end do 
   end subroutine EVEC_SORT_PICKETT
   
   pure subroutine EVEC_SORT_PICKETT_C(n,evec,eval)
      Integer,intent(in) :: n
      double complex,intent(inout) :: evec(n,n),eval(n)
      INTEGER i,idx(2)
      double precision absevec(n,n)
      logical mask2d(n,n)
      
      mask2d=.true.
      do i = 1,n
         absevec=abs(evec)
         idx=maxloc(absevec,mask2d)
         if(idx(1)==idx(2))goto 2222
         call ARRAY_SWAP_COLUMNS_C(n,evec,idx(1),idx(2),mask2d(:,idx(2)))
         call ARRAY1D_SWAP_COLUMNS_C(n,eval,idx(1),idx(2))
2222     mask2d(idx(1),:)=.false.
         mask2d(:,idx(1))=.false.
      end do 
   end subroutine EVEC_SORT_PICKETT_C
   
   pure subroutine ARRAY_SWAP_COLUMNS_C(n,arr,i,j,mask)
      Integer,intent(in) :: i,j,n
      logical, intent(in) :: mask(n)
      double complex,intent(inout) :: arr(n,n)
      double complex buf(n)
      INTEGER k
      
      buf=arr(:,j)
      do k = 1,n
         !if(mask(k))then
            arr(k,j)=arr(k,i)
            arr(k,i)=buf(k)
         !end if
      end do
      
      ! buf=arr(:,i)
      ! arr(:,i)=arr(:,j)
      ! arr(:,j)=buf
   end subroutine ARRAY_SWAP_COLUMNS_C
   
   pure subroutine ARRAY_SWAP_COLUMNS(n,arr,i,j,mask)
      Integer,intent(in) :: i,j,n
      logical, intent(in) :: mask(n)
      double precision,intent(inout) :: arr(n,n)
      double precision buf(n)
      INTEGER k
      
      buf=arr(:,j)
      do k = 1,n
         !if(mask(k))then
            arr(k,j)=arr(k,i)
            arr(k,i)=buf(k)
         !end if
      end do
      
      ! buf=arr(:,i)
      ! arr(:,i)=arr(:,j)
      ! arr(:,j)=buf
   end subroutine ARRAY_SWAP_COLUMNS
   
   pure subroutine ARRAY1D_SWAP_COLUMNS(n,arr,i,j)
      Integer,intent(in) :: i,j,n
      double precision,intent(inout) :: arr(n)
      double precision buf
      
      buf=arr(i)
      arr(i)=arr(j)
      arr(j)=buf
   end subroutine ARRAY1D_SWAP_COLUMNS
   
   pure subroutine ARRAY1D_SWAP_COLUMNS_C(n,arr,i,j)
      Integer,intent(in) :: i,j,n
      double complex,intent(inout) :: arr(n)
      double complex buf
      
      buf=arr(i)
      arr(i)=arr(j)
      arr(j)=buf
   end subroutine ARRAY1D_SWAP_COLUMNS_C
      
   !like factorial except addition
   function summ(num) result (res)
      INTEGER num,res,i
      res = 0
      if(num<=0)return
      
      do i = 1,num
         res = res + i
      end do
   end function summ
   
   !number of possible transitions for a pair of Js
   function Num_transitions(J1,J2) result (res)
      INTEGER res,J1,J2,cou1,cou2
      
      if(j1==j2)then
         res = summ(2*j1)
         return
      end if
      
      cou1=2*J1+1
      cou2=2*J2+1
      if(J1==0 .and. J2==0)then
         res = 0
         return
      end if
      res=cou2*cou1
   end function Num_transitions
   
   !calculate RCD from 0 to Jend
   subroutine CALCULATE_RCD_ALL()
      INTEGER J1,T1,T2,Jamount,i,m,n,idx,Js(dime_full),Km1_1,Km1_2,K1_1,K1_2,denominator,curJ,bufi,newJ
      integer,allocatable :: order(:)
      double precision buf(10),time1,time2,evals(dime_full),j2_r,j1_r,time3
      double precision, allocatable :: freqs(:)
      character(80) chars,FN
      type(TRANS_NTSJ) tr
      type(TRANS_JT),allocatable :: TRANS_JT_col_help(:)
      type(TRANS_NTSJ),allocatable :: TRANS_NTSJ_col_help(:)
      
      
      RCDCalcGoing=.true.
      if(worktype == 'SR')then !Spin rotation calculation, the basis is |N K S J M >
         
         cou = 0
         car = 1
         do J1=startJ,endJ-1
               cou = cou + Num_transitions(J1,J1+1) + Num_transitions(J1,J1)
         end do
         cou=cou+summ(2*endj)
         
         !ie. multiplicity
         Jamount=NINT(2*S+1)
         allocate(TRANS_NTSJ_col(cou*Jamount**2*20)) !¯\_(ツ)_/¯
         !allocate(TRANS_NTSJ_col(80000)) !¯\_(ツ)_/¯
         if(PreviousCalcFlag)then
            open(unitt,file=outputfile//'.trn')
7891        tr=TRANS_NTSJ_col(car)
            read(unitt,6003)tr%freq,tr%dipstr,tr%rotstr,buf(1),tr%N1,Km1_1,K1_1,tr%n2,Km1_2,K1_2,tr%en2,buf(2)
            tr%t1=Km1_1-K1_1
            tr%t2=Km1_2-K1_2
            tr%en1=tr%freq+tr%en2
            tr%dipstr=tr%dipstr*SI2AU_debye**2
            tr%rotstr=tr%rotstr*SI2AU_debye**2
            car=car+1
            goto 7891
7892        close(unitt)
6003        format(E15.8E2,1X,3(E15.8E2,1X),2(3I4,1A6),E15.8E2,1X,F8.4)
         else
            if(file_transitions(1).and.file_transitions(7))then
               open(unitt,file=outputfile//'.trn')
               write(unitt,6000)'Freq','D','R','4R/D','N1','T1','S1','J1','N2','T2','S2','J2'
               write(unitt,6001)'MHz','Debye^2','Debye^2'
               write(unitt,5999)
5999           format(100(1h-))
6000           format(4(A15,1X),2(3A4,1A6))
6001           format(3(A15,1X))
            end if
            j2_r=S
            write(6,*)'  J'
            car =0
            do while(j2_r<=endj)
               call CPU_TIME(time2)
               write(6,601,advance='no')R2HI(j2_r)
               call CALCULATE_RCD_Jpair_SR(j2_r,j2_r)
               call CALCULATE_RCD_Jpair_SR(j2_r,j2_r+1)
               j2_r=j2_r+1
               call CPU_TIME(time1)
               write(6,600)time1-time2
600            format(3X,'...',1F8.3,' seconds')
601            format(1A5)
            end do
            

            if(file_transitions(1).and.file_transitions(7))then
               close(unitt)
            end if
         end if
         if(dipstr_cutoff_percentual)then
            write(6,*)'Will use percentual dipole strength cutoff: ',dipstr_cutoff*100d0,'%'
            dipstr_cutoff=dipstr_cutoff*maxval(TRANS_NTSJ_col(:)%dipstr)*AU2SI_debye**2
         end if
         if(orderBy=='FREQ')then
            allocate(order(car),freqs(car))
            do i = 1,car
               !order(i)=i
               freqs(i)=TRANS_NTSJ_col(i)%freq
            end do
            call Insertion_sort(freqs,car,order)
            TRANS_NTSJ_col_help=TRANS_NTSJ_col
            do i = 1,car
               TRANS_NTSJ_col(i)=TRANS_NTSJ_col_help(order(i))
            end do
            deallocate(TRANS_NTSJ_col_help)
         end if
      else if(worktype == 'RR')then
      
615   format(A127)
      car=1
      cou=0
      idx=0
      do i = startJ,endJ
         do m = 1,2*i+1
            idx=idx+1
            Js(idx)=i
         end do
      end do
      allocate(Fracpop(dime_full))
      n=1
      do i = startJ,endJ
         do m = 1,2*i+1
            evals(n)=eval_full_list(i-startJ+1)%arr(m)
            n=n+1
         end do
      end do
      Fracpop=Calc_Frac_Populations(evals*1d9,Js,dime_full,Q,temp)
      
      do J1=startJ,endJ-1
            cou = cou + Num_transitions(J1,J1+1) + Num_transitions(J1,J1)
      end do
      cou=cou+summ(2*endj)
      allocate(TRANS_JT_col(cou))
      if(separateElNuc)then
         do i = 1,cou
            allocate(TRANS_JT_col(i)%rots_contr(4),TRANS_JT_col(i)%dips_contr(4))
         end do
      end if
      
      write(6,*)
      write(6,2018)'J1  - J2'
      denominator=NINT(dble(endJ-startJ)/min(5d0,dble(endj)))
      call CPU_TIME(time3)
      call CPU_TIME(time1)
      curJ=startJ
      call CALCULATE_RCD_Jpair(0,1)
      do i = startJ,endJ
         if(mod(i,denominator)==0 .and. i/=0)then
            J1=i
            exit
         end if
      end do
      write(6,2015,advance='no')curJ,' - ',min(J1,endj),' ...'
      !the paralellization would probably look like this:
      !make a list of transition quantum numbers
      !hand them out to CPUs one by one as they go on
      do J1=1,endJ-1
         call CALCULATE_RCD_JPAiR(J1,J1)
         call CALCULATE_RCD_Jpair(J1,J1+1)
         if(mod(J1,denominator)==0 .and. J1/=0)then
            write(6,2016)(GET_CPU_TIME()-time1),' seconds'
            do i = J1+1,endJ+1
               if(mod(i,denominator)==0)then
                  newJ=i
                  exit
               end if
            end do
            write(6,2015,advance='no')J1,' - ',newJ,' ...'
            
            call CPU_TIME(time1)
         end if
      end do
      call CALCULATE_RCD_JPAiR(endJ,endJ)
      call CPU_TIME(time2)
      write(6,2016)(time2-time3),' seconds'
      close(30)
2015  format(I3,A4,I3,A4)     
2016  format(1F7.0,A8) 
2014  format(A1,I3,A1)
2017  format(A40)      
2018  format(A8)
      close(667)
         if(dipstr_cutoff_percentual.and.(file_spectrum(2).or.file_transitions(6)))then
            write(6,*)'Will use percentual dipole strength cutoff: ',dipstr_cutoff*100d0,'%'
            dipstr_cutoff=dipstr_cutoff * maxval(trans_jt_col(:)%dipstr)*AU2SI_debye**2
         end if
         if(orderBy=='FREQ')then
            allocate(order(car),freqs(car))
            do i = 1,car
               !order(i)=i
               freqs(i)=trans_jt_col(i)%freq
            end do
            call Insertion_sort(freqs,car,order)
            TRANS_JT_col_help=trans_jt_col
            do i = 1,car
               trans_jt_col(i)=trans_jt_col_help(order(i))
            end do
            deallocate(TRANS_JT_col_help)
         end if
      end if
      RCDCalcGoing=.false.
      rotstr_cutoff=rotstr_cutoff*SI2AU_debye**2
      dipstr_cutoff=dipstr_cutoff*SI2AU_debye**2
      if((file_spectrum(2).or.file_transitions(6)))then
         write(6,'(A35,E11.5E2,A8)')'The threshold for D_ij is ',dipstr_cutoff*au2si_debye**2,' debye^2'
         write(6,'(A35,E11.5E2,A8)')'The threshold for R_ij is ',rotstr_cutoff*au2si_debye**2,' debye^2'
         write(6,'(A35,E11.5E2)')'The threshold for 4R_ij/D_ij is ',kuhnpar_cutoff
      end if
      call WriteOutput(0)
   end subroutine CALCULATE_RCD_ALL
      
   function GET_CPU_TIME()result(res)
      double precision res
      call CPU_TIME(res)
   end function GET_CPU_TIME
      
   !I DONT EVEN WANNA EXPLAIN THIS ONE
   !BECAUSE SPCAT OUTPUT USES LETTERS INSTEAD OF 100,110,120, etc...
   function CheckJPLLetters(FN) result(res)
      character(24) FN
      character(25) checkStr
      INTEGER i
      logical res
      
      checkStr='ABCDEFGHIJKLMNOPQRSTUVWYZ'
      
      res = .true.
      do i=1,25
         if(index(FN,checkStr(i:i))/=0)then
            res = .false.
            return
         end if
      end do
      
      
   end function CheckJPLLetters
            
   function IARR_areEqual(arr1,arr2,n)result(res)
      INTEGER n,i
      INTEGER arr1(n),arr2(n)
      logical res
      
      res=.true.
      do i = 1,n
         if(arr1(i)/=arr2(i))then
            res = .false.
            return
         end if
      end do
      continue
   end function IARR_areEqual
   
   function Number2JPL(num)result(res)
      integer num
      character(1),parameter :: chars100(26) = ['A','B','C','D','E','F','G','H','I','J',&
      'K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
      character(2) res
      character(1) help
      
      
      if(num>359 .or. num<-259)then
         res='**'
      else if(num>100)then
         write(res,'(A1,I1)')chars100((num-100)/10+1),mod(num,10)
      else if(num<=-10)then
         help = chars100((num+10)/10+1)
         call to_Lower(help)
         write(res,'(A1,I1)')help,mod(-num,10)
      else 
         write(res,'(I2)')num
      end if
   end function Number2JPL
   
   subroutine WriteOutput_SR(num)
      INTEGER i,num,k2m,k2p,k1m,k1p,bufI(3),deg
      type(TRANS_NTSJ) tr
      double precision RD,rs,ds,dipstrmax,popfac,jplinten,fp1,fp2
      logical fileOpened
      character(2) catGround(4),catExcited(4)
      character(1) catSign
      
      
      if(num>0)then
         inquire(unitt,opened=fileOpened)
         if(.not.fileOpened)then
            open(unitt,file=outputfile//'.trn')
            write(unitt,6000)'Freq','D','R','4R/D','N1','T1','J1','N2','T2','J2'
            write(unitt,5999)
         end if
         tr=TRANS_NTSJ_col(num)
         ds=tr%dipstr*AU2SI_debye**2
         rs=tr%rotstr*AU2SI_debye**2
         bufI=OldToNewNotation(tr%n1,tr%t1)
         k1m=bufI(2)
         k1p=bufI(3)
         bufI=OldToNewNotation(tr%n2,tr%t2)
         k2m=bufI(2)
         k2p=bufI(3)
         if(.not.file_transitions(8).and.(abs(rs)<outtol.or. ds<outtol))return
         if(ds/=0d0)then
            RD=4d0*rs/ds
         end if
         
         write(unitt,6001)tr%freq,ds,rs,RD,tr%N1,k1m,k1p,R2HI(tr%J1),tr%N2,k2m,k2p,R2HI(tr%J2)
         return
      end if
      if(file_transitions(1))then
         if(file_transitions(9))then
            open(unitt,file=outputfile//'_S'//spinMultChar//'.trn')
         else
            open(unitt,file=outputfile//'.trn')
         end if
      end if
      
      if(file_cat(1))then
         if(MWFlag)then
            open(unitt+1,file=outputfile//'.cat')
         end if
         open(unitt+2,file=outputfile//'_RCD'//'.cat')
      end if
         
      if(file_transitions(5))then
         write(unitt,6000)'Freq','D','R','4R/D','N1','Km1','Kp1','J1','N2','Km2','Kp2','J2'
         write(unitt,6002)'GHz','Debye^2','Debye^2'
         write(unitt,5999)
      end if
      dipstrMax = maxval(TRANS_NTSJ_col(:)%dipstr)*AU2SI_debye**2
      do i = 1,car-1
         tr=TRANS_NTSJ_col(i)
         ds=tr%dipstr*AU2SI_debye**2
         rs=tr%rotstr*AU2SI_debye**2
         if(ds/=0d0)then
            RD=4d0*rs/ds
         end if
         fp2=exp(-dble(tr%en2)*1d9*h*beta)
         fp1=exp(-dble(tr%en1)*1d9*h*beta)
         popfac=(2*tr%N2+1)*(fp2-fp1)/Q
         if(file_transitions(6))then
            if(ds<dipstr_cutoff)cycle
            if(dabs(rs) < rotstr_cutoff .or. dabs(RD) < kuhnpar_cutoff)cycle
         end if
         if(.not.file_transitions(8).and.(abs(rs)<outtol.or. abs(ds)<outtol))cycle
         !if(abs(ds)<1d-4)cycle !!TODO REMOVE
         jplinten=4.16231D-5*(tr%freq)*1D3*ds*(3d0/2d0)*popfac
         if(jplinten==0d0)then
            jplinten=1000
         else
            jplInten=log10(abs(jplinten))
         end if
         
         bufI=OldToNewNotation(tr%n1,tr%t1)
         k1m=bufI(2)
         k1p=bufI(3)
         bufI=OldToNewNotation(tr%n2,tr%t2)
         k2m=bufI(2)
         k2p=bufI(3)
         if(file_transitions(2))then
            write(unitt,6003)tr%freq,ds,rs,RD,tr%N1,k1m,k1p,tr%J1,tr%N2,k2m,k2p,tr%J2,jplinten
         else
            write(unitt,6004)tr%freq,ds,rs,RD,tr%N1,k1m,k1p,tr%J1,tr%N2,k2m,k2p,tr%J2
         end if
         
         !write to .cat format file
         if(file_cat(1))then
            catGround=[Number2JPL(tr%N2),Number2JPL(k2m),Number2JPL(k2p),Number2JPL(NINT(tr%J2*2))]
            catExcited=[Number2JPL(tr%N1),Number2JPL(k1m),Number2JPL(k1p),Number2JPL(NINT(tr%J1*2))]
            deg=NINT(2*tr%j1+1)
            if(MWFlag)then
               if(jplInten>-99.0 .and. abs(jplInten) >= 0.0001)then
                  write(unitt+1,5468)tr%freq*1d3,0d0,jplInten,3,dble(tr%en2*GHZ_2_cm),deg,catTag,314,catExcited,catGround
               end if
            end if
            if(rs<0d0)then
               catSign='-'
            else
               catSign='+'
            end if
            jplInten=4.16231D-5*dabs(tr%freq)*1D3*rs*(3d0/2d0)*(2*tr%N2+1)*dabs(fp1-fp2)/Q!This is the JPL intensity used in spcat
            jplInten=log10(dabs(jplInten))
            write(unitt+2,5469)tr%freq*1d3,0d0,jplInten,3,dble(tr%en2*GHZ_2_cm),deg,catTag,314,catExcited,catGround,catSign
         end if
5468     format(F13.4, F8.4, F8.4, I2, F10.4, I3, A7, I4, 4A2,6X, 4A2,6X) 
5469     format(F13.4, F8.4, F8.4, I2, F10.4, I3, A7, I4, 4A2,6X, 4A2,6X,A1)      
      end do
      
      close(unitt)
      
      !write(6,*)'Written '//outputfile//'.trn'
      if(file_transitions(1))then
         if(file_transitions(9))then
            write(6,*)'Written '//outputfile//'_S'//spinMultChar//'.trn'
         else
            write(6,*)'Written '//outputfile//'.trn'
         end if
      end if
      
      if(file_cat(1))then
         if(MWFlag)then
            write(6,*)'Written '//outputfile//'.cat'
         end if
         write(6,*)'Written '//outputfile//'_RCD'//'.cat'
      end if
      
6000  format(4(A15,1X),2(3A4,1A6))
6002  format(3(A15,1X))
6001  format(E15.8E2,1X,3(E14.8E2,1X),2(3I4,1A6))
6003  format(E15.8E2,1X,3(E14.8E2,1X),2(3I4,F6.1),1X,E15.8E2,1X,F8.4)
6004  format(E15.8E2,1X,3(E14.8E2,1X),2(3I4,F6.1),1X,E15.8E2)
5999  format(180(1h-))
   end subroutine WriteOutput_SR
      
   !writes .trn file
   subroutine WriteOutput(num)
      INTEGER num
      double precision dips,rots,ener,rbuf,buf,transProb,Rot_strs_(cou),Dip_strs_(cou),trans_en_arr_(cou),fp1,fp2,jplInten,en1,en2,dipstrMax,abscoeff,abscoeffMW
      double precision,allocatable :: eval1(:),eval2(:)
      INTEGER i,j1,j2,t1,t2,idx1,idx2,bufI(3),Km1_1,Km1_2,K1_1,K1_2,j,bufi2
      character(:),allocatable :: numform1,numform2,numform3,numform4
      character(7) energyUn,strengthsUn
      character(1) catSign
      character(2) catGround(3),catExcited(3)
      character(14) rotStrengthUn
      type(trans_jt) tr
      
      energyUn='GHz'
      strengthsUn='Debye^2'
      rotStrengthUn='Debye^2'
      if(.not.file_transitions(1) .or. file_transitions(7).and.worktype/='RR')return
      if(file_transitions(8))then
         write(6,*)'Showing D_ij and R_ij of zero value.'
      end if
      if(worktype == 'SR')then
         call WriteOutput_SR(num)
      else if(worktype == 'RR')then
         
         if(file_transitions(1))then
            if(file_transitions(9))then
               open(unitt,file=outputfile//'_S'//spinMultChar//'.trn')
            else
               open(unitt,file=outputfile//'.trn')
            end if
         end if
         
         if(file_cat(1))then
            if(MWFlag)then
               open(unitt+1,file=outputfile//'.cat')
            end if
            open(unitt+2,file=outputfile//'_RCD'//'.cat')
         end if
         
         ! if(index(JPLFile,'.cat')/=0)then
            ! open(unitt+2,file=outputfile//'.jpl')
            ! write(unitt+2,658)'Frequency (GHz)','JPL Frequency','Int (nm^2 MHz)','JPL Intensity','J1','K-1','K1','J2','K-1','K1'
         ! end if
         if(.not.PreviousCalcFlag)then
   658   format(4A17,6A3)
         if(file_transitions(5))then
            write(unitt,6621)energyUn,strengthsUn,strengthsUn
6621         format('           freq              D              R           4R/D  J1 K-1  K1  J2 K-1  K1       JPL_INTEN',/,'           ',A7,'    ',A7,' ',A14)
         end if
   
         dipstrMax = maxval(trans_jt_col(:)%dipstr)
         do i = 1,cou
         tr=trans_jt_col(i)
         dips=tr%dipstr*au2si_debye**2
         rots=tr%rotstr*au2si_debye**2
         ener=tr%freq
         if(.not.file_transitions(8).and.(dips<outTol.or.dabs(rots)<outTol))cycle
         if(.not.(abs(dips)<=0.0d0))then
            buf = 4.0d0*rots/dips
         else
            buf = 0.0d0
         end if
                  
         if(file_transitions(6))then
            if(dips<dipstr_cutoff)cycle
            if((dabs(rots) < rotstr_cutoff .or. dabs(buf) < kuhnpar_cutoff))cycle
         end if
         j1=tr%j1
         t1=tr%t1
         j2=tr%j2
         t2=tr%t2
         
         bufI=OldToNewNotation(j1,t1)
         Km1_1=bufI(2)
         K1_1=bufI(3)
         bufI=OldToNewNotation(j2,t2)
         Km1_2=bufI(2)
         K1_2=bufI(3)
         
         idx1=Calc_Frac_index(j1,t1)
         idx2=Calc_Frac_index(j2,t2)
         
         allocate(eval1(2*j1+1),eval2(2*j2+1))
         call SLICE_EVAL(eval1,J1)
         call SLICE_EVAL(eval2,J2)
         
         en1=eval1(T1+J1+1)
         en2=eval2(T2+J2+1)
         
         fp1=exp(-en1*1D9*h/kb/temp)
         fp2=exp(-en2*1D9*h/kb/temp)
         
         jplInten=4.16231D-5*dabs(ener)*1D3*dips*(3d0/2d0)*dble(2*j1+1)*dabs(fp1-fp2)/Q!This is the JPL intensity used in spcat
         !To get the same value as SPCAT, I recommend setting temperature to 300 (K) and partition function to 1 
         if(jplInten < 1d-99)then!small intensity, dont bother with it
            jplInten = 0d0
         else
            jplInten = log10(jplInten)
         end if

         if(file_transitions(4))then
            numform1="(F15.8,3(' ',F14.8),6i4,' ')"
            numform2="(F8.4)"
            numform3="(2(1X,F15.8))"
            numform4="(4(1X,F15.8))"
         else
            numform1="(E15.8E2,3(' ',E14.8E2),6i4,' ')"
            numform2="(F8.4)"
            numform3="(2(1X,E15.8E2))"
            numform4="(4(1X,E15.8E2))"
         end if
         
         write(unitt,numform1,advance='no')ener,dips,rots,buf,j2,Km1_2,K1_2,j1,Km1_1,K1_1
         
         if(file_transitions(2))then
            write(unitt,numform2,advance='no')jplInten
         end if
         
         !write to .cat format file
         if(file_cat(1))then
            catGround=[Number2JPL(j1),Number2JPL(Km1_1),Number2JPL(K1_1)]
            catExcited=[Number2JPL(j2),Number2JPL(Km1_2),Number2JPL(K1_2)]
            if(MWFlag)then
               if(jplInten>-99.0 .and. abs(jplInten) >= 0.0001)then
                  write(unitt+1,5468)ener*1d3,0d0,jplInten,3,en1*GHZ_2_cm,2*j2+1,catTag,303,catExcited,catGround
               end if
            end if
            if(rots<0d0)then
               catSign='-'
            else
               catSign='+'
            end if
            jplInten=4.16231D-5*dabs(ener)*1D3*dabs(rots)*(3d0/2d0)*dble(2*j1+1)*dabs(fp1-fp2)/Q!This is the JPL intensity used in spcat
            if(jplInten>1d-99)then
               jplInten=log10(jplInten)
               write(unitt+2,5469)ener*1d3,0d0,jplInten,3,en1*GHZ_2_cm,2*j2+1,catTag,303,catExcited,catGround,catSign
            end if
         end if
         
5468     format(F13.4, F8.4, F8.4, I2, F10.4, I3, A7, I4, 3A2,6X, 3A2,6X) 
5469     format(F13.4, F8.4, F8.4, I2, F10.4, I3, A7, I4, 3A2,6X, 3A2,6X,A1)      
         
         if(separateElNuc)then
            write(unitt,numform4,advance='no')tr%dips_contr(1),tr%dips_contr(2),tr%dips_contr(3),tr%dips_contr(4)
            write(unitt,numform4,advance='no')tr%rots_contr(1),tr%rots_contr(2),tr%rots_contr(3),tr%rots_contr(4)
         end if
         
         write(unitt,*)
            
   664   format(F15.8,3(' ',F14.8),6i4,' ')
   665   format(F8.4)
   666   format(2F15.8)

   661   format(E15.8E2,3(' ',E14.8E2),6i4,' ')
   662  format(F8.4)
   663  format(2E15.8E2)
         ! if(index(JPLFile,'.cat')/=0)then
            ! do j=1,size(trans_en_arr_JPL)
               ! if(IARR_areEqual([j1,km1_1,k1_1,j2,km1_2,k1_2],trans_arr_JPL(:,j),6))then
                  ! write(unitt+2,659)ener,trans_en_arr_JPL(j),jplInten,inten_JPL(j),j1,km1_1,k1_1,j2,km1_2,k1_2
               ! else
               
               ! end if
            ! end do
         ! end if
         
         deallocate(eval1,eval2)
         end do
   659   format(2(1X,F16.8),2(1X,E16.10),6I3)      
         end if
         if(file_transitions(1))then
            close(unitt)
            write(6,*)'WRITTEN '//outputfile//'.trn'
         end if
         
         if(file_cat(1))then
            if(MWFlag)then
               close(unitt+1)
               write(6,*)'WRITTEN '//outputfile//'.cat'
            end if
            close(unitt+2)
            write(6,*)'WRITTEN '//outputfile//'_RCD.cat'
         end if
         ! if(index(JPLFile,'.cat')/=0)then
            ! close(unitt+2)
         ! end if
      end if
      continue
   end subroutine WriteOutput
     
   subroutine CALCULATE_RCD_Jpair_SR(J2,J1)
      INTEGER i,j,N1,N2,T1,T2,fj1,fj2,idx2,idx1,p,idxhelp,nthr,idx,thr_idx,start_idx,it,jt,carr,car_copy
      Integer,allocatable :: idx_col(:)
      INTEGER :: Km1,Kp1,Km2,Kp2,bufarr(3),F,m2_idx,m1_idx,n2n,n1n
      double precision J1,J2,dipstr,rotstr,ener,m2,m1,help,cgs(3),j2n,j1n,cg_p,dipstr_new,rotstr_new
      type(EVECS_SPINROT) J2_col,J1_col,J_help !collection
      type(TUPLE_EVEC_SPINROT) c2,c1,chelp
      type(TRANS_NTSJ),allocatable :: trans_arr(:,:)
      double complex :: dipsc,rotsc,dipsc_new,rotsc_new,dircoseM,small_td
      double complex,allocatable :: trans_dips1(:)
      
      if(IsHalfInt(S))then
         fj1=floor(j1)+1
         fj2=floor(j2)+1
      else
         fj1=NINT(j1)
         fj2=NINT(j2)
      end if
      nthr=-1000
      !allocate(evec2(dims2,dims2),evec1(dims1,dims1),eval2(dims2),eval1(dims1))
      5698 format(A10,1X,I2)
      5699 format(A14,1X,I2)
      J2_col=EVECS_MAT_SR(fj2)
      J1_col=EVECS_MAT_SR(fj1)
      carr=0
      dircoseM=0d0
      do F = 1,2
         do M2_idx = 0,NINT(2*J2)
            m2=m2_idx-J2
            do m1_idx = 0,NINT(2*J1)
               m1=m1_idx-J1
               dircoseM=dircoseM + DIR_COSE_JM_FROM_SPH(J2,M2,J1,M1,axes(F))*DIR_COSE_JM_FROM_SPH(J1,M1,J2,M2,axes(F))
            end do
         end do
      end do
      do i = 1,size(J2_col%N,dim=1)
         N2=J2_col%N(i)
         do j = 1,size(J1_col%N,dim=1)
            N1=J1_col%N(j)
            if(N2>N1)cycle
            do T2=-N2,N2
               do T1=-N1,N1
                  if(N1==N2 .and. J1==J2 .and. t1==t2 .or. n1==n2 .and. t2>t1)cycle
                  idx1=N1+T1+1
                  if(j/=1)then
                     idx1=SumNs(J1_col%N(1),N1-1)+idx1
                  end if
                  idx2=N2+T2+1
                  if(i/=1)then
                     idx2=SumNs(J2_col%N(1),N2-1)+idx2
                  end if
                  
                  ener=dble(J1_col%Evals(idx1)-J2_col%Evals(idx2))
                  if(ener==0d0 .or. SkipTransFlag.and.(abs(ener)<fmin .or. abs(ener)>fmax))cycle

                  !I did it like this I got confused by all the quantum numbers
                  if(ener<0d0)then
                     c2=J1_col%Evecs(idx1)
                     c1=J2_Col%Evecs(idx2)
                     J2n=J1
                     J1n=J2
                     n2n=n1
                     n1n=n2
                  else
                     c2=J2_Col%Evecs(idx2)
                     c1=J1_col%Evecs(idx1)
                     J2n=J2
                     J1n=J1
                     n2n=n2
                     n1n=n1
                  end if
                  allocate(trans_dips1(NINT(2*J2n+1)))
                  
                  !call TRANS_DIPSTR_SR_CART_FROM_SPH(c2,c1,J2n,J1n,S,NINT(2*S+1),dipsc,dip,trans_dips1)
                  !call TRANS_ROTSTR_SR_CART_FROM_SPH(c2,c1,J2n,J1n,S,NINT(2*S+1),rotsc,dip,gt,trans_dips1)
                  call TRANS_DIPSTR_SR_CART_FROM_SPH_NEW(c2,c1,J2n,J1n,S,NINT(2*S+1),dipsc_new,dip,small_td,dircoseM)
                  call TRANS_ROTSTR_SR_CART_FROM_SPH_NEW(c2,c1,J2n,J1n,S,NINT(2*S+1),rotsc_new,dip,gt,small_td,dircoseM)
                  !dipstr=4d0*dble(dipsc)/(2*J2n+1)
                  !rotstr=4d0*aimag(rotsc)/(2*J2n+1)
                  dipstr=dble(dipsc_new)/(2*n2n+1)
                  rotstr=aimag(rotsc_new)/(2*n2n+1)
                  if(ener<0d0)then
                     call NKSJ_FillTrans(TRANS_NTSJ_col(car+1),J1_col%Evals(idx1),J2_col%Evals(idx2),N1,N2,T1,T2,S,S,J1,J2,dipstr,rotstr,abs(dble(ener)),car)
                  else
                     call NKSJ_FillTrans(TRANS_NTSJ_col(car+1),J2_col%Evals(idx2),J1_col%Evals(idx1),N2,N1,T2,T1,S,S,J2,J1,dipstr,rotstr,abs(dble(ener)),car)
                  end if
                  deallocate(trans_dips1)
               end do
            end do
         end do
      end do
   end subroutine CALCULATE_RCD_Jpair_SR
      
      
   subroutine TRANS_DIPSTR_SR_CART_FROM_SPH_NEW(c2,c1,J2,J1,S,mult,dipsc,dip,trans_dips1,dircoseM)
      double precision, intent(in) :: J1,J2
      type(TUPLE_EVEC_SPINROT), intent(in) :: c2,c1
      double precision, intent(in) :: S,dip(3)
      double complex,intent(out) :: dipsc,trans_dips1
      double complex,intent(in) :: dircoseM
      integer, intent(in) :: mult
      
      integer :: N2_idx,N1_idx,N2,N1,Mn2,Mn1,Ms2_idx,M2_idx,M1_idx,K2,K1,a,F,diffJ
      double precision :: Ms2,M2,M1,cg1,cg2
      double complex,allocatable :: evec2(:),evec1(:)
      character :: labaxis
      character, parameter :: axes(3) = ['X','Y','Z']
      double complex :: bufm1,bufm2,bufk1,bufk2,bufn1,bufn2,contrs_F(3),bufa1,bufa2,bufm
      
      dipsc=0d0
      diffJ=abs(NINT(J2-J1))
      contrs_F=0d0
      labaxis=axes(3)
      
      bufn2=0d0
      bufn1=0d0
      do N2_idx = 1,mult
         N2=c2%N_arr(N2_idx)
         evec2=conjg(c2%Ns(N2_idx)%arr)
         do N1_idx = 1,mult
            N1=c1%N_arr(N1_idx)
            if(N1+N2==0)cycle
            evec1=c1%Ns(N1_idx)%arr
            bufk1=0d0
            bufk2=0d0
            do K2 = -N2,N2
               do K1 = max(-N1,K2-1),min(N1,K2+1)
                  bufa1=0d0
                  bufa2=0d0
                  do a = 1,3
                     bufa1 = bufa1+DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1,axes(a))*dip(a)
                     bufa2 = bufa2+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2,axes(a))*dip(a)
                  end do
                  bufk1 = bufk1 + bufa1*evec1(N1+K1+1)*evec2(N2+K2+1)
                  bufk2 = bufk2 + bufa2*evec1(N1+K1+1)*evec2(N2+K2+1)
               end do
            end do
            
            bufn1 = bufn1+DIR_COSE_NSJ(N2,S,J2,N1,S,J1)*bufk1
            bufn2 = bufn2+DIR_COSE_NSJ(N1,S,J1,N2,S,J2)*bufk2
         end do
      end do
      trans_dips1=bufn1
      dipsc=(bufn1*bufn2)*dircoseM
   end subroutine TRANS_DIPSTR_SR_CART_FROM_SPH_NEW
      
   subroutine TRANS_DIPSTR_SR_CART_FROM_SPH(c2,c1,J2,J1,S,mult,dipsc,dip,trans_dips1)
      double precision, intent(in) :: J1,J2
      type(TUPLE_EVEC_SPINROT), intent(in) :: c2,c1
      double precision, intent(in) :: S,dip(3)
      double complex,intent(out) :: dipsc,trans_dips1(NINT(2*J2+1))
      integer, intent(in) :: mult
      
      integer :: N2_idx,N1_idx,N2,N1,Mn2,Mn1,Ms2_idx,M2_idx,M1_idx,K2,K1,a,F,diffJ
      double precision :: Ms2,M2,M1,cg1,cg2
      double complex,allocatable :: evec2(:),evec1(:)
      character :: labaxis
      character, parameter :: axes(3) = ['X','Y','Z']
      double complex :: bufm1,bufm2,bufk1,bufk2,bufn1,bufn2,contrs_F(3),bufa1,bufa2,bufm
      
      dipsc=0d0
      diffJ=abs(NINT(J2-J1))
      contrs_F=0d0
      !do F = 1,3
      labaxis=axes(3)
      
      do M2_idx = 0,NINT(2*J2)
         M2=M2_idx-J2
      !do m1_idx = max(m2_idx-1+diffJ,0),min(m2_idx+1+diffJ,NINT(2*J1))
         M1=m2
         !M1=m2+1
         bufn2=0d0
         bufn1=0d0
         do N2_idx = 1,mult
            N2=c2%N_arr(N2_idx)
            evec2=conjg(c2%Ns(N2_idx)%arr)
            do N1_idx = 1,mult
               N1=c1%N_arr(N1_idx)
               if(N1+N2==0)cycle
               evec1=c1%Ns(N1_idx)%arr
               
               bufM1=0d0
               bufm2=0d0
               !do Mn2 = -N2,N2
                 ! Mn1=NINT(M1-M2)+Mn2
                  do Ms2_idx=0,NINT(2*S)
                     Ms2=Ms2_idx-S
                     Mn2=NINT(M2-Ms2)
                     Mn1=NINT(M1-Ms2)
                     cg1=CG_SR(N2,S,J2,Mn2,Ms2,M2)
                     cg2=CG_SR(N1,S,J1,Mn1,Ms2,M1)
                     bufm1=bufm1+cg1*cg2*DIR_COSE_NM_FROM_SPH(N2,Mn2,N1,Mn1,labaxis)
                     bufm2=bufm2+cg1*cg2*DIR_COSE_NM_FROM_SPH(N1,Mn1,N2,Mn2,labaxis)
                  end do
              ! end do
               
               bufK1=0d0
               bufK2=0d0
               do K2 = -N2,N2
                  do K1 = max(-N1,K2-1),min(N1,K2+1)
                     bufa1=0d0
                     bufa2=0d0
                     do a = 1,3
                        bufa1 = bufa1+DIR_COSE_NK_FROM_SPH(N2,K2,N1,K1,axes(a))*dip(a)
                        bufa2 = bufa2+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2,axes(a))*dip(a)
                     end do
                     bufk1 = bufk1 + bufa1*evec1(N1+K1+1)*evec2(N2+K2+1)
                     bufk2 = bufk2 + bufa2*evec1(N1+K1+1)*evec2(N2+K2+1)
                  end do
               end do
               
               bufn1 = bufn1+DIR_COSE_N_FROM_SPH(N2,N1)*bufk1*bufm1
               bufn2 = bufn2+DIR_COSE_N_FROM_SPH(N1,N2)*bufk2*bufm2
            end do
         end do
         dipsc = dipsc + bufn1*bufn2
         trans_dips1(M2_idx+1)=bufn1
      end do
      !end do
      !end do
      !dipsc = sum(contrs_F)
   end subroutine TRANS_DIPSTR_SR_CART_FROM_SPH
   
   subroutine TRANS_ROTSTR_SR_CART_FROM_SPH_NEW(c2,c1,J2,J1,S,mult,rotsc,dip,gt,trans_dips1,dircoseM)
      double precision, intent(in) :: J1,J2
      type(TUPLE_EVEC_SPINROT), intent(in) :: c2,c1
      double precision, intent(in) :: S,dip(3),gt(3,3)
      double complex,intent(out) :: rotsc
      double complex,intent(in) :: trans_dips1,dircoseM
      integer, intent(in) :: mult
      
      integer :: N2_idx,N1_idx,N2,N1,mn3,Mn2,Mn1,Ms2_idx,M2_idx,M1_idx,K2,K1,a,i,n3l,n3u,k3,n3,ks1_idx,N_idx,F,diffJ
      double precision :: ms1,Ms2,M2,M1,cg1,cg2,kj1,ks1,ks3
      double complex,allocatable :: evec2(:),evec1(:)
      character :: labaxis 
      character, parameter :: axes(3) = ['X','Y','Z']
      double complex :: bufm1,bufm2,bufk1,bufk2,bufk3,bufn1,bufn2,bufsK,bufa,bufs,bufm3,magel,magel2,magel2m,dircosjm
      double complex :: ga(3,3),gcor(3),mag1,mag2,contrs_F(3),contrs_F_magel(3),contrs_F_magnuc(3)
      
      
      do i = 1,3
         ga(i,1)=(gt(i,1)-iu*gt(i,2))*0.5d0
         ga(i,2)=gt(i,3)
         ga(i,3)=(gt(i,1)+iu*gt(i,2))*0.5d0
      end do
       
      gcor(1)=0.5d0*iu*(gt(2,3)-gt(3,2))
      gcor(2)=0.5d0*iu*(gt(3,1)-gt(1,3))
      gcor(3)=0.5d0*iu*(gt(1,2)-gt(2,1))
      rotsc=0d0
      diffJ=abs(NINT(J2-J1))
      labaxis=axes(3)
      magel2=0d0
      if(magdip_contr(2))then
         do F = 1,2
            do M2_idx = 0,NINT(2*J2)
               M2=M2_idx-J2
               do m1_idx = max(m2_idx-1+diffJ,0),min(m2_idx+1+diffJ,NINT(2*J1))
                  M1=m1_idx-J1
                  dircosjm=DIR_COSE_JM_FROM_SPH(J2,M2,J1,M1,axes(F))
                  if(abs(dircosjm)==0d0)cycle
                  magel=0d0
                  do N2_idx = 1,mult
                     N2=c2%N_arr(N2_idx)
                     evec2=c2%Ns(N2_idx)%arr
                     N_idx=findloc(c1%N_arr,N2,dim=1)
                     if(N_idx>0)then
                        evec1=conjg(c1%Ns(N_idx)%arr)
                        bufs=0d0
                        do Ms2_idx=0,NINT(2*S)
                           Ms2=Ms2_idx-S
                           mn2=M2-Ms2
                           Ms1=M1-M2+Ms2
                           cg1=CG_SR(N2,S,J2,Mn2,Ms2,M2)
                           cg2=CG_SR(N2,S,J1,Mn2,Ms1,M1)
                           bufs=bufs+cg1*cg2*SAlpha_matel_space(S,Ms1,S,Ms2,axes(F))
                        end do
                        
                        bufsk=0d0
                        do K2=-N2,N2
                           bufsk = bufsk+evec2(N2+K2+1)*evec1(N2+K2+1)
                        end do
                        magel=magel+bufsk*bufs
                     end if
                  end do
                  magel2=magel2+magel*dircosjm
               end do
            end do
         end do
         magel2=magel2*bohrm_cgs*gs
      end if
      
      
      bufn2=0d0
      if(magdip_contr(1))then
         do N2_idx = 1,mult
            N2=c2%N_arr(N2_idx)
            evec2=conjg(c2%Ns(N2_idx)%arr)
            do N1_idx = 1,mult
               N1=c1%N_arr(N1_idx)
               if(N1+N2==0)cycle
               evec1=c1%Ns(N1_idx)%arr
               bufk3=0d0
               do K2 = -N2,N2
                  do K1 = max(-N1,K2-2),min(N1,K2+2)
                     bufk2=0d0
                     !nuclear magnetic moment in molecule-fixed axes
                     do a = 1,3
                        bufK2=bufk2+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2+1,axes(a))*Jminus_eval(N2,K2)*ga(a,3)+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2-1,axes(a))*Jplus_eval(N2,K2)*ga(a,1)+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2,axes(a))*K2*ga(a,2)
                        bufk2=bufk2+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2,axes(a))*gcor(a)
                     end do
                     bufk3=bufk3+bufk2*evec1(N1+K1+1)*evec2(N2+K2+1)
                  end do
               end do
               
               bufn2 = bufn2+DIR_COSE_NSJ(N1,S,J1,N2,S,J2)*bufk3
            end do
         end do
         mag1=bufn2*nucm_cgs
      end if
      
      rotsc = trans_dips1*mag1*dircoseM+trans_dips1*magel2
   end subroutine TRANS_ROTSTR_SR_CART_FROM_SPH_NEW
   
   
   subroutine TRANS_ROTSTR_SR_CART_FROM_SPH(c2,c1,J2,J1,S,mult,rotsc,dip,gt,trans_dips1)
      double precision, intent(in) :: J1,J2
      type(TUPLE_EVEC_SPINROT), intent(in) :: c2,c1
      double precision, intent(in) :: S,dip(3),gt(3,3)
      double complex,intent(out) :: rotsc
      double complex,intent(in) :: trans_dips1(NINT(2*J2+1))
      integer, intent(in) :: mult
      
      integer :: N2_idx,N1_idx,N2,N1,mn3,Mn2,Mn1,Ms2_idx,M2_idx,M1_idx,K2,K1,a,i,n3l,n3u,k3,n3,ks1_idx,N_idx,F,diffJ
      double precision :: ms1,Ms2,M2,M1,cg1,cg2,kj1,ks1,ks3
      double complex,allocatable :: evec2(:),evec1(:)
      character :: labaxis 
      character, parameter :: axes(3) = ['X','Y','Z']
      double complex :: bufm1,bufm2,bufk1,bufk2,bufk3,bufn1,bufn2,bufsK,bufa,bufs,bufm3,magel,magel2,magel2m
      double complex :: ga(3,3),gcor(3),mag1,mag2,contrs_F(3),contrs_F_magel(3),contrs_F_magnuc(3)
      
      
      do i = 1,3
         ga(i,1)=(gt(i,1)-iu*gt(i,2))*0.5d0
         ga(i,2)=gt(i,3)
         ga(i,3)=(gt(i,1)+iu*gt(i,2))*0.5d0
      end do
       
      gcor(1)=0.5d0*iu*(gt(2,3)-gt(3,2))
      gcor(2)=0.5d0*iu*(gt(3,1)-gt(1,3))
      gcor(3)=0.5d0*iu*(gt(1,2)-gt(2,1))
      rotsc=0d0
      magel2m=0d0
      contrs_F=0d0
      contrs_F_magel=0d0
      contrs_F_magnuc=0d0
      diffJ=abs(NINT(J2-J1))
      !bufsk=0d0
      !do F=1,3
      labaxis=axes(3)
      do M2_idx = 0,NINT(2*J2)
         M2=M2_idx-J2
      ! do m1_idx = max(m2_idx-1,0),min(m2_idx+1,NINT(2*J1))
      !do m1_idx = 0,NINT(2*J1)
      !do m1_idx = max(m2_idx-1+diffJ,0),min(m2_idx+1+diffJ,NINT(2*J1))
         M1=m2
         
         bufn2=0d0
         bufn1=0d0
         magel=0d0
         magel2m=0d0
         do N2_idx = 1,mult
            N2=c2%N_arr(N2_idx)
            evec2=conjg(c2%Ns(N2_idx)%arr)
            N_idx=findloc(c1%N_arr,N2,dim=1)
            if(N_idx>0 .and. magdip_contr(2))then
               evec1=c1%Ns(N_idx)%arr
               bufs=0d0
               !do mn2=-N2,N2
                  do Ms2_idx=0,NINT(2*S)
                     Ms2=Ms2_idx-S
                     mn2=M2-Ms2
                     Ms1=M1-M2+Ms2
                     cg1=CG_SR(N2,S,J2,Mn2,Ms2,M2)
                     cg2=CG_SR(N2,S,J1,Mn2,Ms1,M1)
                     bufs=bufs+cg1*cg2*SAlpha_matel_space(S,Ms1,S,Ms2,labaxis)
                  end do
               !end do
               
               bufsk=0d0
               do K2=-N2,N2
                  bufsk = bufsk+conjg(evec2(N2+K2+1))*conjg(evec1(N2+K2+1))
               end do
               magel=magel+bufsk*bufs
            end if
            
            do N1_idx = 1,mult
               N1=c1%N_arr(N1_idx)
               if(N1+N2==0)cycle
               evec1=c1%Ns(N1_idx)%arr
               
               bufm2=0d0
               !do Mn2 = -N2,N2
                  !Mn1=NINT(M1-M2)+Mn2
                  do Ms2_idx=0,NINT(2*S)
                     Ms2=Ms2_idx-S
                     Mn2=NINT(M2-Ms2)
                     Mn1=NINT(M1-Ms2)
                     cg1=CG_SR(N2,S,J2,Mn2,Ms2,M2)
                     cg2=CG_SR(N1,S,J1,Mn1,Ms2,M1)
                     bufm2=bufm2+cg1*cg2*DIR_COSE_NM_FROM_SPH(N1,Mn1,N2,Mn2,labaxis)
                  end do
               !end do
               
               
               bufk3=0d0
               do K2 = -N2,N2
                  do K1 = max(-N1,K2-2),min(N1,K2+2)
                     bufk2=0d0
                     !nuclear magnetic moment in molecule-fixed axes
                     do a = 1,3
                        bufK2=bufk2+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2+1,axes(a))*Jminus_eval(N2,K2)*ga(a,3)+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2-1,axes(a))*Jplus_eval(N2,K2)*ga(a,1)+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2,axes(a))*K2*ga(a,2)
                        bufk2=bufk2+DIR_COSE_NK_FROM_SPH(N1,K1,N2,K2,axes(a))*gcor(a)
                     end do
                     bufk3=bufk3+bufk2*evec1(N1+K1+1)*evec2(N2+K2+1)
                  end do
               end do
               
               bufn2 = bufn2+DIR_COSE_N_FROM_SPH(N1,N2)*bufk3*bufm2
            end do
         end do
         
         mag1=bufn2*nucm_cgs
         mag2=magel*bohrm_cgs*gs
         
        ! contrs_F(F) = contrs_F(F) + trans_dips1(F,M2_idx+1,M1_idx+1)*(mag1+mag2)
         rotsc = rotsc + trans_dips1(M2_idx+1)*(mag1+mag2)
         !contrs_F_magnuc(F) = contrs_F_magnuc(F) + trans_dips1(F,M2_idx+1,M1_idx+1)*(mag1)
         !contrs_F_magel(F) = contrs_F_magel(F) + trans_dips1(F,M2_idx+1,M1_idx+1)*(mag2)
      end do
      !end do
      !end do
      !rotsc=sum(contrs_F_magel+contrs_F_magnuc)
   end subroutine TRANS_ROTSTR_SR_CART_FROM_SPH
   
   subroutine NKSJ_FillTrans(tr,ener2,ener1,N2,N1,T2,T1,S2,S1,J2,J1,dipstr,rotstr,freq,carrr)
      type(TRANS_NTSJ) tr
      INTEGER N2,N1,T2,T1
      Integer,intent(inout) :: carrr
      double precision S2,S1,J2,J1,dipstr,rotstr,freq,ds,rs,rd
      double complex :: ener2,ener1
      
      carrr=carrr+1
      if(freq>0d0)then
         tr%N2=N2
         tr%N1=N1
         tr%J2=J2
         tr%J1=J1
         tr%S2=S2
         tr%S1=S1
         tr%T2=T2
         tr%T1=T1
         tr%dipstr=dipstr
         tr%rotstr=rotstr
         tr%freq=freq
         tr%en2=ener2
         tr%en1=ener1
      else
         tr%N2=N1
         tr%N1=N2
         tr%J2=J1
         tr%J1=J2
         tr%S2=S1
         tr%S1=S2
         tr%T2=T1
         tr%T1=T2
         tr%dipstr=dipstr
         tr%rotstr=rotstr
         tr%freq=abs(freq)
         tr%en2=ener1
         tr%en1=ener2
      end if
      if(file_transitions(7).and.file_transitions(1))then
         ds=tr%dipstr*AU2SI_debye**2
         rs=tr%rotstr*AU2SI_debye**2
         if(.not.file_transitions(8).and.(abs(rs)<outtol.or. ds<outtol))return
         if(ds/=0d0)then
            RD=4d0*rs/ds
         end if
6000     format(4(A15,1X),2(3A4,1A6))
6001     format(4(E15.8E2,1X),2(2I4,1A4,1A6))
         write(unitt,6001)tr%freq,ds,rs,RD,tr%N1,tr%T1,R2HI(tr%S1),R2HI(tr%J1),tr%N2,tr%T2,R2HI(tr%S2),R2HI(tr%J2)
      end if
   end subroutine NKSJ_FillTrans
   
   subroutine CALCULATE_RCD_Jpair(J1,J2)
      INTEGER J1,J2,T1,T2,idx1,idx2,m1,m2,groundKs(3),excitedKs(3),p,F
      double precision,allocatable :: evec1(:,:),evec2(:,:),eval1(:),eval2(:)
      double precision rots,dips,dipsn,rotsn,ener,rotss,rotsx,dipsx,rotsy,rotsqx,rotsqy,rotsq,rotsmsum,rotsqsum
      double complex carr(6),iu,quadred,rotx,dipx
      double precision darr(9),bufr,avgMM,cosM,cosM_sph,cosJ,dips_contr(4),rots_contr(4),rots_contr_sum,dips_contr_sum,bufm,dips_sph,rots_sph
      double complex dipel1,dipel2,dipnuc1,dipnuc2,dipsc,rotsc,magel2,magnuc2,dipsc_sph,rotsc_sph,trans_dip
      double precision cosj_sw
      double complex cosm_sw
      type(TRANS_JT),pointer :: tr
      character(1),parameter :: axes(3)=['X','Y','Z']
      
      
      !1-lower
      !2-higher
      allocate(evec1(2*J1+1,2*J1+1),evec2(2*J2+1,2*J2+1),eval1(2*J1+1),eval2(2*J2+1))
      evec1=0.0d0
      evec2=0.0d0
      eval1=0.0d0
      eval2=0.0d0
      call SLICE_EVEC(evec1,J1)
      call SLICE_EVEC(evec2,J2)
      call SLICE_EVAL(eval1,J1)
      call SLICE_EVAL(eval2,J2)
            
      !avgMM=1.0d0
      cosJ = dble(DIR_COSE_N_FROM_SPH(J1,J2)*DIR_COSE_N_FROM_SPH(J2,J1))
      cosM=0d0
      do F = 1,2
         do m1=-J1,J1
            do m2=max(m1-1,-J2),min(m1+1,J2)
               cosM = cosM + dble(DIR_COSE_NM_FROM_SPH(J1,M1,J2,M2,axes(F))*DIR_COSE_NM_FROM_SPH(J2,M2,J1,M1,axes(F)))
            end do
         end do
      end do
      
      
      !lower state
      if(separateElNuc)then
         do T1=-J1,J1
            !higher state
            do T2=-J2,J2
               if((J1==J2).AND.(T1>=T2))cycle
               tr=>trans_jt_col(car)
               ener=eval2(T2+J2+1)-eval1(T1+J1+1)
               if(SkipTransFlag.and.(abs(ener)<fmin .or. abs(ener)>fmax))then
                  cycle
               end if
               if(skipTransFlagLowPop)then
                  if(((2*J2+1)*exp(-h*eval2(T2+J2+1)*beta)-(2*J1+1)*exp(-h*eval1(T1+J1+1)*beta))/Q < LowPopTol)cycle
               end if
               rotsc=0d0
               dipsc=0d0
               if(ener>0d0)then
                  avgMM=1.0d0/dble(2*J1+1)
                  dipel1 =SmallDips(evec1,evec2,J1,J2,T1,T2,dip_el)
                  dipel2 =SmallDips(evec2,evec1,J2,J1,T2,T1,dip_el)
                  dipnuc1=SmallDips(evec1,evec2,J1,J2,T1,T2,dip_nuc)
                  dipnuc2=SmallDips(evec2,evec1,J2,J1,T2,T1,dip_nuc)
                  dips_contr=dble([dipel1*dipel2,dipel1*dipnuc2,dipnuc1*dipel2,dipnuc1*dipnuc2])*cosJ*cosM*avgMM
                  
                  magel2=SmallRots(evec2,evec1,J2,J1,T2,T1,ge)
                  magnuc2=SmallRots(evec2,evec1,J2,J1,T2,T1,gn)
                  rots_contr=aimag([dipel1*magel2,dipel1*magnuc2,dipnuc1*magel2,dipnuc1*magnuc2])*nucm_cgs*cosJ*cosM*avgMM
                  
                  call TRANS_DIPSTR_CART_FROM_SPH(evec1(:,T1+J1+1),evec2(:,T2+J2+1),J1,J2,dipsc,dip,trans_dip)
                  call TRANS_ROTSTR_CART_FROM_SPH(evec1(:,T1+J1+1),evec2(:,T2+J2+1),J1,J2,rotsc,dip,gt,trans_dip)
                  rots=nucm_cgs*cosJ*cosM*avgMM*aimag(rotsc)
                  dips=cosJ*cosM*avgMM*dble(dipsc)
               else
                  avgMM=1.0d0/dble(2*J2+1)
                  dipel1 =SmallDips(evec2,evec1,J2,J1,T2,T1,dip_el)
                  dipel2 =SmallDips(evec1,evec2,J1,J2,T1,T2,dip_el)
                  dipnuc1=SmallDips(evec2,evec1,J2,J1,T2,T1,dip_nuc)
                  dipnuc2=SmallDips(evec1,evec2,J1,J2,T1,T2,dip_nuc)
                  dips_contr=dble([dipel1*dipel2,dipel1*dipnuc2,dipnuc1*dipel2,dipnuc1*dipnuc2])*cosJ*cosM*avgMM
                  
                  magel2=SmallRots(evec1,evec2,J1,J2,T1,T2,ge)
                  magnuc2=SmallRots(evec1,evec2,J1,J2,T1,T2,gn)
                  rots_contr=aimag([dipel1*magel2,dipel1*magnuc2,dipnuc1*magel2,dipnuc1*magnuc2])*nucm_cgs*cosJ*cosM*avgMM
                  
                  call TRANS_DIPSTR_CART_FROM_SPH(evec2(:,T2+J2+1),evec1(:,T1+J1+1),J2,J1,dipsc,dip,trans_dip)
                  call TRANS_ROTSTR_CART_FROM_SPH(evec2(:,T2+J2+1),evec1(:,T1+J1+1),J2,J1,rotsc,dip,gt,trans_dip)
                  rots=nucm_cgs*cosJ*cosM*avgMM*aimag(rotsc)
                  dips=cosJ*cosM*avgMM*dble(dipsc)
               end if
               dips_contr_sum=sum(dips_contr)
               rots_contr_sum=sum(rots_contr)
               if(ener>0d0)then
                  tr%J1=j1
                  tr%J2=j2
                  tr%t1=t1
                  tr%t2=t2
                  tr%en1=eval1(T1+J1+1)
                  tr%en2=eval2(T2+J2+1)
               else
                  tr%J1=j2
                  tr%J2=j1
                  tr%t1=t2
                  tr%t2=t1
                  tr%en1=eval2(T2+J2+1)
                  tr%en2=eval1(T1+J1+1)
               end if
               tr%freq=abs(ener)
               tr%dipstr=dips
               tr%rotstr=rots
               tr%dips_contr=dips_contr
               tr%rots_contr=rots_contr
               car=car+1
      20       continue
            end do
         end do
      else
         do T1=-J1,J1
            !higher state
            do T2=-J2,J2
               if((J1==J2).AND.(T1>=T2))cycle
               tr=>trans_jt_col(car)
               if(skipTransFlagLowPop)then
                  if(((2*J2+1)*exp(-h*eval2(T2+J2+1)*beta)-(2*J1+1)*exp(-h*eval1(T1+J1+1)*beta))/Q < LowPopTol)cycle
               end if
               rots=0d0
               dips=0d0
               ener=eval2(T2+J2+1)-eval1(T1+J1+1)
               
               if(SkipTransFlag.and.(abs(ener)<fmin .or. abs(ener)>fmax))then
                  cycle
               end if

               rotsc=(0d0,0d0)
               !rotquadsc1=(0d0,0d0)
               dipsc=(0d0,0d0)
               rotx=0d0
               dipx=0d0
               if(ener>0d0)then
                  avgMM=1.0d0/dble(2*J1+1)
                  !call TRANS_DIPSTR_SPH(evec1(:,T1+J1+1),evec2(:,T2+J2+1),J1,J2,dipsc,T_dip1)
                  call TRANS_DIPSTR_CART_FROM_SPH(evec1(:,T1+J1+1),evec2(:,T2+J2+1),J1,J2,dipsc,dip,trans_dip)
                  
                  !call TRANS_ROTSTR(evec1,evec2,J1,J2,T1,T2,rotsc,dip,gt)
                  call TRANS_ROTSTR_CART_FROM_SPH(evec1(:,T1+J1+1),evec2(:,T2+J2+1),J1,J2,rotsc,dip,gt,trans_dip)
               else
                  avgMM=1.0d0/dble(2*J2+1)
                  !call TRANS_DIPSTR_SPH(evec2(:,T2+J2+1),evec1(:,T1+J1+1),J2,J1,dipsc,T_dip1)
                  call TRANS_DIPSTR_CART_FROM_SPH(evec2(:,T2+J2+1),evec1(:,T1+J1+1),J2,J1,dipsc,dip,trans_dip)
                  
                  !call TRANS_ROTSTR(evec1,evec2,J2,J1,T2,T1,rotsc,dip,gt)
                  call TRANS_ROTSTR_CART_FROM_SPH(evec2(:,T2+J2+1),evec1(:,T1+J1+1),J2,J1,rotsc,dip,gt,trans_dip)
               end if
               !dipsx=2*dble(dipx)
               !rotsx=2*aimag(rotx)
               dips=cosJ*cosM*avgMM*dble(dipsc)
               rots=nucm_cgs*cosJ*cosM*avgMM*aimag(rotsc)!-3d0/2d0*(2*J1+1)**(-1)*ener*1d9*SI2AU_Hz*2*pi*dble(rotquadsc1)
6733           continue    
               !car = Calc_storage_index(J1,T1,J2,T2)
               !!$OMP CRITICAL
               if(ener>0d0)then
                  tr%J1=j1
                  tr%J2=j2
                  tr%t1=t1
                  tr%t2=t2
                  tr%en1=eval1(T1+J1+1)
                  tr%en2=eval2(T2+J2+1)
               else
                  tr%J1=j2
                  tr%J2=j1
                  tr%t1=t2
                  tr%t2=t1
                  tr%en1=eval2(T2+J2+1)
                  tr%en2=eval1(T1+J1+1)
               end if
               tr%freq=abs(ener)
               tr%dipstr=dips
               tr%rotstr=rots
               car=car+1
               !!$OMP END CRITICAL
            end do
         end do
      end if
      !!$OMP END DO
      !!$OMP END PARALLEL
      deallocate(evec1,evec2,eval1,eval2)
   end subroutine CALCULATE_RCD_Jpair
   
   function SmallDips(cg,ce,J1,J2,T1,T2,dip)result(res)
      INTEGER j1,j2,t1,t2,i,k1,k2
      double complex res,buf1,buf2
      double precision cg(2*J1+1,2*J1+1),ce(2*J2+1,2*J2+1),cgb(2*J1+1),ceb(2*J2+1),dip(3)
      character(1) axes(3)
      parameter(axes=['X','Y','Z'])
      
      cgb=cg(:,T1+J1+1)
      ceb=ce(:,T2+J2+1)
      buf1=0d0
      do K1 = -J1,J1
         do K2 = max(K1-1,-J2),min(K1+1,J2)
            do i = 1,3
               buf1=buf1+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2,axes(i))*dip(i)*cgb(K1+J1+1)*ceb(K2+J2+1)
            end do
         end do
      end do 
      
      res = buf1
   end function SmallDips
   
   function SmallRots(cg,ce,J1,J2,T1,T2,gten)result(res)
      INTEGER j1,j2,t1,t2,i,k1,k2,j,k
      double complex res,buf2,buf3,gcor(3),ga(3,3)
      double precision cg(2*J1+1,2*J1+1),ce(2*J2+1,2*J2+1),cgb(2*J1+1),ceb(2*J2+1),gten(3,3)
      character(1) axes(3)
      parameter(axes=['X','Y','Z'])
      
      cgb=cg(:,T1+J1+1)
      ceb=ce(:,T2+J2+1)
      
      do i = 1,3
         ga(i,1)=(gten(i,1)-iu*gten(i,2))*0.5d0
         ga(i,2)=gten(i,3)
         ga(i,3)=(gten(i,1)+iu*gten(i,2))*0.5d0
      end do
       
      gcor(1)=0.5d0*iu*(gten(2,3)-gten(3,2))
      gcor(2)=0.5d0*iu*(gten(3,1)-gten(1,3))
      gcor(3)=0.5d0*iu*(gten(1,2)-gten(2,1))
      
      buf2=(0d0,0d0)
      buf3=(0d0,0d0)
      do K1 = -J1,J1
         do K2 = max(K1-2,-J2),min(K1+2,J2)
            buf2=0d0
            do i = 1,3
               buf2=buf2+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2+1,axes(i))*Jminus_eval(J2,K2)*ga(i,3)+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2-1,axes(i))*Jplus_eval(J2,K2)*ga(i,1)+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2,axes(i))*K2*ga(i,2)
               buf2=buf2+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2,axes(i))*gcor(i)
            end do
            buf3=buf3+buf2*cgb(K1+J1+1)*ceb(K2+J2+1)
         end do
      end do
      
      res=buf3
   end function SmallRots
   
   !get eigenvector matrix for given J
   subroutine SLICE_EVEC(curHam,J)
      Integer,INTENT(IN) :: J
      double precision, INTENT(OUT) :: curHAM(2*J+1,2*J+1)

      
      if(J>endJ) STOP 'J higher than maximum allowed J'

      
      curHAM=Evec_full_list(J+1)%arr
   end subroutine SLICE_EVEC

   !get eigenvalues for given J
   subroutine SLICE_EVAL(eval,J)
      Integer,INTENT(IN) :: J
      double precision, INTENT(OUT) :: eval(2*J+1)
      
      if(J>endJ) STOP 'J higher than maximum allowed J'
      
      eval=eval_full_list(J+1)%arr
   end subroutine SLICE_EVAL

   subroutine TRANS_DIPSTR_CART_FROM_SPH(cg,ce,J1,J2,res,dip,trans_dip)
      INTEGER j1,j2,k1,k2,alpha
      double complex res,buf1,buf2,buf3,buf4,trans_dip
      double precision cg(2*J1+1),ce(2*J2+1),dip(3)
      character(1),parameter :: axes(3)=['X','Y','Z']
      
      buf1=(0d0,0d0)
      buf2=(0d0,0d0)
      do k1=-J1,J1
         do k2=max(-J2,k1-1),min(J2,k1+1)
            buf3=0d0
            buf4=0d0
            do alpha = 1,3
               buf3=buf3+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2,axes(alpha))*dip(alpha)
               buf4=buf4+DIR_COSE_NK_FROM_SPH(J2,K2,J1,K1,axes(alpha))*dip(alpha)
            end do
            buf1=buf1+cg(K1+J1+1)*ce(K2+J2+1)*buf3
            buf2=buf2+cg(K1+J1+1)*ce(K2+J2+1)*buf4
         end do
      end do
      
      trans_dip=buf1
      res=buf1*buf2
   end subroutine TRANS_DIPSTR_CART_FROM_SPH
   
   subroutine TRANS_ROTSTR_CART_FROM_SPH(cgb,ceb,J1,J2,res,dip,gt,trans_dip)
      INTEGER j1,j2,t1,t2,i,k1,k2,j,k
      double complex res,buf1,buf2,buf3,gcor(3),ga(3,3),trans_dip
      double precision cgb(2*J1+1),ceb(2*J2+1),gt(3,3),dip(3)
      character(1) axes(3)
      parameter(axes=['X','Y','Z'])
      
      
      do i = 1,3
         ga(i,1)=(gt(i,1)-iu*gt(i,2))*0.5d0
         ga(i,2)=gt(i,3)
         ga(i,3)=(gt(i,1)+iu*gt(i,2))*0.5d0
      end do
       
      gcor(1)=0.5d0*iu*(gt(2,3)-gt(3,2))
      gcor(2)=0.5d0*iu*(gt(3,1)-gt(1,3))
      gcor(3)=0.5d0*iu*(gt(1,2)-gt(2,1))
      
      buf1=(0d0,0d0)
     ! buf2=(0d0,0d0)
      buf3=(0d0,0d0)
      do K1 = -J1,J1
         do K2 = max(K1-2,-J2),min(K1+2,J2)
            buf2=0d0
            do i = 1,3
               !buf1=buf1+DIR_COSE_NK_FROM_SPH(J1,K1,J2,K2,axes(i))*dip(i)*cgb(K1+J1+1)*ceb(K2+J2+1)
               buf2=buf2+DIR_COSE_NK_FROM_SPH(J2,K2,J1,K1+1,axes(i))*Jminus_eval(J1,K1)*ga(i,3)+DIR_COSE_NK_FROM_SPH(J2,K2,J1,K1-1,axes(i))*Jplus_eval(J1,K1)*ga(i,1)+DIR_COSE_NK_FROM_SPH(J2,K2,J1,K1,axes(i))*K1*ga(i,2)
               buf2=buf2+DIR_COSE_NK_FROM_SPH(J2,K2,J1,K1,axes(i))*gcor(i)
               
            end do
            buf3=buf3+buf2*cgb(K1+J1+1)*ceb(K2+J2+1)
         end do
      end do
      
      res=trans_dip*buf3
   end subroutine TRANS_ROTSTR_CART_FROM_SPH
           
end program main 
