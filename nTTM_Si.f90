!    Copyright 2020 Vladimir Lipp, Baerbel Rethfeld, Martin Garcia, Dmitry Ivanov
!    Contact: v.p.lipp@gmail.com

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.

Program nTTM

!       This code solves the system of continuum equations implementing the nTTM model of short pulse laser interaction with Si
!	The implemented finite-difference scheme is a modified Crank-Nicolson with a predictor-corrector. The scheme, chosen set of parameters, and material properties are published in [Lipp, Rethfeld, Garcia, Ivanov, Appl. Sci., 2020]
!	Equation numbers correspond to the above publication

!	The code has been written as a part of doctoral research of Vladimir Lipp between 2009 and 2015
!	For more details, see his dissertation "Atomistic-continuum modeling of ultrafast laser-induced melting of silicon targets" (University of Kassel): https://kobra.uni-kassel.de/handle/123456789/2016031650027

!	The  code is designed for GNU/Linux 
!	To compile with gfortran, enter ./compile.sh
!	After compilation, to run the program, enter ./nTTM.x
!	After the calculation is done, run ./plot_em_all.sh to generate plots in png format using gnuplot


        IMPLICIT NONE

        real(8) dI_dz,Source_t,k_l,K_c,C_l,tau_e,teta,Refl, &    ! functions
             alpha,theta,EnergyGap,EnergyGap_n,EnergyGap_T_l

        EXTERNAL dI_dz,Source_t,K_l,K_c,C_l,tau_e,teta,Refl, & ! functions
             alpha,theta,EnergyGap,EnergyGap_n,EnergyGap_T_l
        
!       Material and laser input parameters
        real(8), parameter :: gamma = 3.8d-43                         ! [m6 s-1]           Auger recombination coefficient
        !real(8), parameter :: beta  = 15.d-11                           ! [m W-1]            two-phonon absorption coefficent 
		real(8), parameter :: beta  = 2.0d-11
        real(8), parameter :: me    = 9.10956d-31                     ! [kg]                 electron's mass
        real(8), parameter :: mee   = 0.36d0*me                       ! [kg]                 effective electron mass
        real(8), parameter :: meh   = 0.81d0*me                       ! [kg]                 effective hole mass
        real(8), parameter :: pi    = 3.1415926535897932384626433832795d0 
        real(8), parameter :: hP    = 1.05459d-34                       ! [J s]                Planck's constant
        real(8), parameter :: kb    = 1.3806504d-23   ! [J K-1]                                 Boltzmann constant
        real(8), parameter :: el    = 1.602176487d-19   ! [C]
        real(8), parameter :: mue   = 0.0085d0      ! [m2 V-1 s-1]         ! mobility of electrons
        real(8), parameter :: muh   = 0.0019d0      ! [m2 V-1 s-1]         ! mobility of holes
        real(8), parameter :: las   = 1.55d0             ! [eV] for 800 nm               photon energy in eV
        real(8), parameter :: hV    = las*el      ! [J]                                         photon energy

        integer, parameter :: MP = 200                                                  !      maximum Number of connections between cells
        integer, parameter :: SP = MP + 1                                              !       maximum Number of cells 
        integer, parameter :: FD = 1000000                                              !       Maximum lenght of Fermi functions

        real(8), parameter :: energ_den_melting    = 3859199787.d0      ! [J/m3]                          		! Found "experimentally" by running this code
	
        !real(8), parameter :: tpulse = 130.d-15      ! [s]       pulse duration defined as the full width at the half maximum
		real(8), parameter :: tpulse = 500.d-15 		! [s] pulse duration defined as full width at half maximum in Chen paper

!       Arrays of Variables
        real(8) k_carriers_left, k_carriers_right ! carriers heat conductivity
        real(8) Den_new(SP),Den_old(SP),T_ll(SP),Ta_old(SP),T_l_new(SP),X(SP)
        real(8) T_e(SP),Te_old(SP),T_e_new(SP),IL(SP)                                    ! electron temperature, the one at the next time step, 
        real(8) Egap_old(SP),Egap_n(SP),Egap_T_l(SP),Egap_new(SP),Egap_n_new(SP),Egap_T_l_new(SP)
        real(8) W(MP),Ju_old(MP),Ju_new(MP),D(MP)                                             ! ambipolar energy current, particles current, diffusion coefficient
        real(8) dTe(MP),dTl(MP),dDn(MP)
        real(8) ATe(MP),ATl(MP),ADn(MP),AKl(MP)

        real(8) F_32(0:FD),F_1(0:FD),F_12(0:FD),F0(0:FD),F12(0:FD),F1(0:FD),F32(0:FD)
        real(8) muni(0:FD)
        real(8) theta_e(SP),H10e(SP),H12_12e(SP),H012e(SP),H3212e(SP),H3212old(SP)    ! reduced chemical potential of electrons, 
	real(8) H3212e_old(sp), H3212h_old(SP), H12_12e_old(sp), H12_12h_old(sp)
        real(8) theta_h(SP),H10h(SP),H12_12h(SP),H012h(SP),H3212h(SP)    ! reduced chemical potential of holes,
        real(8) AH012e(MP),AH012h(MP),AH10(MP),AH10_pre(MP),AH12_12(MP)

        real(8) muni0,munif,dmuni,Fhalfe,Fhalfh,kh1,mstar
        real(8) eta_emin,eta_emax,eta_hmin,eta_hmax                            ! maximum and minimum of reduced Fermi levels for electrons and holes

        real(8) F_total,LaserTTM,difz,t_curr,dz2,dtz,temm,denn,Istar,IL0
        real(8) F_las(SP),F_las_old(SP),F_las_new(SP),E_total,E_total_0,E_total_e,E_total_l,E_base,e_energ_den(SP),l_energ_den(SP),tot_energ_den(SP),tot_energ_dens
        real(8) dn_dt(SP),dn_dt_old(SP),dn_dt_new(SP),dT_l_dt(SP),dT_l_dt_old(SP),dT_l_dt_new(SP),dT_e_dt,C_e_h(SP),C_e_h_old(SP),C_e_h_new(SP),C_lat,t_rel
	real(8) dtheta_e_dn_old(SP), dtheta_h_dn_old(SP),dtheta_e_dn_new(SP), dtheta_h_dn_new(SP), theta_e_left, theta_e_right, theta_h_left, theta_h_right, delta_den
	real(8) dtheta_e_dTe_old(SP), dtheta_h_dTe_old(SP),dtheta_e_dTe_new(SP), dtheta_h_dTe_new(SP)
        real(8) FF12e,FF_12e,FF1e,FF0e,FF32e		!  F(1/2,e), F(-1/2,e), F(1,e), F(3/2,e)
        real(8) FF12h,FF_12h,FF1h,FF0h,FF32h		!  F(1/2,h), F(-1/2,h), F(1,h), F(3/2,h)
	real(8) E_nonthermal(SP)		! see Korfiatis, J Appl Phys 40 6803-6808 (2007) and Lipp, Ivanov, Rethfeld, Garcia, Journal of Optical Technology 81, 254 (2014)

        real(8) timen,timef,t1,t2,E_add
	
	real(8) Etrans(SP), E_pere

        integer EKk(SP),HK(SP)              !     array of indexes of F12 function which correspond to current Fhalfe value for every cell , the same for Fhalfh

        integer N,I,J,Nini,IK,JK,JMIN
	integer(8) NstepTTM

        integer(8) q, corr_step  ! every q steps records in file 22 and screen printings are being done
	real(8) F_12e(SP), F_12h(SP),F_12e_new(SP), F_12h_new(SP)

        real(8) L,T_ini,D_ini,t_exit,dift
    	real(8) psi, corrector_precision
        integer Lcouple
	
	real(8) :: almost_zero=1.d-50
    
!===========Crank-Nicolson scheme variables================================        
        real(8) a_cn(SP),b_cn(SP),c_cn(SP),r_cn(SP)
        real(8) F_cn, f_old(SP), summ
!============/Crank-Nicolson scheme variables===============================

	logical ND_case								! if the system is supposed to be non-degenerate (with Boltzmann distribution of electrons)

	ND_case=.false. 	! non-degenrate case is not thoroughly tested
	
	psi = 0.5d0 		! explicit method: psi = 0; implicit method: psi = 1; semi-implicit: psi = 0.5

!       Initial parameters of the material
!		changing the parameters to plot figure 2 from Chen and Beraun 
        !L       =  0.80d-06                                     ! Thickness of the target film   [m]
		L = 20.0d-06											!thickness of target film [Chen and Beraun]
        T_ini   = 300.d0                                      ! Initial temperature     [K] 
        !D_ini   =   1.d+16                                   ! Initial density         [m-3]
		D_ini = 1.d+16
        N       = 160                                             ! Number of grid points in 1d-space
        dift      =   1.d-17                                    ! Time step               [s]
        !t_exit =  10.d-12                                    ! Time of the simulation  [s]
		!t_exit = 5.0d-12										! time duration of simulation for case 1 and case 2 in chen and Beraun
		t_exit = 20.0d-12										!time duration for case 2 simulation in chen and beraun 
        Lcouple = 1

       	corrector_precision = 1.d-6
	

!       Some additional parameters
        kh1   = 2.d0**0.3333333333333333333d0*pi*hP*hP/kb    ! part of eq. (2)
        mstar = (mee/meh)**1.5d0                          ! ratio of effective e- and h-masses in the power of 3/2
	
	    t1 = 0.d0                                                       ! variable for estimation of real calculation time
	    q = int(1.d-15/dift)
        if (q.lt.1) q=1 

	Etrans(1:sp) = 0.d0	! energy transferred between electrons and lattice

        call CPU_TIME(timen)

        open(UNIT = 1,FILE = 'contour.dat')              ! 3d graphs of Te,Tl,Den depending on (t_curr,x)
	write(1,*) "#time (ps),    coordinate (mkm),    electron temperature (K), lattice temperature (K), carriers density (1e26 1/m^3)"
        open(UNIT = 112,FIlE = 'surf.dat')                    ! enetgies and temperatures depending on time
        open(UNIT = 4,FILE = 'fermi.tab')                  ! the file with calculated Fermi-integrals
        open(UNIT = 5,FILE = 'potential.dat')              ! for reduced Fermi levels of e and h
        open (UNIT = 22,FILE='res.dat')
	write (22,'(17(A,4x))') "time [s]","Te [K]", "Ta [K]","n [m^-3]" &
						, "Energ_dens_at_surf", "energ_dens_e_surf", "Etot [m^-2]" &
						, "Energ_dens_tot_surf", "Energ_dens_tot"  &
						,"Te(2) [kK]", "Ta(2) [kK]" &
						, "n(2) [1e26 m^-3]" &
						, "Energ_nonthermal" &
						, "Tot_fluence [m^-2]","E_tot_e [m^-2]", "E_tot_a [m^-2]", "E_gap_surf [J]"

!       Read the table of Fermi-Dirack Integral values        
        READ(4,*) muni0,munif,Nini

	if (Nini.GT.FD) then
           print *,"Dimension of Fermi Tab",Nini, &          ! if the program arrays are smaller than than number of calculated functions
                "exceeds the array limit",FD                     ! one must exceed the limits of the arrays
           stop
        endif
        
        read(4,*) (muni(I),F_32(I),F_1(I),F_12(I),F0(I),F12(I),F1(I),F32(I),I=0,Nini)  ! reading of precalculated Fermi-integrals       F = F (muni)
        eta_emin = munif         ! this will be real minimum of reduced Fermi level for electrons
        eta_emax = muni0        ! this will be real maximum of reduced Fermi level for electrons
        eta_hmin = munif         ! this will be real minimum of reduced Fermi level for holes
        eta_hmax = muni0        ! this will be real maximum of reduced Fermi level for holes

!       Determine the points for finite differences method
        difz = L/DBLE(N-1)                                     ! space grid step
        dz2 = 1.d0/(difz*difz)
        dtz = dift*difz

!       Set the initial temperature, density, and energies distribution and
!       Determine the initial value for Fermi-Dirac functions in all cells
        do I=1, N
              X(I)    = (I-1)*difz                       ! 1d space linear grid
              Ta_old(I)  = T_ini                            ! lattice temperature linear grid
              Te_old(I)  = T_ini                           ! electron temperature linear grid
              den_old(I)  = D_ini                           ! carriers density linear grid

           Fhalfe = den_old(I)*(kh1/(mee*Te_old(I)))**1.5d0                ! eq. (2)                   Fermi-Dirac integral 1/2 for electrons
           
           if (Fhalfe.LT.F12(0).OR.Fhalfe.GT.F12(Nini)) then
              print *,"chemical potential of electrons is out of bound"				! the precalculated Fermi tables do not cover this value of the chemical potential
              print *,I,den_old(I),Te_old(I),Fhalfe,F12(0),F12(Nini)
              stop
           endif

           if (I.EQ.1) then           ! finding the index J of array F12 knowing the value of the function at the first cell

              chem1: Do J = 0,Nini-1
                 if (Fhalfe.ge.F12(J).and.Fhalfe.lt.F12(J+1)) then
                    IK = J               ! storing the found index
                    exit chem1
                 endif
              enddo chem1

           else
           
! in the next cell the index should not change significantly so we search it near the previous value
                   do while (Fhalfe.GE.F12(IK+1))
                      IK = IK + 1
                   enddo
                   do while (Fhalfe.LT.F12(IK))
                      IK = IK - 1
                   enddo
           endif
           
           
		EKk(I) = IK

		theta_e(I) = muni(IK) + &
		(Fhalfe-F12(IK))*(muni(IK+1)-muni(IK))/(F12(IK+1)-F12(IK))       ! linear approximation fot evaluation of eta(e)

			if (ND_case) theta_e(i) = dlog(Fhalfe)

		FF12e  = Fhalfe
		! linear approximations for evaluation of F-functions
		FF_12e = F_12(IK) + (F_12(IK+1)-F_12(IK))*(theta_e(I) - &
		muni(IK))/(muni(IK+1)-muni(IK))
		FF0e   = F0(IK) + (F0(IK+1)-F0(IK))*(theta_e(I) - &
		muni(IK))/(muni(IK+1)-muni(IK))
		FF1e   = F1(IK) + (F1(IK+1)-F1(IK))*(theta_e(I) - &
		muni(IK))/(muni(IK+1)-muni(IK))
		FF32e  = F32(IK) + (F32(IK+1)-F32(IK))*(theta_e(I) - &
		muni(IK))/(muni(IK+1)-muni(IK))

			if (ND_case) then
				FF_12e=Fhalfe
				FF0e=Fhalfe
				FF1e=Fhalfe
				FF32e=Fhalfe
			endif

		F_12e(I)=ff_12e

		H10e(I)    = FF1e/FF0e                ! definitions of fuctions H
		H12_12e(I) = FF12e/FF_12e         !
		H012e(I)   = FF0e/FF12e             !
		H3212e(I)  = FF32e/FF12e           !

		Fhalfh = Fhalfe*mstar                          ! eq. (2)

           if (Fhalfh.LT.F12(0).OR.Fhalfh.GT.F12(Nini)) then
              print *,"chemical potential of holes is out of bound"				! the precalculated Fermi tables do not cover this value of the chemical potential
              print *,I,den_old(I),Te_old(I),Fhalfh,F12(0),F12(Nini)
              stop
           endif

           if (I.EQ.1) then
              chem2: Do J = 0,Nini-1
                 if (Fhalfh.ge.F12(J).and.Fhalfh.lt.F12(J+1)) then
                    JK = J
                    exit chem2
                 endif
              enddo chem2
           else
   
! in the next cell the index should not change significantly so we search it near the previous value
	   do while (Fhalfh.GE.F12(JK+1))
	      JK = JK + 1
	   enddo
	   do while (Fhalfh.LT.F12(JK))
	      JK = JK - 1
	   enddo
                   
           endif

		HK(I) = JK

		theta_h(I) = muni(JK) + &
		(Fhalfh-F12(JK))*(muni(JK+1)-muni(JK))/(F12(JK+1)-F12(JK))       ! linear approximation fot evaluation of eta(h)

		if (ND_case) theta_h(i) = dlog(Fhalfh)

		FF12h  = Fhalfh
		! linear approximations for evaluation of F-functions
		FF_12h = F_12(JK) + (F_12(JK+1)-F_12(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		FF0h   = F0(JK) + (F0(JK+1)-F0(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		FF1h   = F1(JK) + (F1(JK+1)-F1(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		FF32h  = F32(JK) + (F32(JK+1)-F32(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		
			if (ND_case) then
				FF_12h=Fhalfh
				FF0h=Fhalfh
				FF1h=Fhalfh
				FF32h=Fhalfh
			endif
		
	F_12h(I)=ff_12h
	
		H10h(I)    = FF1h/FF0h
		H12_12h(I) = FF12h/FF_12h
		H012h(I)   = FF0h/FF12h
		H3212h(I)  = FF32h/FF12h

		if (theta_e(I).LE.eta_emin) eta_emin = theta_e(I)
		if (theta_e(I).GE.eta_emax) eta_emax = theta_e(I)
		if (theta_h(I).LE.eta_hmin) eta_hmin = theta_h(I)
		if (theta_h(I).GE.eta_hmax) eta_hmax = theta_h(I)

		Egap_old(I)       = EnergyGap(Ta_old(I),den_old(I))
		Egap_n(I)    = EnergyGap_n(Ta_old(I),den_old(I))
		Egap_T_l(I)  = EnergyGap_T_l(Ta_old(I),den_old(I))

        enddo

!       Initial settings for the time loop
           t_curr = 0.d0		! current time
           F_total   = 0.d0		! total absorbed fluence
           E_add     = 0.d0		! energy transferred to atoms so far
	   l_energ_den(1:N) = C_l(T_ini)*Ta_old(1:N)	! only valid for homogeneous initial temperature
	   e_energ_den(1:N) = den_old(1:N)*(Egap_old(1:N) + 1.5d0*kb*Te_old(1:N)*H3212old(1:N))         ! (16)
	   tot_energ_dens = sum(l_energ_den(1:N)) + sum(e_energ_den(1:N))	! total energy density
	   
        NstepTTM = 0	! current time step

	! Initial energy of the electrons
	E_base = 0.5d0*den_old(1)*(Egap_old(1) + &
			1.5d0*kb*Te_old(1)*(H3212e(1)+H3212h(1)))*difz 
	E_base = E_base + 0.5d0*den_old(N)*(Egap_old(N) + &
			1.5d0*kb*Te_old(N)*(H3212e(N)+H3212h(N)))*difz 
	do i=2,N
	E_base = E_base + den_old(i)*(Egap_old(i) + &
			1.5d0*kb*Te_old(i)*(H3212e(i)+H3212h(i)))*difz 
	enddo
	
        E_total_0 = E_base + sum(C_l(T_ini)*Ta_old(1:N)*difz)		! only valid for homogeneous initial temperature
	
	
      write (*,'(7(A,4x))')  &
			"time [s]","Te [K]", "Ta [K]" &
			, "Energ_dens_at", "energ_dens_e", "n [m^-3]" &
			, "Etot [m^-2]"
			
					
        Time_loop: do

        NstepTTM = NstepTTM + 1
	E_total_e = 0.d0
	tot_energ_den(1:N) = 0.d0

	do i=1, N
	      Egap_old(I)     = EnergyGap(Ta_old(I),den_old(I))
	      Egap_n(I)   = EnergyGap_n(Ta_old(I),den_old(I))
	      Egap_T_l(I) = EnergyGap_T_l(Ta_old(I),den_old(I))
	enddo

!          Perform calculations of J and W for cell-connecting points
	DWJ_loop: do I = 1, N-1

!             Auxiliary multipyers
              dTe(I) = Te_old(I+1) - Te_old(I)
              dTl(I) = Ta_old(I+1) - Ta_old(I)
              dDn(I) = den_old(I+1) - den_old(I)

              ATe(I) = 0.5d0*(Te_old(I) + Te_old(I+1))
              ATl(I) = 0.5d0*(Ta_old(I) + Ta_old(I+1))
              ADn(I) = 0.5d0*(den_old(I) + den_old(I+1))
              AKl(I) = 0.5d0*(K_l(Ta_old(I)) + K_l(Ta_old(I+1)))

              AH012e(I)  = 0.5d0*(H012e(I) + H012e(I+1))
              AH012h(I)  = 0.5d0*(H012h(I) + H012h(I+1))
              AH10(I)    = 0.5d0*(H10e(I)  + H10e(I+1) + &
                   H10h(I) + H10h(I+1))
              AH12_12(I) = 0.5d0*(H12_12e(I) + H12_12e(I+1) + &
                   H12_12h(I) + H12_12h(I+1))

!             Diffusion coefficient
              D(I) = ATe(I)*kb*mue*muh*AH012e(I)*AH012h(I)*AH12_12(I)/ &	! eq. (23)
                   (mue*AH012e(I)+muh*AH012h(I))/el

!             carrier pair current
              Ju_old(I) = -D(I)*(dDn(I) + &
                   ADn(I)*(Egap_old(I+1)-Egap_old(I))/ATe(I)/kb/AH12_12(I) + &	! eq. (22)
                   ADn(I)*dTe(I)*(2.d0*AH10(I)/AH12_12(I)-1.5d0)/ATe(I))/difz
              
!             ambipolar energy current
              W(I) = (0.5d0*(Egap_old(I)+Egap_old(I+1)) + &
                   2.d0*kb*ATe(I)*AH10(I))*Ju_old(I) - &						! eq. (21)
                   (K_c(Te_old(I)) + K_c(Te_old(I+1)))*dTe(I)/difz
	      
           end do DWJ_loop


d_theta_by_d_nt: do I=1, N

	dtheta_e_dn_old(I) = .5d0 * ( 2.d0*pi*hp*hp/(mee*kb*Te_old(i)) )**(1.5d0) / F_12e(I)							! eq. (25)
	dtheta_h_dn_old(I) = .5d0 * ( 2.d0*pi*hp*hp/(meh*kb*Te_old(i)) )**(1.5d0) / F_12h(I)							! eq. (25)
	dtheta_e_dTe_old(I) = -3.d0/dsqrt(2.d0)*den_old(i)*Te_old(I)**(-2.5d0) * (pi*hp*hp/(mee*kb))**(1.5d0) / F_12e(I)		! eq. (26)
	dtheta_h_dTe_old(I) = -3.d0/dsqrt(2.d0)*den_old(i)*Te_old(I)**(-2.5d0) * (pi*hp*hp/(meh*kb))**(1.5d0) / F_12h(I)		! eq. (26)

enddo d_theta_by_d_nt



!          Perform calculations for the 1st node

!          Source term
	   if (t_curr .le. tpulse*10.d0) then
		IL0   = Source_t(t_curr+0.5d0*dift,Ta_old(1),tpulse)				! eq. (14)
		Istar = IL0    + 0.25d0*difz*dI_dz(IL0,Ta_old(1),den_old(1))
		denn  = den_old(1) + 0.25d0*dDn(1)
		temm  = Ta_old(1) + 0.25d0*dTl(1)
		IL(1) = IL0    +  0.5d0*difz*dI_dz(Istar,temm,denn)           
		LaserTTM = 0.5d0*(IL0 + IL(1))
	   else 
		IL0=0.d0
		Istar=0.d0
		denn  = 0.d0
		temm  = 0.d0
		IL(1) = 0.d0
		LaserTTM = 0.d0
	   endif

           C_e_h_old(1) = 1.5d0*kb*den_old(1)*(	H3212e(1) + H3212h(1) &                       ! eq. (10)
                   + Te_old(1)*dtheta_e_dTe_old(1)*(1.d0-H3212e(1)/H12_12e(1)) &              ! H12_12 = 1/( H_1212 )
                   + Te_old(1)*dtheta_h_dTe_old(1)*(1.d0-H3212h(1)/H12_12h(1)) 	)  ! + Den * dE_gap/dT_e

           C_lat = C_l(Ta_old(1))
           t_rel = tau_e(den_old(1))             ! electron relaxation time
           f_las_old(1) = (alpha(temm) + theta(temm)*denn)*LaserTTM + beta*LaserTTM*LaserTTM	! eqs. (13)

!          Density of the carriers
           dn_dt_old(1) = alpha(temm)*LaserTTM/hV + &
                0.5d0*beta*LaserTTM*LaserTTM/hV - &                   ! eq. (27)
                gamma*den_old(1)*den_old(1)*den_old(1) + &
                teta(Egap_old(1),Te_old(1))*den_old(1) - &
                2.d0*Ju_old(1)/difz

           Den_new(1) = den_old(1) + dn_dt_old(1)*dift

!          Lattice temperature
           dT_l_dt_old(1) = (2.d0*AKl(1)*dTl(1)*dz2 + &          		 ! eq. (28)
                C_e_h_old(1)*(Te_old(1)-Ta_old(1))/t_rel)/C_lat*Lcouple
                
           T_l_new(1) = Ta_old(1) + dT_l_dt_old(1)*dift



!          Perform calculations for T_e, T_l, and Den for internal cells
           Int_Loop: DO I = 2, N-1

!             Source term
	   if (t_curr .le. tpulse*10.d0) then
		Istar = IL(I-1) + 0.5d0*difz*dI_dz(IL(I-1),ATl(I-1),ADn(I-1))
		IL(I) = IL(I-1) + difz*dI_dz(Istar,Ta_old(I),den_old(I))
		LaserTTM = 0.5d0*(IL(I-1) + IL(I))
	   else
		Istar = 0.d0
		IL(I) = 0.d0
		LaserTTM = 0.d0
	   endif

!             Auxiliary definitions
              C_e_h_old(i) = 1.5d0*kb*den_old(I)*(H3212e(I) + H3212h(I) &                       ! eq. (10)
                   + Te_old(i)*dtheta_e_dTe_old(i)*(1.d0-H3212e(i)/H12_12e(i)) &              ! H12_12 = 1/( H_1212 )
                   + Te_old(i)*dtheta_h_dTe_old(i)*(1.d0-H3212h(i)/H12_12h(i)) 	)  ! + Den * dE_gap/dT_e
		   
              C_lat = C_l(Ta_old(I))
              t_rel = tau_e(den_old(I))
              f_las_old(I) = (alpha(Ta_old(I)) + theta(Ta_old(I))*den_old(I))*LaserTTM + &	! eqs. (13)
                   beta*LaserTTM*LaserTTM

!             Density of the carriers
              dn_dt_old(I) = alpha(Ta_old(I))*LaserTTM/hV + &
                   0.5d0*beta*LaserTTM*LaserTTM/hV - &                    ! eq. (17)
                   gamma*den_old(I)*den_old(I)*den_old(I) + &
                   teta(Egap_old(I),Te_old(I))*den_old(I) - &
                   (Ju_old(I)-Ju_old(I-1))/difz

              Den_new(I) = den_old(I) + dn_dt_old(I)*dift

!             Lattice temperature
              dT_l_dt_old(I) = ((AKl(I)*dTl(I)-AKl(I-1)*dTl(I-1))*dz2 + &          ! eq. (18)
                   C_e_h_old(i)*(Te_old(I)-Ta_old(I))/t_rel)/C_lat* Lcouple
                   
              T_l_new(I) = Ta_old(I) + dT_l_dt_old(I)*dift


           enddo Int_Loop

!          Perform calculations for the last node

!          Source term
	   if (LaserTTM .ge. almost_zero) then
		   Istar   = IL(N-1)    + 0.25d0*difz*dI_dz(IL(N-1),ATl(N-1),ADn(N-1))
		   denn    = den_old(N) - 0.25d0*dDn(N-1)
		   temm    = Ta_old(N) - 0.25d0*dTl(N-1)
		   IL(N) = IL(N-1)    +  0.5d0*difz*dI_dz(Istar,temm,denn)
		   LaserTTM   = 0.5d0*(IL(N-1) + IL(N))
	   else
		   Istar   = 0.d0
		   denn    = 0.d0
		   temm    = 0.d0
		   IL(N) = 0.d0
		   LaserTTM   = 0.d0
	   endif
	   
!          Auxiliary definitions
	C_e_h_old(N) = 1.5d0*kb*den_old(N)*(H3212e(N) + H3212h(N) &                       ! eq. (10)
			   + Te_old(N)*dtheta_e_dTe_old(N)*(1.d0-H3212e(N)/H12_12e(N)) &              ! H12_12 = 1/( H_1212 )
			   + Te_old(N)*dtheta_h_dTe_old(N)*(1.d0-H3212h(N)/H12_12h(N)) 	)  ! + Den * dE_gap/dT_e

           C_lat = C_l(Ta_old(N))                                 ! lattice heat conductivity
           t_rel = tau_e(den_old(N))                             ! relaxation time
           f_las_old(N) = (alpha(temm) + theta(temm)*denn)*LaserTTM + beta*LaserTTM*LaserTTM	! eqs. (13)

!          Density of the carriers
           dn_dt_old(N) = alpha(temm)*LaserTTM/hV + &
                0.5d0*beta*LaserTTM*LaserTTM/hV - &
                gamma*den_old(N)*den_old(N)*den_old(N) + &		! eq. (31)
                teta(Egap_old(N),Te_old(N))*den_old(N) + &
                2.d0*Ju_old(N-1)/difz

           Den_new(N) = den_old(N) + dn_dt_old(N)*dift

!          Lattice temperature
           dT_l_dt_old(N) = (-2.d0*AKl(N-1)*dTl(N-1)*dz2 + &		! eq. (32)
                C_e_h_old(N)*(Te_old(N)-Ta_old(N))/t_rel)/C_lat*Lcouple
                
           T_l_new(N) = Ta_old(N) + dT_l_dt_old(N)*dift




do i=1, N
	T_e_new(i)=Te_old(i)
	!~ Den_new(i)=Den_old(i)
	dn_dt_new(i)=dn_dt_old(i)
	dT_l_dt_new(i)=dT_l_dt_old(i)
	f_las_new(i) = f_las_old(i)
	C_e_h_new(i)=C_e_h_old(i)
	!~ T_l_new(i)=Ta_old(i)
	Ju_new(i)=Ju_old(i)
	H3212old(i) = H3212e(i) + H3212h(i)
	H3212e_old(I) = H3212e(I)
	H3212h_old(I) = H3212h(I)
	H12_12e_old(i) = H12_12e(i)
	H12_12h_old(i) = H12_12h(i)
	dtheta_e_dn_new(I)=dtheta_e_dn_old(I)
	dtheta_h_dn_new(I)=dtheta_h_dn_old(I)
	dtheta_e_dTe_new(I)=dtheta_e_dTe_old(I)
	dtheta_h_dTe_new(I)=dtheta_h_dTe_old(I)
enddo


corr_step = 0
! ************* Crank-Nicolson *****************************

!   Boundary conditions: left

222 do i=1, N
	Egap_new(i)     = EnergyGap(T_l_new(i),Den_new(i))
	Egap_n_new(i)   = EnergyGap_n(T_l_new(i),Den_new(i))
	Egap_T_l_new(i) = EnergyGap_T_l(T_l_new(i),Den_new(i))
       enddo


i=1

	f_old(i) =    F_las_old(i) - 2.d0*W(i)/difz - C_e_h_old(i)/tau_e(den_old(i))*(Te_old(i)-Ta_old(i))*Lcouple &	! eq. (30)
			- dn_dt_old(i)*( Egap_old(i)+1.5d0*kb*Te_old(i)*H3212old(i) ) - Den_old(I)*(Egap_n(I)*dn_dt_old(I) + Egap_T_l(I)*dT_l_dt_old(i))  &
			- 1.5d0*kb*den_old(i)*Te_old(i)*dn_dt_old(i) * &
														( (1.d0-H3212e_old(I)/H12_12e_old(I)) * dtheta_e_dn_old(I) &
														+(1.d0-H3212h_old(I)/H12_12h_old(I)) * dtheta_h_dn_old(I) )

	k_carriers_right =  k_c(T_e_new(1))+k_c(T_e_new(2))
	F_cn = psi*dift/C_e_h_new(1)

	a_cn(1) = 0.d0	! eq. (38)
	
	c_cn(i) = - F_cn/difz * 2.d0 * (  - kb*AH10(i)*Ju_new(i) + k_carriers_right/difz  )	! eq. (44)

	b_cn(i) = 1.d0 - F_cn *  ( - 2.d0*kb/difz * AH10(i)*Ju_new(i) - 2.d0*k_carriers_right/(difz*difz) &
						    - C_e_h_new(i)/tau_e(Den_new(i))*Lcouple &			! eq. (43)
						    - dn_dt_new(i)*1.5d0*kb*(H3212e(i)+H3212h(i))  &
		
		- 1.5d0*kb*den_new(i)*dn_dt_new(i) * &
														( (1.d0-H3212e(I)/H12_12e(I)) * dtheta_e_dn_new(I) &
														+(1.d0-H3212h(I)/H12_12h(I)) * dtheta_h_dn_new(I) )       )
    
	r_cn(i) = Te_old(i) + (1.d0-psi)*dift*f_old(i)/C_e_h_old(i) + F_cn *  (   f_las_new(i) - (Egap_new(i)+Egap_new(i+1)) * Ju_new(i) / difz &
			+ C_e_h_new(i) / tau_e( Den_new(i) ) * T_l_new(i) *Lcouple&		! eq. (45)
			- dn_dt_new(i) * Egap_new(i) &
			- Den_new(I)*( Egap_n_new(I)*dn_dt_new(I) + Egap_T_l_new(I)*dT_l_dt_new(i) )   )

!   Boundary conditions: right

i=N
	f_old(i) =  F_las_old(i) + 2.d0*W(N-1)/difz - C_e_h_old(i)/tau_e(den_old(i))*(Te_old(i)-Ta_old(i))*Lcouple &	! eq. (34)
			- dn_dt_old(i)*( Egap_old(i)+1.5d0*kb*Te_old(i)*H3212old(i) ) - Den_old(I)*(Egap_n(I)*dn_dt_old(I) + Egap_T_l(I)*dT_l_dt_old(i))   &
			- 1.5d0*kb*den_old(i)*Te_old(i)*dn_dt_old(i) * &
														( (1.d0-H3212e_old(I)/H12_12e_old(I)) * dtheta_e_dn_old(I) &
														+(1.d0-H3212h_old(I)/H12_12h_old(I)) * dtheta_h_dn_old(I) )
			

	k_carriers_left = k_c(T_e_new(N-1))+k_c(T_e_new(N))
	F_cn = psi*dift/C_e_h_new(N)

	a_cn(N) = - F_cn * ( 2.d0/difz*kb*AH10(N-1)*Ju_new(N-1) + 2.d0*k_carriers_left/(difz*difz) ) 		! eq. (46)

	b_cn(N) = 1.d0 - F_cn * ( 2.d0/difz * kb*AH10(N-1)*Ju_new(N-1) - 2.d0*k_carriers_left/(difz*difz) &	! eq. (47)
						    - C_e_h_new(i)/tau_e(Den_new(i))*Lcouple &
						    - dn_dt_new(i)*1.5d0*kb*(H3212e(i)+H3212h(i))  &
		
		- 1.5d0*kb*den_new(i)*dn_dt_new(i) * &
														( (1.d0-H3212e(I)/H12_12e(I)) * dtheta_e_dn_new(I) &
														+(1.d0-H3212h(I)/H12_12h(I)) * dtheta_h_dn_new(I) )       )
						
	c_cn(N) = 0.d0

	r_cn(N) = Te_old(i) + (1.d0-psi)*dift*f_old(i)/C_e_h_old(i) + F_cn *  (   f_las_new(i) + (Egap_new(i)+Egap_new(i-1)) * Ju_new(N-1) / difz &
					+ C_e_h_new(i) / tau_e( Den_new(i) ) * T_l_new(i) *Lcouple&		! eq. (48)
					- dn_dt_new(i) * Egap_new(i) &
					- Den_new(I)*( Egap_n_new(I)*dn_dt_new(I) + Egap_T_l_new(I)*dT_l_dt_new(i) )   )  


! Internal nodes:

        CNpre: do i=2, N-1

	f_old(i) =    F_las_old(i) - (W(i)-W(i-1))/difz - C_e_h_old(i)/tau_e(den_old(i))*(Te_old(i)-Ta_old(i))*Lcouple &	! eq. (20)
			- dn_dt_old(i)*( Egap_old(i)+1.5d0*kb*Te_old(i)*H3212old(i) ) - Den_old(I)*(Egap_n(I)*dn_dt_old(I) + Egap_T_l(I)*dT_l_dt_old(i))   &
			- 1.5d0*kb*den_old(i)*Te_old(i)*dn_dt_old(i) * &
														( (1.d0-H3212e_old(I)/H12_12e_old(I)) * dtheta_e_dn_old(I) &
														+(1.d0-H3212h_old(I)/H12_12h_old(I)) * dtheta_h_dn_old(I) )

	k_carriers_left = k_c(T_e_new(i))+k_c(T_e_new(i-1))
	k_carriers_right = k_c(T_e_new(i))+k_c(T_e_new(i+1))
	
	F_cn = psi * dift/C_e_h_new(i)
			
	a_cn(i) = - F_cn/difz *  (  kb*AH10(i-1)*Ju_new(i-1) + k_carriers_left/difz  )		! eq. (39)
	c_cn(i) = F_cn/difz *  (   kb*AH10(i)*Ju_new(i) - k_carriers_right/difz  )			! eq. (41)

	b_cn(i) = 1.d0 - F_cn *  ( - kb/difz * AH10(i)*Ju_new(i) + kb/difz*AH10(i-1)*Ju_new(i-1) &	! eq. (40)
						    - (k_carriers_left+k_carriers_right)/(difz*difz) &
						    - C_e_h_new(i)/tau_e(Den_new(i))*Lcouple &
						    - dn_dt_new(i)*1.5d0*kb*(H3212e(i)+H3212h(i))  &
		
		-  1.5d0*kb*den_new(i)*dn_dt_new(i) * &
														( (1.d0-H3212e(I)/H12_12e(I)) * dtheta_e_dn_new(I) &
														+(1.d0-H3212h(I)/H12_12h(I)) * dtheta_h_dn_new(I) )       )

	r_cn(i) = Te_old(i) + (1.d0-psi)*dift*f_old(i)/C_e_h_old(i) &					! eq. (42)
					+ F_cn *  (  f_las_new(i) - .5d0*(Egap_new(i)+Egap_new(i+1)) * Ju_new(i) / difz &
					+ .5d0*(Egap_new(i)+Egap_new(i-1)) * Ju_new(i-1) / difz &

		+ C_e_h_new(i) / tau_e( Den_new(i) ) * T_l_new(i) *Lcouple&
		- dn_dt_new(i) * Egap_new(i) &
		- Den_new(I)*( Egap_n_new(I)*dn_dt_new(I) + Egap_T_l_new(I)*dT_l_dt_new(i) )   )
		
		
        enddo CNpre

	T_e(1:N) = T_e_new(1:N)
	T_e_new(1:N) = Te_old(1:N)

        call tridag(a_cn,b_cn,c_cn,r_cn,T_e_new,N)

! *************/Crank-Nicolson predictor***********************************



! **************Preparation of new values for corrector*****************


!=============Fermi-integrals for corrector==============================================================
!          Perform calculations for Fermi-Dirak Integrals in all cells - prediction based on Te from CN      
                               Fermi_Loop_pre: do I=1, N

                                  Fhalfe = Den_new(I)*(kh1/(mee*T_e_new(I)))**1.5d0                      ! eq. (2)

                                  if (Fhalfe.LT.F12(0).OR.Fhalfe.GE.F12(Nini)) then
                                     print *,"chemical potential of electrons is out of bound (corrector)"				! the precalculated Fermi tables do not cover this value of the chemical potential
                                     print *,I,Den_new(I),T_e_new(I),Fhalfe
                                     stop
                                  endif

                                  IK = EKk(I)    ! restoring previous value of the index for certain cell
                                                   ! new value of IK shoul be near the old one

                                       do while (Fhalfe.GE.F12(IK+1))
                                          IK = IK + 1
                                       enddo
                                       do while (Fhalfe.LT.F12(IK))
                                          IK = IK - 1
                                       enddo
                                              
                                 EKk(I) = IK           ! saving new found value for the F12 index for every cell to use it at the next time step

                                  theta_e(I) = muni(IK) + &
                                       (Fhalfe-F12(IK))*(muni(IK+1)-muni(IK))/(F12(IK+1)-F12(IK))       ! linear approximation fot evaluation of eta(e)
				       
			if (ND_case) theta_e(i) = dlog(Fhalfe)

                               ! linear approximations for evaluation of F-functions
                                  FF12e  = F12(IK) + (F12(IK+1)-F12(IK))*(theta_e(I) - &
                                       muni(IK))/(muni(IK+1)-muni(IK))
                                  FF_12e = F_12(IK) + (F_12(IK+1)-F_12(IK))*(theta_e(I) - &
                                       muni(IK))/(muni(IK+1)-muni(IK))
                                  FF0e   = F0(IK) + (F0(IK+1)-F0(IK))*(theta_e(I) - &
                                       muni(IK))/(muni(IK+1)-muni(IK))
                                  FF1e   = F1(IK) + (F1(IK+1)-F1(IK))*(theta_e(I) - &
                                       muni(IK))/(muni(IK+1)-muni(IK))
                                  FF32e  = F32(IK) + (F32(IK+1)-F32(IK))*(theta_e(I) - &
                                       muni(IK))/(muni(IK+1)-muni(IK))
				       
			if (ND_case) then
				FF12e=Fhalfe
				FF_12e=Fhalfe
				FF0e=Fhalfe
				FF1e=Fhalfe
				FF32e=Fhalfe
			endif
				       
			F_12e_new(i) =FF_12e
                                  
			  H10e(I)    = FF1e/FF0e               ! definitions of H-functions
			  H12_12e(I) = FF12e/FF_12e
			  H012e(I)   = FF0e/FF12e
			  H3212e(I)  = FF32e/FF12e

			  Fhalfh = Den_new(I)*(kh1/(meh*T_e_new(I)))**1.5d0                      ! eq. (2)

			  if (Fhalfh.LT.F12(0).OR.Fhalfh.GT.F12(Nini)) then
			     print *,"chemical potential of holes is out of bound (corrector)"				! the precalculated Fermi tables do not cover this value of the chemical potential
			     print *,I,Den_new(I),T_e_new(I),F12(0),Fhalfh,F12(Nini)
			     stop
			  endif

			  JK = HK(I) 

		       do while (Fhalfh.GE.F12(JK+1))
			  JK = JK + 1
		       enddo
		       do while (Fhalfh.LT.F12(JK))
			  JK = JK - 1
		       enddo

			  theta_h(I) = muni(JK) + &
			       (Fhalfh-F12(JK))*(muni(JK+1)-muni(JK))/(F12(JK+1)-F12(JK))       ! linear approximation fot evaluation of eta(h)

			if (ND_case) theta_h(i)=dlog(Fhalfh)

                                  FF12h  = Fhalfh
                                ! linear approximations for evaluation of F-functions
                                 FF_12h = F_12(JK) + (F_12(JK+1)-F_12(JK))*(theta_h(I) - &
                                       muni(JK))/(muni(JK+1)-muni(JK))
                                  FF0h   = F0(JK) + (F0(JK+1)-F0(JK))*(theta_h(I) - &
                                       muni(JK))/(muni(JK+1)-muni(JK))
                                  FF1h   = F1(JK) + (F1(JK+1)-F1(JK))*(theta_h(I) - &
                                       muni(JK))/(muni(JK+1)-muni(JK))
                                  FF32h  = F32(JK) + (F32(JK+1)-F32(JK))*(theta_h(I) - &
                                       muni(JK))/(muni(JK+1)-muni(JK))

			if (ND_case) then
				FF_12h=Fhalfh
				FF0h=Fhalfh
				FF1h=Fhalfh
				FF32h=Fhalfh
			endif

			F_12h_new(i) =FF_12h

                                  HK(I) = JK

                                  H10h(I)    = FF1h/FF0h
                                  H12_12h(I) = FF12h/FF_12h
                                  H012h(I)   = FF0h/FF12h
                                  H3212h(I)  = FF32h/FF12h

                               enddo Fermi_Loop_pre
!          End of Fermi-Dirac Loop                                      
!============/Fermi-integrals for corrector==============================================================


! predicted heat capacity of e-h pairs
do i=1, N
	dtheta_e_dn_new(I) = .5d0 * ( 2.d0*pi*hp*hp/(mee*kb*T_e_new(i)) )**(1.5d0) / F_12e_new(I)		! eq. (25)
	dtheta_h_dn_new(I) = .5d0 * ( 2.d0*pi*hp*hp/(meh*kb*T_e_new(i)) )**(1.5d0) / F_12h_new(I)		! eq.( 25)
	dtheta_e_dTe_new(I) = -3.d0/dsqrt(2.d0)*den_new(i)*T_e_new(I)**(-2.5d0) * (pi*hp*hp/(mee*kb))**(1.5d0) / F_12e_new(I)	! eq. (26)
	dtheta_h_dTe_new(I) = -3.d0/dsqrt(2.d0)*den_new(i)*T_e_new(I)**(-2.5d0) * (pi*hp*hp/(meh*kb))**(1.5d0) / F_12h_new(I)	! eq. (26)
enddo


	do i=1, N
		C_e_h_new(i) = 1.5d0*kb*Den_new(I)*(H3212e(I) + H3212h(I) &                      ! eq. (10)
			   + T_e_new(i)*dtheta_e_dTe_new(i)*(1.d0-H3212e(i)/H12_12e(i)) &              ! H12_12 = 1/( H_1212 )
			   + T_e_new(i)*dtheta_h_dTe_new(i)*(1.d0-H3212h(i)/H12_12h(i)) 	)  ! + Den * dE_gap/dT_e
        enddo


! predicted values:
	do i=1, N-1
		AKl(I) = 0.5d0*(K_l(T_l_new(I)) + K_l(T_l_new(I+1)))
		dTl(I) = T_l_new(I+1) - T_l_new(I)
	enddo

 !dT_l_dt 
	C_lat = C_l(T_l_new(1))                                 ! lattice heat capacity
	t_rel = tau_e(Den_new(1))                             ! relaxation time

	dT_l_dt_new(1) = (2.d0*AKl(1)*dTl(1)*dz2 + C_e_h_new(1)*(T_e_new(1)-T_l_new(1))/t_rel*Lcouple)/C_lat		! eqs. (28) and (49)
	
	do i=2,N
		C_lat = C_l(T_l_new(i))                                    ! lattice heat capacity
		t_rel = tau_e(Den_new(i))                             ! relaxation time
		dT_l_dt_new(I) = ((AKl(I)*dTl(I)-AKl(I-1)*dTl(I-1))*dz2 + C_e_h_new(i)*(T_e_new(I)-T_l_new(I))/t_rel*Lcouple)/C_lat   ! eqs. (18) and (49)
	enddo
	
	C_lat = C_l(T_l_new(N))                                    ! lattice heat capacity
	t_rel = tau_e(Den_new(N))                             ! relaxation time
	dT_l_dt_new(N) = (-2.d0*AKl(N-1)*dTl(N-1)*dz2 + C_e_h_new(N)*(T_e_new(N)-T_l_new(N))/t_rel*Lcouple)/C_lat	! eqs. (32) and (49)




!          Perform calculations of J for cell-connecting points
           DWJ_loop_test: do I = 1, N-1

!             Auxiliary multipyers
              dTe(I) = T_e_new(I+1) - T_e_new(I)
              dTl(I) = T_l_new(I+1) - T_l_new(I)
              dDn(I) = Den_new(I+1) - Den_new(I)

              ATe(I) = 0.5d0*(T_e_new(I) + T_e_new(I+1))
              ATl(I) = 0.5d0*(T_l_new(I) + T_l_new(I+1))
              ADn(I) = 0.5d0*(Den_new(I) + Den_new(I+1))
              AKl(I) = 0.5d0*(K_l(T_l_new(I)) + K_l(T_l_new(I+1))) 

              AH012e(I)  = 0.5d0*(H012e(I) + H012e(I+1))
              AH012h(I)  = 0.5d0*(H012h(I) + H012h(I+1))
              AH10(I)    = 0.5d0*(H10e(I)  + H10e(I+1) + &
                   H10h(I) + H10h(I+1))
              AH12_12(I) = 0.5d0*(H12_12e(I) + H12_12e(I+1) + &
                   H12_12h(I) + H12_12h(I+1))

!             Diffusion coefficient
              D(I) = ATe(I)*kb*mue*muh*AH012e(I)*AH012h(I)*AH12_12(I)/ &		! eqs. (23) and (49)
                   (mue*AH012e(I)+muh*AH012h(I))/el

!             carrier pair current
              Ju_new(I) = -D(I)*(dDn(I) + &
                   ADn(I)*(Egap_new(I+1)-Egap_new(I))/ATe(I)/kb/AH12_12(I) + &		! eqs. (22) and (49)
                   ADn(I)*dTe(I)*(2.d0*AH10(I)/AH12_12(I)-1.5d0)/ATe(I))/difz
              
           end do DWJ_loop_test



!          Perform calculations for the 1st node

!          Source term
	   if (t_curr .le. tpulse*10.d0) then
		IL0   = Source_t(t_curr+0.5d0*dift,T_l_new(1),tpulse)
		Istar = IL0    + 0.25d0*difz*dI_dz(IL0,T_l_new(1),Den_new(1))
		denn  = Den_new(1) + 0.25d0*dDn(1)
		temm  = T_l_new(1) + 0.25d0*dTl(1)
		IL(1) = IL0    +  0.5d0*difz*dI_dz(Istar,temm,denn)           
		LaserTTM = 0.5d0*(IL0 + IL(1))
	   else 
		IL0=0.d0
		Istar = 0.d0
		denn  = 0.d0
		temm  = 0.d0
		IL(1) = 0.d0         
		LaserTTM = 0.d0
	   endif


           f_las_new(1) = (alpha(temm) + theta(temm)*denn)*LaserTTM + beta*LaserTTM*LaserTTM	! eqs. (13) and (49)


!          Density of the carriers
           dn_dt_new(1) = alpha(temm)*LaserTTM/hV + &
                0.5d0*beta*LaserTTM*LaserTTM/hV - &					! eqs. (27) and (49)
                gamma*Den_new(1)*Den_new(1)*Den_new(1) + &
                teta(Egap_new(1),T_e_new(1))*Den_new(1) - &
                2.d0*Ju_new(1)/difz

           Den_new(1) = Den_old(1) + ((1.d0-psi)*dn_dt_old(1)+psi*dn_dt_new(1))*dift		! eq. (50)

!          Lattice temperature
            T_l_new(1) = Ta_old(1) + ((1.d0-psi)*dT_l_dt_old(1)+psi*dT_l_dt_new(1))*dift	! eq. (51)
    



!          Perform calculations for T_e, T_l, and Den for internal cells
           Int_Loop_test: do I = 2,N

!             Source term
	   if (LaserTTM .ge. almost_zero) then
	      Istar = IL(I-1) + 0.5d0*difz*dI_dz(IL(I-1),ATl(I-1),ADn(I-1))
              IL(I) = IL(I-1) + difz*dI_dz(Istar,T_l_new(I),Den_new(I))
              LaserTTM = 0.5d0*(IL(I-1) + IL(I))
	   else
              Istar = 0.d0
              IL(I) = 0.d0
              LaserTTM = 0.d0
	   endif
                   
              f_las_new(I) = (alpha(T_l_new(I)) + theta(T_l_new(I))*Den_new(I))*LaserTTM + &	! eqs. (13) and (49)
                   beta*LaserTTM*LaserTTM


!             Density of the carriers

              dn_dt_new(I) = alpha(T_l_new(I))*LaserTTM/hV + &
                   0.5d0*beta*LaserTTM*LaserTTM/hV - &					! eqs. (17) and (49)
                   gamma*Den_new(I)*Den_new(I)*Den_new(I) + &
                   teta(Egap_new(I),T_e_new(I))*Den_new(I) - &
                   (Ju_new(I)-Ju_new(I-1))/difz

              Den_new(I) = Den_old(I) + ((1.d0-psi)*dn_dt_old(i)+psi*dn_dt_new(i))*dift		! eq. (50)

!             Lattice temperature
               T_l_new(I) = Ta_old(I) + ((1.d0-psi)*dT_l_dt_old(i)+psi*dT_l_dt_new(i))*dift	! eq. (51)



           enddo Int_Loop_test

!          Perform calculations for the last node

!          Source term
           Istar   = IL(N-1)    + 0.25d0*difz*dI_dz(IL(N-1),ATl(N-1),ADn(N-1))
           denn    = Den_new(N) - 0.25d0*dDn(N-1)
           temm    = T_l_new(N) - 0.25d0*dTl(N-1)
           IL(N) = IL(N-1)    +  0.5d0*difz*dI_dz(Istar,temm,denn)
           LaserTTM   = 0.5d0*(IL(N-1) + IL(N))


           f_las_new(N) = (alpha(temm) + theta(temm)*denn)*LaserTTM + &	! eqs. (13) and (49)
                beta*LaserTTM*LaserTTM


!          Density of the carriers
           dn_dt_new(N) = alpha(temm)*LaserTTM/hV + &
                0.5d0*beta*LaserTTM*LaserTTM/hV - &
                gamma*Den_new(N)*Den_new(N)*Den_new(N) + &				! eqs. (31) and (49)
                teta(Egap_new(N),T_e_new(N))*Den_new(N) + &
                2.d0*Ju_new(N-1)/difz

          Den_new(N) = Den_old(N) + ((1.d0-psi)*dn_dt_old(N)+psi*dn_dt_new(N))*dift		! eq. (50)


!          Lattice temperature
            T_l_new(N) = Ta_old(N) + ((1.d0-psi)*dT_l_dt_old(N)+psi*dT_l_dt_new(N))*dift	! eq. (51)

! **************/Preparation of new values for corrector*****************

summ=0.d0
do i=1, N
	summ = summ + dabs(T_e_new(i)-T_e(i))
enddo 
			i=1


corr_step = corr_step+1

			if (corr_step.ge.5000000) then
					write (*,*) t_curr, Te_old(1), T_e(1),  T_e_new(1), summ, corr_step
					stop "Could not reach the neccessary precision, 5 000 000 corrections attempted"
			endif

! **************Corrector*****************
if ( summ.ge.corrector_precision ) goto 222
! **************/Corrector*****************


!          Perform calculations for Fermi-Dirak Integrals in all cells        
           Fermi_Loop: do I=1, N

              Fhalfe = Den_new(I)*(kh1/(mee*T_e_new(I)))**1.5d0                      ! eq. (2)

              if (Fhalfe.LT.F12(0).OR.Fhalfe.GE.F12(Nini)) then
                 print *,"chemical potential of electrons is out of bound 2"
                 print *,I,Den_new(I),T_e_new(I),F12(0),Fhalfe,F12(nini)
                 stop
              endif

              IK = EKk(I)    ! restoring previous value of the index for certain cell
                               ! new value of IK shoul be near the old one

                   do while (Fhalfe.GE.F12(IK+1))
                      IK = IK + 1
                   enddo
                   do while (Fhalfe.LT.F12(IK))
                      IK = IK - 1
                   enddo
                           
			EKk(I) = IK           ! saving new found value for the F12 index for every cell to use it at the next time step

			theta_e(I) = muni(IK) + &
			(Fhalfe-F12(IK))*(muni(IK+1)-muni(IK))/(F12(IK+1)-F12(IK))       ! linear approximation fot evaluation of eta(e)
			
		if (ND_case) theta_e(i) = dlog(Fhalfe)

			! linear approximations for evaluation of F-functions
			FF12e  = F12(IK) + (F12(IK+1)-F12(IK))*(theta_e(I) - &
			muni(IK))/(muni(IK+1)-muni(IK))
			FF_12e = F_12(IK) + (F_12(IK+1)-F_12(IK))*(theta_e(I) - &
			muni(IK))/(muni(IK+1)-muni(IK))
			FF0e   = F0(IK) + (F0(IK+1)-F0(IK))*(theta_e(I) - &
			muni(IK))/(muni(IK+1)-muni(IK))
			FF1e   = F1(IK) + (F1(IK+1)-F1(IK))*(theta_e(I) - &
			muni(IK))/(muni(IK+1)-muni(IK))
			FF32e  = F32(IK) + (F32(IK+1)-F32(IK))*(theta_e(I) - &
			muni(IK))/(muni(IK+1)-muni(IK))
			
			if (ND_case) then
				FF12e=Fhalfe
				FF_12e=Fhalfe
				FF0e=Fhalfe
				FF1e=Fhalfe
				FF32e=Fhalfe
			endif
			
		F_12e(I)=ff_12e
			
			H10e(I)    = FF1e/FF0e               ! definitions of H-functions
			H12_12e(I) = FF12e/FF_12e
			H012e(I)   = FF0e/FF12e
			H3212e(I)  = FF32e/FF12e

			Fhalfh = Den_new(I)*(kh1/(meh*T_e_new(I)))**1.5d0                      ! eq. (2)


              if (Fhalfh.LT.F12(0).OR.Fhalfh.GT.F12(Nini)) then
                 print *,"chemical potential of holes is out of bound"
                 print *,I,Den_new(I),T_e_new(I),Fhalfh,F12(0),F12(Nini)
                 stop
              endif

              JK = HK(I) 

                   do while (Fhalfh.GE.F12(JK+1))
                      JK = JK + 1
                   enddo
                   do while (Fhalfh.LT.F12(JK))
                      JK = JK - 1
                   enddo

		theta_h(I) = muni(JK) + &
		(Fhalfh-F12(JK))*(muni(JK+1)-muni(JK))/(F12(JK+1)-F12(JK))       ! linear approximation fot evaluation of eta(h)
		
			if (ND_case) theta_h(i) = dlog(Fhalfh)

		FF12h  = Fhalfh
		! linear approximations for evaluation of F-functions
		FF_12h = F_12(JK) + (F_12(JK+1)-F_12(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		FF0h   = F0(JK) + (F0(JK+1)-F0(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		FF1h   = F1(JK) + (F1(JK+1)-F1(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))
		FF32h  = F32(JK) + (F32(JK+1)-F32(JK))*(theta_h(I) - &
		muni(JK))/(muni(JK+1)-muni(JK))

			if (ND_case) then
				FF_12h=Fhalfh
				FF0h=Fhalfh
				FF1h=Fhalfh   
				FF32h=Fhalfh
			endif

	F_12h(I)=ff_12h
	
		HK(I) = JK

		H10h(I)    = FF1h/FF0h
		H12_12h(I) = FF12h/FF_12h
		H012h(I)   = FF0h/FF12h
		H3212h(I)  = FF32h/FF12h

		if (theta_e(I).LE.eta_emin) eta_emin = theta_e(I)             ! serching minimum and maximum of reduced Fermi level fot printing it later
		if (theta_e(I).GE.eta_emax) eta_emax = theta_e(I)
		if (theta_h(I).LE.eta_hmin) eta_hmin = theta_h(I)
		if (theta_h(I).GE.eta_hmax) eta_hmax = theta_h(I)

           enddo Fermi_Loop
!          End of Fermi-Dirack Loop                                      

! **************/Preparation of new values for energy calculation at the new step*****************
	


   ! ---------------Energy control---------------------------------------------------------
	F_total   = F_total   + 0.5d0*f_las_old(1)*dtz                                      ! total fluence
	E_total_l = E_total_l + 0.5d0*dT_l_dt_old(1)*dtz*C_l(Ta_old(1))
	E_total_e = E_total_e + 0.5d0*den_old(1)*(Egap_old(1) + &		! eq. (11)
	1.5d0*kb*Te_old(1)*H3212old(1))*difz 

	e_energ_den(1) =  den_old(1)*(Egap_old(1) + 1.5d0*kb*Te_old(1)*H3212old(1))
	l_energ_den(1) = l_energ_den(1) + dT_l_dt_old(1)*dift*C_l(Ta_old(1))
	tot_energ_den(1) = e_energ_den(1) + l_energ_den(1)
	
	E_add = E_add + 0.5d0*C_e_h_old(1)*(Te_old(1)-Ta_old(1))/tau_e(den_old(1))*dtz*Lcouple
	Etrans(1)=Etrans(1)-0.5d0*C_e_h_old(1)*(Te_old(1)-Ta_old(1))/tau_e(den_old(1)) * Lcouple *dift 
	E_nonthermal(1) = den_old(1)*Egap_old(1)

do i=2, N-1
	F_total   = F_total   + f_las_old(I)*dtz                                      ! total fluence
	E_total_l = E_total_l + dT_l_dt_old(I)*dtz*C_l(Ta_old(i))
	E_total_e = E_total_e + den_old(I)*(Egap_old(I) + &		! eq. (11)
		1.5d0*kb*Te_old(I)*H3212old(i))*difz

	e_energ_den(i) = den_old(I)*(Egap_old(I) + 1.5d0*kb*Te_old(I)*H3212old(i))		! eq. (11)
	l_energ_den(i) = l_energ_den(i) + dT_l_dt_old(i)*dift*C_l(Ta_old(i))
	tot_energ_den(i) = e_energ_den(i) + l_energ_den(i)

	E_add = E_add + C_e_h_old(I)*(Te_old(I)-Ta_old(I))/tau_e(den_old(i))*dtz*Lcouple
	Etrans(I)=Etrans(I)-C_e_h_old(I) *(Te_old(I)-Ta_old(I))/tau_e(den_old(i)) * Lcouple * dift
	E_nonthermal(i) = den_old(I)*Egap_old(I)
enddo

	F_total   = F_total   + 0.5d0*f_las_old(N)*dtz                                      ! total fluence
	E_total_l = E_total_l + 0.5d0*dT_l_dt_old(N)*dtz*C_l(Ta_old(N))
	E_total_e = E_total_e + 0.5d0*den_old(N)*(Egap_old(N) + &		! eq. (11)
		1.5d0*kb*Te_old(N)*H3212old(N))*difz          ! (16)
	E_add = E_add + 0.5d0*C_e_h_old(N)*(Te_old(N)-Ta_old(N))/tau_e(den_old(N))*dtz*Lcouple
	Etrans(N)=Etrans(N)-0.5d0*C_e_h_old(N)*(Te_old(N)-Ta_old(N))/tau_e(den_old(N)) * Lcouple * dift

	e_energ_den(N) = den_old(N)*(Egap_old(N) + 1.5d0*kb*Te_old(N)*H3212old(N)) 		! eq. (11)
	l_energ_den(N) = l_energ_den(N) + dT_l_dt_old(N)*dift*C_l(Ta_old(N))
	tot_energ_den(N) = e_energ_den(N) + l_energ_den(N)

	E_nonthermal(N) = den_old(N)*Egap_old(N)


! searching for the melting energy. the result is  3859199787 for Tm=1688 K
!~ if (Ta_old(1).ge.1688.d0) then
	!~ print *, Ta_old(1), E_nonthermal(1), l_energ_den(1)
	!~ stop
!~ endif


!          End of calculations for the edge cells
	E_total = E_total_0 + E_total_e + E_total_l
	tot_energ_dens = E_total / L
	

! --------------/Energy control---------------------------------------------------------

do i=1, N
	den_old(1:N) = den_new(1:N)
	Ta_old(1:N) = T_l_new(1:N)
	Te_old(1:N) = T_e_new(1:N)
enddo

           t_curr = t_curr + dift                                                                  ! time step ended, the time advances by dt

           if (MOD(NstepTTM,q).EQ.0.OR.NstepTTM.EQ.1) then
	   
			if (MOD(NstepTTM,q*100).EQ.0.OR.NstepTTM.EQ.1) then
                                do j=1, N
				      write(1,111) t_curr*1.d+12,X(J)*1.d+06,Te_old(J), &            ! time [ps], coordinate [mkm], electron temperature,
						Ta_old(J),den_old(J)*1.d-26                             ! lattice temperature, carriers density []
                                enddo
                                write (1,*) ''
                        endif
              write(112,112) t_curr*1.d+12,E_total,E_total+E_base, &
                   E_total_e+E_add,E_total_l-E_add,E_add, &
                   Te_old(1)*1.d-03,Ta_old(1)*1.d-03,den_old(1)*1.d-26, &
                   Te_old(2)*1.d-03,Ta_old(2)*1.d-03,den_old(2)*1.d-26, &
                   Te_old(N)*1.d-03,Ta_old(N)*1.d-03,den_old(N)*1.d-26, &
                   F_total,F_total+E_base,IL0*dift, &
		   l_energ_den(1), e_energ_den(1), &
		   l_energ_den(2), e_energ_den(2), &
		   l_energ_den(3), e_energ_den(3), &
		   l_energ_den(N), e_energ_den(N)
		   
              write(5,115) t_curr*1.d+12,theta_e(1),theta_h(1),eta_emin, &        ! reduced Fermi levels of e and h
                   eta_emax,eta_hmin,eta_hmax
		   
           endif

              if (MOD(NstepTTM,q*10).EQ.0.OR.NstepTTM.EQ.1) then
           i=1
                      write (*,'(E11.4, 1x, F10.3, 1x, F10.3, 2x, F10.3, 1x, F10.3, 3x, E16.8, F19.9,i5)')  &
					t_curr,Te_old(1), Ta_old(1) &
					, l_energ_den(1)/energ_den_melting, e_energ_den(1)/energ_den_melting, den_old(1) &
					, E_total
		endif

        if (MOD(NstepTTM,q).EQ.0.OR.NstepTTM.EQ.1) then
	   
		   if (t_curr .gt.(t_exit)) exit Time_loop        ! checking if simulation time has ended
   
   
				call CPU_TIME(t2)
				!                      write (*,'(4(E22.12,1x))') t_curr,T_e_new(1),T_e_new(1), E_total*.5d0
				write (22,'(17(E22.12,1x))') t_curr*1.d12,Te_old(1)/1.d3, Ta_old(1)/1.d3, den_old(1)/1.d26 &
									, l_energ_den(1), e_energ_den(1), E_total &
									, tot_energ_den(1), tot_energ_dens  &
									,Te_old(2)/1.d3, Ta_old(2)/1.d3 &
									, den_old(2)/1.d26 &
									, E_nonthermal(1)/energ_den_melting &
									, F_total,E_total_e, E_total_l, Egap_old(1)

				t1 = t2
              
           endif


enddo Time_loop


!============================================================



	open(UNIT = 3,FILE = 'prof.out')
        
        write(3,113)(X(I)*1.d+06,T_e_new(I)*1.d-03, &
             Ta_old(I)*1.d-03,den_old(I)*1.d-26,I=1, N)

109     FORMAT(3(E22.12,1x))
110     FORMAT(E22.12,1x,E22.12,1x,I6)
111     FORMAT(5(1x,E22.12))
112     FORMAT(26(1x,E22.12))
113     FORMAT(4(1x,E22.12))
115     FORMAT(7(1x,E22.12))

        close(UNIT = 1)
        close(UNIT = 2)
        close(UNIT = 3)
        close(UNIT = 4)
        close(UNIT = 5)
        close(UNIT = 22)

        call CPU_TIME(timef)

        print *,"# Total CPU time: ",timef-timen,"seconds"

      end Program nTTM

!***************************************************************************

	real(8) Function Source_t(time,tem,tpulse)
	implicit none 
	real(8) time,tem,tpulse
	real(8) omega
	real(8) Refl
	EXTERNAL Refl
	real(8), parameter :: pi     = 3.1415926535897932384626433832795d0  !
	!real(8), parameter :: PHI    = 2600.d0             ! [J m-2]            laser energy per unit of area (incident fluence)
	!real(8), parameter :: PHI = 50.d0 
	!real(8), parameter :: PHI = 150.d0 					![J m-2]
	real(8), parameter :: PHI = 1500.d0					! [J m-2]
	real(8) :: tend     ! [s]
		omega  = 4.d0*dlog(2.d0)                          ! so integration gives Fluence
		tend   = 3.d0*tpulse    ! [s]
		Source_t = dsqrt(omega/pi) * (1.d0-Refl(tem))*PHI*dexp(-omega*((time-tend)/tpulse)**2.d0)/tpulse    ! eq. (14)
	end Function Source_t



	real(8) Function dI_dz(St,tem,dens)    !  dI/dz
	implicit none 
	real(8) St,tem,dens                              ! intensity, temperature, density
	real(8) alpha,theta
	EXTERNAL alpha,theta
	real(8), parameter :: beta = 15.d-11        ! [m W-1]
		dI_dz = -(alpha(tem) + theta(tem)*dens)*St - beta*St*St       ! eq. (15)
	end Function dI_dz


	real(8) Function C_l(Tc)           ! [J m-3]                      lattice specific heat capacity
	implicit none 
	real(8) Tc
		C_l = 1.978d+06 + 3.54d+02*Tc - 3.68d+06/(Tc*Tc)                     ! table 1
	end Function C_l


	real(8) Function K_l(Tc)           ! [W m-1 K-1]                   lattice thermal conductivity  
	implicit none 
	real(8) Tc
		K_l = 1.585d+05*Tc**(-1.23d0)                     ! table 1
	end Function K_l


	real(8) Function K_c(Tc)           ! [W m-1 K-1]                  carriers thermal conductivity  
	implicit none 
	real(8) Tc
	real(8), parameter :: en  = 1.602176487d-09   ! [A J eV-1 M-1]
		K_c = -3.47d+08 + 4.45d+06*Tc                                        ! table 1  
		K_c = K_c*en                    !    changing of units
	end Function K_c


	real(8) Function alpha(Tc)         ! [m-1]       one-photon absorption coefficient
	implicit none 
	real(8) Tc
		!alpha = 1.34d5*dexp(Tc/427.d0)						! table 1
		alpha = 5.02d5*dexp(Tc/427.d0)						!value taken from Chen and Beraun
	end Function alpha


	real(8) Function Refl(Tc)                ! surface reflectivity
	implicit none 
	real(8) Tc
		Refl = 0.329d0 + 5.d-05*(Tc - 300.d0)					! table 1   
	end Function Refl


	real(8) Function teta(Ezone,Tc)    ! [s-1]              impact ionization coefficient
	implicit none 
	real(8) Ezone,Tc
	real(8), parameter :: kb = 1.3806504d-23 ! [J K-1]
		teta = 3.6d+10*dexp(-1.5d0*Ezone/(kb*Tc))				 ! table 1                
	end Function teta


	real(8) Function theta(Tc)           ! [m2]           free-carrier absorption cross-section
	real(8) Tc
	real(8), parameter :: Trm = 300.d0       ! room temperature
		!theta = 2.91d-22*Tc/Trm								! table 1
		theta = 5.1d-22*Tc/Trm					!value taken from Chen and Beraun
	end Function theta 


	real(8) Function tau_e(densit)     ! [s]                relaxation time of electrons     
	implicit none 
	real(8) densit
		!tau_e = 0.5d-12*(1.d0 + densit/(2.d+27))    ! [density] = 1/m^-3			! table 1
		tau_e = 2.4d-13*(1.d0 + densit/(2.d+27))    ! [density] = 1/m^-3			! table 1
	end Function tau_e


	real(8) Function EnergyGap(Tc,densit)   ! [J]
	implicit none 
	real(8) Tc,densit                               ! temperature, K ; density, m^-3
	real(8), parameter :: en = 1.602176487d-19                   ! changing eV to J
		EnergyGap = 1.170d0-4.73d-04*Tc*Tc/(Tc + 636.d0) - &       				! table 1
				  1.5d-10*(densit**0.33333333333333333d0)            ! cm^-3 --> m^-3
		if (EnergyGap .lt. 0.d0) then
			EnergyGap = 0.d0 ! you also need to make two other functions 0
		else
			EnergyGap = EnergyGap*en                   ! changing eV to J
		endif
	end Function EnergyGap


	real(8) Function EnergyGap_n(Tc,densit)             ! [J cm3]   
	implicit none 
	real(8) Tc,densit                               ! temperature, K ; density, m^-3
	real(8) EnergyGap
	EXTERNAL EnergyGap
	real(8), parameter :: en = 1.602176487d-19

		EnergyGap_n = -0.5d-10/densit**(0.6666666666666666666d0)     !   dEg/dn
		if (EnergyGap(Tc,densit) .lt. 1.d-50) then
			EnergyGap_n = 0.d0
		else
			EnergyGap_n = EnergyGap_n*en
		endif
	end Function EnergyGap_n


	real(8) Function EnergyGap_T_l(Tc,densit)               ! [J K-1]         
	implicit none 
	real(8) Tc,densit                               ! temperature, K ; density, m^-3
	real(8) EnergyGap
	EXTERNAL EnergyGap
	real(8), parameter :: en = 1.602176487d-19
		EnergyGap_T_l = -4.73d-04*(Tc*Tc + 1272.d0*Tc)/ &     !   dEg/dTl
		((Tc + 636.d0)*(Tc + 636.d0))
		if (EnergyGap(Tc,densit) .lt. 1.d-50) then
			EnergyGap_T_l = 0.d0
		else
			EnergyGap_T_l = EnergyGap_T_l*en
		endif
	end Function EnergyGap_T_l




!!!!!!!!!!!!!!!!!!!!!!!!!TRIDAG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Press, W. H., Teukolsky, S. A., Flannery, B. P., & Vetterling, W. T. (1992). Numerical recipes in Fortran 77: volume 1, volume 1 of Fortran numerical recipes: the art of scientific computing. Cambridge university press.
subroutine tridag(a,b,c,r,T,n)
!Solves tridiagonal matrix equation for T(1:n)
intent(in) :: a, b, c, r, n
intent(inout) :: T
integer :: n,NMAX
real(8) :: a(n),b(n),c(n),r(n),T(n)
parameter (NMAX=500)

integer :: j
real(8) :: ki,temp(NMAX)
	if (b(1).lt.1.d-100) stop "tridag: equations must be rewritted"
	ki=b(1)
	T(1)=r(1)/ki
	do j=2,n ! Decomposition and forward substitution
		temp(j)=c(j-1)/ki
		ki=b(j)-a(j)*temp(j)
		if(ki.lt.1.d-100) stop "tridag: equations must be rewritted"
		T(j)=(r(j)-a(j)*T(j-1))/ki
	enddo
	do j=n-1,1,-1 ! Backsubstitution
	T(j)=T(j)-temp(j+1)*T(j+1)
enddo

end subroutine tridag

 
