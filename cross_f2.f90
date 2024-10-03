Program cross_f2
	implicit none


	real (kind=8) :: s, t, t1, t2, mf, mp, Ampl1_rho, Ampl2_rho, Ampl3_rho
	real (kind=8) :: Ampl4_rho, Ampl_rho, Ampl, Ampl1_w, Ampl2_w, Ampl3_w
	real (kind=8) :: ga_rho, gb_rho, gt_rho, gv_rho, gs_rho, alpha_rho, alphal_rho
	real (kind=8) :: ga_w, gb_w, gt_w, gv_w, gs_w, alpha_w, alphal_w, Ampl_w 
	real (kind=8) :: cross, conv, E, Ampl4_w, pi, gamma_rho, mrho, tm
	real (kind=8) :: gamma_w, mw, eu, ed
	integer :: i
	complex :: A_rho, A_w, Aconj_rho, Aconj_w, B_rho, Bconj_rho, B_w, Bconj_w, o
!
!	Arquivo de resultados.
!
	open (unit = 50, file = 'cross_t.dat')
!	
!	Parametros gerais.
!
	cross = 0.0
	o = (0,1)
	mp = 0.938
	pi = acos(-1.)
    alphal_w = 0.9
	alphal_rho = 0.8
	mrho = 0.775
	mw = 0.782
!
!	Parametros para f_2(1270)
!
	mf = 1.2754
	E = 3.8
	gv_rho = 3.4
	gt_w = 0.0
	gv_w = 15.
	gt_rho = 11.
	gamma_rho = 0.57e-3
	gamma_w = 0.57e-3
	eu = 2./3.
	ed = -1./3.
!
!
!
    conv=0.3894e6 
!	
	s = mp**2.+2.*E*mp
!	
	t1= 1./(2.*s) * (-1.*(mp**2. - s)**2. + mf**2. *(mp**2. + s) &
	+(mp**2.-s)*sqrt((mp**2. - s)**2. - 2.*mf**2.*(mp**2.+s)+mf**4.))
	t2= 1./(2.*s) * (-1.*(mp**2. - s)**2. + mf**2. *(mp**2. + s) &
	-(mp**2.-s)*sqrt((mp**2. - s)**2. - 2.*mf**2.*(mp**2.+s)+mf**4.))
!	
!	Acoplamento para f_2(1270)
!
	gs_rho = sqrt(32.*pi*gamma_rho*mf**3./(mf**2.-mrho**2.)**3.)
	gs_w = sqrt(32.*pi*gamma_w*mf**3./(mf**2.-mw**2.)**3.)
!
!
!
	ga_rho = gs_rho*(gv_rho+2.*mp*gt_rho)
	gb_rho = 2.*gs_rho*gt_rho
	ga_w= gs_w*(gv_w+2.*mp*gt_w)
	gb_w = 2.*gs_w*gt_w
!
	tm = 0.05
	do 10, i=1,2000
!
	t = -tm
	alpha_rho = 0.55 + alphal_rho * t
	alpha_w = 0.44 + alphal_w * t
!	
	A_rho = ga_rho * s**(alpha_rho -1.)*pi * alphal_rho /(sin(pi*alpha_rho)) &
	*(1.-exp(-o*pi*alpha_rho))/(2.*gamma(alpha_rho))
	Aconj_rho = ga_rho * s**(alpha_rho -1.)*pi * alphal_rho /(sin(pi*alpha_rho)) &
	*(1.-exp(o*pi*alpha_rho))/(2.*gamma(alpha_rho))
	B_rho = -gb_rho/ga_rho * A_rho
	Bconj_rho = -gb_rho/ga_rho * Aconj_rho
	A_w = ga_w * s**(alpha_w -1.)*pi * alphal_w /(sin(pi*alpha_w)) &
	*(1.-exp(-o*pi*alpha_w))/(2.*gamma(alpha_w))
	Aconj_w = ga_w * s**(alpha_w -1.)*pi * alphal_w /(sin(pi*alpha_w)) &
	*(1.-exp(o*pi*alpha_w))/(2.*gamma(alpha_w))
	B_w = -gb_w/ga_w * A_w
	Bconj_w = -gb_w/ga_w * Aconj_w
!
	Ampl1_rho = -1./2.* (A_rho * Aconj_rho)*(s*(t-t1)*(t-t2) &
	+1./2.* t*(t**2. - 2. * (mf**2. + s)* t + mf**4.))
	Ampl2_rho = -1./2.* (A_rho * Bconj_rho)*mp*s*(t-t1)*(t-t2)
	Ampl3_rho = -1./2.* (Aconj_rho * B_rho)*mp*s*(t-t1)*(t-t2)
	Ampl4_rho = -1./8.*(B_rho * Bconj_rho)*s*(4.*mp**2.-t)*(t-t1)*(t-t2)
	Ampl_rho = Ampl1_rho + Ampl2_rho + Ampl3_rho + Ampl4_rho
	Ampl1_w = -1./2.* (A_w * Aconj_w)*(s*(t-t1)*(t-t2) &
	+1./2.* t*(t**2. - 2. * (mf**2. + s)* t + mf**4.))
	Ampl2_w = -1./2.* (A_w * Bconj_w)*mp*s*(t-t1)*(t-t2)
	Ampl3_w = -1./2.* (Aconj_w * B_w)*mp*s*(t-t1)*(t-t2)
	Ampl4_w = -1./8.*(B_w * Bconj_w)*s*(4.*mp**2.-t)*(t-t1)*(t-t2)
	Ampl_w = Ampl1_w + Ampl2_w + Ampl3_w + Ampl4_w
	Ampl = Ampl_rho+Ampl_w
!
	cross = conv * Ampl/(16.*pi*(s-mp**2.)**2.)
!	
	write(50,*) tm, cross
!	write(*,*) sqrt(s)
!	
	tm = tm + 0.001
 10 continue	
!
	End program cross_f2