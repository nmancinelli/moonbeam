# RTCOEF calculates P/SV reflection/transmission coefficients
#   for an interface between two solid layers, based on the 
#   equations on p. 149-150 of Aki and Richards.  This version 
#   is modified from an older routine provided by Tom Sereno.
#
#  Inputs:    vp1     =  P-wave velocity of layer 1 (top layer)
#  (real)     vs1     =  S-wave velocity of layer 1
#             den1    =  density of layer 1
#             vp2     =  P-wave velocity of layer 2 (bottom layer)
#             vs2     =  S-wave velocity of layer 2
#             den2    =  density of layer 2
#             hslow   =  horizontal slowness (ray parameter)
#  Returns:   rt(1)   =  down P to P up     (refl)
#  (complex)  rt(2)   =  down P to S up     (refl)
#             rt(3)   =  down P to P down   (tran)
#             rt(4)   =  down P to S down   (tran)
#             rt(5)   =  down S to P up     (refl)
#             rt(6)   =  down S to S up     (refl)
#             rt(7)   =  down S to P down   (tran)
#             rt(8)   =  down S to S down   (tran)
#             rt(9)   =    up P to P up     (tran)
#             rt(10)  =    up P to S up     (tran)
#             rt(11)  =    up P to P down   (refl)
#             rt(12)  =    up P to S down   (refl)
#             rt(13)  =    up S to P up     (tran)
#             rt(14)  =    up S to S up     (tran)
#             rt(15)  =    up S to P down   (refl)
#             rt(16)  =    up S to S down   (refl)
#
# NOTE:  All input variables are real.  
#        All output variables are complex!
#        Coefficients are not energy normalized.
#
def rtcoef(vp1,vs1,den1,vp2,vs2,den2,hslow):
    #implicit complex (a-h,o-z)
    #complex rt(16)
    #real vp1,vs1,den1,vp2,vs2,den2,hslow
    from numpy import sqrt

    alpha1=vp1+0j
    beta1=vs1+0j
    rho1=den1+0j
    alpha2=vp2+0j
    beta2=vs2+0j
    rho2=den2+0j
    p=hslow+0j

    cone=1+0j
    ctwo=2+0j
    
    rt={}

    term1=(cone-ctwo*beta1**2*p**2)
    term2=(cone-ctwo*beta2**2*p**2)
    a=rho2*term2-rho1*term1
    b=rho2*term2+ctwo*rho1*beta1**2*p**2
    c=rho1*term1+ctwo*rho2*beta2**2*p**2
    d=ctwo*(rho2*beta2**2-rho1*beta1**2)

    # compute signs and cosines, allowing for complex incidence angles
    si1=alpha1*p 
    si2=alpha2*p          
    sj1=beta1*p
    sj2=beta2*p          
    ci1=sqrt(cone-si1**2)
    ci2=sqrt(cone-si2**2)
    cj1=sqrt(cone-sj1**2)
    cj2=sqrt(cone-sj2**2)         

    E=b*ci1/alpha1+c*ci2/alpha2
    F=b*cj1/beta1+c*cj2/beta2
    G=a-d*ci1*cj2/(alpha1*beta2)
    H=a-d*ci2*cj1/(alpha2*beta1)
    DEN=E*F+G*H*p**2
      
    trm1=b*ci1/alpha1-c*ci2/alpha2          
    trm2=a+d*ci1*cj2/(alpha1*beta2)
    rt[1]=(trm1*F-trm2*H*p**2)/DEN         #refl down P to P up

    trm1=a*b+c*d*ci2*cj2/(alpha2*beta2)       
    rt[2]=(-ctwo*ci1*trm1*p)/(beta1*DEN)   #refl down P to S up

    rt[3]=ctwo*rho1*ci1*F/(alpha2*DEN)     #trans down P to P down

    rt[4]=ctwo*rho1*ci1*H*p/(beta2*DEN)    #trans down P to S down

    trm1=a*b+c*d*ci2*cj2/(alpha2*beta2)       
    rt[5]=(-ctwo*cj1*trm1*p)/(alpha1*DEN)  #refl down S to P up

    trm1=b*cj1/beta1-c*cj2/beta2               
    trm2=a+d*ci2*cj1/(alpha2*beta1)
    rt[6]=-(trm1*E-trm2*G*p**2)/DEN        #refl down S to S up

    rt[7]=-ctwo*rho1*cj1*G*p/(alpha2*DEN)  #trans down S to P down 

    rt[8]=ctwo*rho1*cj1*E/(beta2*DEN)      #trans down S to S down

    rt[9]=ctwo*rho2*ci2*F/(alpha1*DEN)     #trans up P to P up

    rt[10]=-ctwo*rho2*ci2*G*p/(beta1*DEN)  #trans up P to S up

    trm1=b*ci1/alpha1-c*ci2/alpha2          
    trm2=a+d*ci2*cj1/(alpha2*beta1)
    rt[11]=-(trm1*F+trm2*G*p**2)/DEN       #refl up P to P down

    trm1=a*c+b*d*ci1*cj1/(alpha1*beta1)       
    rt[12]=(ctwo*ci2*trm1*p)/(beta2*DEN)   #refl up P to S down

    rt[13]=ctwo*rho2*cj2*H*p/(alpha1*DEN)  #trans up S to P up

    rt[14]=ctwo*rho2*cj2*E/(beta1*DEN)     #trans up S to S up

    trm1=a*c+b*d*ci1*cj1/(alpha1*beta1)       
    rt[15]=(ctwo*cj2*trm1*p)/(alpha2*DEN)  #refl up S to P down

    trm1=b*cj1/beta1-c*cj2/beta2               
    trm2=a+d*ci1*cj2/(alpha1*beta2)
    rt[16]=(trm1*E+trm2*H*p**2)/DEN        #refl up S to S down

    return rt

def rtcoefSH(vs1,den1,vs2,den2,hslow):
    from numpy import arcsin, cos
    j1=arcsin(hslow*vs1)
    j2=arcsin(hslow*vs2)
    a=vs1*den1*cos(j1) - vs2*den2*cos(j2)
    b=vs1*den1*cos(j1) + vs2*den2*cos(j2)

    c=2*den1*vs1*cos(j1)

    return a/b, c/b

def rcoef_free(vp,vs,hslow):
    from numpy import arcsin, cos
    i = arcsin(hslow*vp)
    j = arcsin(hslow*vs)
    p=hslow
    alpha=vp
    beta=vs

    r={}

    b=(1./beta**2-2*p**2)**2 + 4.*p**2*cos(i)*cos(j)/alpha/beta
    a=4.*beta/alpha*p*cos(j)/beta*(1./beta**2 - 2*p**2)
    SupPdn = a/b

    b=(1./beta**2-2*p**2)**2 + 4.*p**2*cos(i)*cos(j)/alpha/beta
    a=(1./beta**2-2*p**2)**2 - 4.*p**2*cos(i)*cos(j)/alpha/beta
    SupSdn = a/b

    b= (1./beta**2-2*p**2)**2 + 4.*p**2*cos(i)*cos(j)/alpha/beta
    a=-(1./beta**2-2*p**2)**2 + 4.*p**2*cos(i)*cos(j)/alpha/beta
    PupPdn = a/b

    b= (1./beta**2-2*p**2)**2 + 4.*p**2*cos(i)*cos(j)/alpha/beta
    a=4.*alpha/beta*p*cos(i)/alpha*(1./beta**2 - 2*p**2)
    PupSdn = a/b

    r['PP'] = PupPdn
    r['PS'] = PupSdn
    r['SP'] = SupPdn
    r['SS'] = SupSdn

    return r


if __name__ == "__main__":
    from numpy import arange
    vp1=10.0
    vs1=5.0
    den1=3.0
    vp2=10.0
    vs2=2.0
    den2=3.0

    for hslow in arange(0.00,0.12,0.01):
        #rt = rtcoef(vp1, vs1, den1, vp2, vs2, den2, hslow)
        #rf = rcoef_free(vp2, vs2, hslow)
        rf = rtcoefSH(vs1,den1,vs2,den2,hslow)
        #a = (rf['PS']-rt[12])/rf['PS']
        #b = (rf['SS']-rt[16])/rf['SS']

        print(hslow, rf)

