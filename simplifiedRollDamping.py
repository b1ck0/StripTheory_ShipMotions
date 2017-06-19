from math import pi,sqrt,exp, log

class SimplifiendDamping:
    def __init__(self, LPP, LB, BD, CB, CMID, OGD, PHI, TW, LBKL, BBKB):
        # initializing constants
        self.RO = 102
        self.KVC = 1.14e-6 #kinematic viscosity coefficient
         
        # initializing the ship particulars
        self.LPP = LPP
        self.LB = LB
        self.BD = BD
        self.CB = CB
        self.CMID = CMID
        self.OGD = OGD
        self.PHI = PHI # in degrees
        self.TW = TW # in seconds
        self.OMEGA = 2*pi/TW
        self.BRTH = LPP/LB
        self.DRAFT = self.BRTH/BD
        self.OMEGAHAT = self.OMEGA*sqrt(self.BRTH/2.0/9.81)
        
        #initializing the bilge-keel particulars
        self.LBKL = LBKL
        self.BBKB = BBKB
        
        #initializing the damping components
        self.BFHAT = 0.0 #frictional damping
        self.BWHAT = 0.0 #wave-making damping
        self.BEHAT = 0.0 #eddy-making damping
        self.BBKHAT = 0.0 #bilge-keel damping
        self.B44HAT = 0.0 #total damping
        
        # initializing errors
        self.error = []
        
        if self.OMEGAHAT > 1.0:
            self.error.append("OMEGAHAT -> Out of Range!")
            
        if self.OGD > 0.2 or self.OGD < -1.5:
            self.error.append("OGD -> Out of Range!")
            
        if self.CMID < 0.9 or self.CMID > 0.99:
            self.error.append("CMID -> Out of Range!")
        
        if self.CB < 0.5 or self.CB > 0.85:
            self.error.append("CB -> Out of Range!")
            
        if self.BD < 2.5 or self.BD > 4.5:
            self.error.append("BD -> Out of Range!")
        
        if self.LBKL < 0.05 or self.LBKL > 0.4:
            self.error.append("LBKL -> Out of Range!")
            
        if self.BBKB < 0.01 or self.BBKB > 0.06:
            self.error.append("BBKB -> Out of Range!")
            
        return self.error
        
    def frictionalDamping(self):
        RF = self.DRAFT*((0.887+0.145*self.CB)*(1.7+self.CB*self.BD)-2.0*self.OGD)/pi
        SF = self.LPP*(1.75*self.DRAFT+self.CB*self.BRTH)
        CF = 1.328*((3.22*RF**2*(self.PHI*pi/180.0)**2)/(self.TW*self.KVC))**-0.5
        BF = 4.0/3.0/pi*self.RO*SF*RF**3*(self.PHI*pi/180.0)*self.OMEGA*CF
        self.BFHAT = BF/(self.RO*self.LPP*self.BRTH**3*self.DRAFT*self.CB)*sqrt(self.BRTH/2.0/9.81)
        
    def waveDamping(self):
        X1 = self.BD
        X2 = self.CB
        X3 = self.CMID
        X5 = self.OMEGAHAT
        X4 = 1.0 - self.OGD

        A111 = -0.002222*X1**3+0.040871*X1**2-0.286866*X1+0.599424
        A112 = 0.010185*X1**3-0.161176*X1**2+0.904989*X1 -1.641389
        A113 = -0.015422*X1**3+0.220371*X1**2-1.084987*X1+1.834167
        A121 = -0.0628667*X1**4+0.4989259*X1**3+0.52735*X1**2-10.7918672*X1+16.616327
        A122 = 0.1140667*X1**4-0.8108963*X1**3-2.2186833*X1**2+25.1269741*X1-37.7729778
        A123 = -0.0589333*X1**4+0.2639704*X1**3+3.1949667*X1**2-21.8126569*X1+31.4113508
        A124 = 0.0107667*X1**4+0.0018704*X1**3-1.2494083*X1**2+6.9427931*X1-10.2018992
        A131 = 0.192207*X1**3-2.787462*X1**2+12.507855*X1-14.764856
        A132 = -0.350563*X1**3+5.222348*X1**2-23.974852*X1+29.007851
        A133 = 0.237096*X1**3-3.535062*X1**2+16.368376*X1-20.539908
        A134 = -0.067119*X1**3+0.966362*X1**2-4.407535*X1+5.894703

        A11=A111*X2**2+A112*X2+A113
        A12=A121*X2**3+A122*X2**2+A123*X2+A124
        A13=A131*X2**3+A132*X2**2+A133*X2+A134

        AA111=17.945*X1**3-166.294*X1**2+489.799*X1-493.142
        AA112=-25.507*X1**3+236.275*X1**2-698.683*X1+701.494
        AA113=9.077*X1**3-84.332*X1**2+249.983*X1-250.787
        AA121=-16.872*X1**3+156.399*X1**2-460.689*X1+463.848
        AA122=24.015*X1**3-222.507*X1**2+658.027*X1-660.665
        AA123=-8.56*X1**3+79.549*X1**2-235.827*X1+236.579

        AA11=AA111*X2**2+AA112*X2+AA113
        AA12=AA121*X2**2+AA122*X2+AA123

        AA1=(AA11*X3+AA12)*(1-X4)+1.0

        A1=(A11*X4**2+A12*X4+A13)*AA1
        A2=-1.402*X4**3+7.189*X4**2-10.993*X4+9.45

        A31=-7686.0287*X2**6+30131.5678*X2**5-49048.9664*X2**4+42480.7709*X2**3-20665.147*X2**2+5355.2035*X2-577.8827
        A32=61639.9103*X2**6-241201.0598*X2**5+392579.5937*X2**4-340629.4699*X2**3+166348.6917*X2**2-43358.7938*X2+4714.7918
        A33=-130677.4903*X2**6+507996.2604*X2**5-826728.7127*X2**4+722677.104*X2**3-358360.7392*X2**2+95501.4948*X2-10682.8619
        A34=-110034.6584*X2**6+446051.22*X2**5-724186.4643*X2**4+599411.9264*X2**3-264294.7189*X2**2+58039.7328*X2-4774.6414
        A35=709672.0656*X2**6-2803850.2395*X2**5+4553780.5017*X2**4-3888378.9905*X2**3+1839829.259*X2**2-457313.6939*X2+46600.823
        A36=-822735.9289*X2**6+3238899.7308*X2**5-5256636.5472*X2**4+4500543.147*X2**3-2143487.3508*X2**2+538548.1194*X2-55751.1528
        A37=299122.8727*X2**6-1175773.1606*X2**5+1907356.1357*X2**4-1634256.8172*X2**3+780020.9393*X2**2-196679.7143*X2+20467.0904
        
        AA311=(-17.102*X2**3+41.495*X2**2-33.234*X2+8.8007)*X4+36.566*X2**3-89.203*X2**2+71.8*X2-18.108

        AA31=(-0.3767*X1**3+3.39*X1**2-10.356*X1+11.588)*AA311
        AA32=-0.0727*X1**2+0.7*X1-1.2818

        XX4=X4-AA32

        AA3=AA31*(-1.05584*XX4**9+12.688*XX4**8-63.70534*XX4**7+172.84571*XX4**6-274.05701*XX4**5+257.68705*XX4**4-141.40915*XX4**3+44.13177*XX4**2-7.1654*XX4-0.0495*X1**2+0.4518*X1-0.61655)

        A3=A31*X4**6+A32*X4**5+A33*X4**4+A34*X4**3+A35*X4**2+A36*X4+A37+AA3
        
        self.BWHAT=A1/X5*exp(-A2*(log(X5)-A3)**2/1.44)