# Ship Motion Analysis Application - Strip Theory
# The application estimates the viscous roll and surge dampings for a range of sea-states

import numpy as np
from math import *
from scipy.spatial import ConvexHull
from sympy.geometry import *
from operator import *
from scipy.integrate import simps

from matplotlib import pyplot
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from prettytable import PrettyTable

#defining rmax function used in eddymaking damping
class funct:
    def rmax(self,psi,B,a1,a3):
        M = np.round(B/(2*(1+a1+a3)),3)
        return np.round(M*sqrt(((1+a1)*sin(psi)-a3*sin(3*psi))**2+((1-a1)*cos(psi)+a3*cos(3*psi))**2),3)
                         

#container for labeled values
class variable:
    def __init__(self,value, label):
        self.value = value
        self.label = label

#container for physical constants        
class constants:
    def __init__(self):
        self.rho = variable(1025,"Sea Water Density [kg/m^3]")
        self.t = variable(10,"Sea Water Temperature [deg C]")
        self.nu = variable(1.307e-5,"Sea Water Kinematic Viscosity [m^2/s]")
        self.g = variable(9.81,"Acceleration due to Gravity [m/s^2]")
    
    def getRho(self):
        return self.rho.value
    
    def getNu(self):
        return self.nu.value
        
    def getG(self):
        return self.g.value

#container for rolling condition
class rollingCondition:
    def __init__(self,const,rollAmplitude,rollPeriod,waveDirection,speed):
        g = const.getG()
        self.fa_d = variable(rollAmplitude,'Roll Amplitude [deg]')
        self.fa_r = variable(np.round(rollAmplitude*pi/180,6),'Roll Amplitude [rad]')
        self.Tr = variable(rollPeriod,'Roll Period [sec]')
        self.X = variable(waveDirection,'Wave Direction [deg]')
        self.V = variable(speed,"Ship Forward Speed in [m/s]")
        self.we = variable(np.round(speed*(((2*pi/rollPeriod)**2)/g)*cos(radians(waveDirection))+(2*pi/rollPeriod),3),'Roll Frequency [rad/s]')
        
#container for bilgekeel data        
class bilgeKeel:
    def __init__(self,L,B,x0):
        self.L = variable(L,"Length of the bilge keel [m]")
        self.B = variable(B,"Breadth of the bilge keel [m]")
        self.x0 = variable(x0,"Bilge keel start coordinate [m]")
        self.x1 = variable(x0+L,"Bilge keel end coordinate [m]")
           
    def getL(self):
        return self.L.value
        
    def getB(self):
        return self.B.value

#container for moonpool data        
class moonpool:
    def __init__(self,L,B,h,rho,g):
        self.L = variable(L,"Length of the Moonpool in X-Direction (Local Ship Coordinates) [m]")
        self.B = variable(B,"Length of the Moonpool in Y-Direction (Local Ship Coordinates) [m]")
        self.h = variable(h,"Height of the Water inside of the moonpool [m]")
        self.w = variable(np.round((pi/B)*sqrt(g*h),3),"Fluid Natural Frequency [rad/s]")

#container for anti-rolling tank data        
class antiRollingTank:
    def __init__(self,L,B,h,rho,g):
        self.L = variable(L,"Length of the Anti-Rolling Tank in X-Direction (Local Ship Coordinates) [m]")
        self.B = variable(B,"Length of the Anti-Rolling Tank  in Y-Direction (Local Ship Coordinates) [m]")
        self.h = variable(h,"Height of the Water inside of the Anti-Rolling Tank [m]")
        self.w = variable(np.round((pi/B)*sqrt(g*h),3),"Fluid Natural Frequency [rad/s]")

#container for hydrostatics stiffness matrix                                
class hydrostaticStiffness:
    def __init__(self):
        self.matrix = np.zeros(shape=(6,6))
        
    def setC44(self,const,particulars):
        self.matrix[4,4] = np.round(const.g.value * particulars.Displ.value * particulars.GM.value,3) # [kg.m^2/s^2] --> [N.m]
        
    def getC44(self):
        return self.matrix[4,4]
        
class damping:
    def __init__(self):
        self.irregularMatrix = np.zeros(shape=(6,6))
        self.regularMatrix =  np.zeros(shape=(6,6))
        self.matrix = np.zeros(shape=(6,6))
        self.B44F = variable(0.0,'Skin-friction Roll Damping Coefficient [N.m.sec/rad]')
        self.B44L = variable(0.0,'Lift Roll Damping Coefficient [N.m.sec/rad]')
        self.B44E = variable(0.0,'Eddy-making Roll Damping Coefficient [N.m.sec/rad]')
        self.B44BKL = variable(0.0,'Bilge-keel Lift Damping Coefficient [N.m.sec/rad]')
        self.B44BKN = variable(0.0,'Bilge-keel Normal Force Damping Coefficient [N.m.sec/rad]')
        self.B44BKH = variable(0.0,'Bilge-keel Pressure Damping Coefficient [N.m.sec/rad]')
        self.B44IWT = variable(0.0,'Anti-rolling Tanks Damping Coefficient [N.m.sec/rad]')
        self.B44IWM = variable(0.0,'Moonpool Damping Coefficient [N.m.sec/rad]')
        self.B44W = variable(0.0,'Wave-making Damping Coefficient [N.m.sec/rad]')
        self.B44 = variable(0.0,"Total Roll Damping Coefficient [N.m.sec/rad]")
        self.B22 = variable(0.0,"Sway Damping Coefficient")
        self.B42 = variable(0.0,"Coupling Damping Coefficient of Roll into Sway")
        self.sect_B44F = []
        self.sect_B44E = []
        self.sect_B44BKN = []
        self.sect_B44BKH = []
    
    def calculateLiftDamping(self,const,L,B,d,Cm,OG,V):       
        rho = const.rho.value
        if Cm <= 0.92:
            k = 0
        elif Cm <= 0.97:
            k = 0.1
        else:
            k = 0.3   
        
        kN = np.round(2*pi*(d/L)*k*(4.1*(B/L)-0.045),3) #lift slope
        L0 = 0.3*d #lever L0*fa/V representing angle of attack of lifting body
        LR = 0.5*d #lever from O to centre of lift force
        
        B44 = np.round(0.5*rho*V*L*d*kN*L0*LR*(1-1.4*(OG/LR)+0.7*OG/(L0*LR)),3)
        
        self.B44L.value = iadd(self.B44L.value,B44)
        self.B44.value = iadd(self.B44.value,B44)
        
    def calculateMoonpoolRollDamping(self,const,moonpool,rollingCondition,OG,Bship):
        L = moonpool.L.value
        Bcomp = moonpool.B.value
        h = moonpool.h.value
        fa = rollingCondition.fa_r.value
        we = rollingCondition.we.value
        g = const.g.value
        rho = const.rho.value

        A = np.round((1.8*h/Bcomp-1.9882*fa+0.429)/(1.2*OG/Bship+1),3)
        B = np.round(40.842*h/Bcomp-10.502*fa+2.1,3)
        C = np.round(1/pi*sqrt(Bship/g)*(we/sqrt(h/Bcomp)),3)

        B44 = np.round(A*(C**(B))*exp(-C**(B))*L*rho*(Bcomp**5)*sqrt(2*g/Bcomp),3) 
        
        self.B44IWM =iadd(self.B44IWM.value,B44)
        self.B44.value = iadd(self.B44.value,B44)

    def calculateAntiRollingTankDamping(self,const,antiRollingTank,rollingCondition,OG,Bship):
        L = antiRollingTank.L.value
        Bcomp = antiRollingTank.B.value
        h = antiRollingTank.h.value
        fa = rollingCondition.fa_r.value
        we = rollingCondition.we.value
        rho = const.rho.value
        g = const.g.value        

        A = np.round((1.8*h/Bcomp-1.9882*fa+0.429)/(1.2*OG/Bship+1),3)
        B = np.round(40.842*h/Bcomp-10.502*fa+2.1,3)
        C = np.round(1/pi*sqrt(Bship/g)*(we/sqrt(h/Bcomp)),3)
        
        B44 = np.round(A*(C**(B))*exp(-C**(B))*L*rho*(Bcomp**5)*sqrt(2*g/Bcomp),3)        
        
        self.B44IWT = iadd(self.B44IWT.value,B44)
        self.B44.value = iadd(self.B44.value,B44)
        
    def calculateHullFrictionDamping(self,sections,particulars,rollingCondition,const):
        OG = particulars.OG.value
        Cb = particulars.Cb.value
        L = particulars.L.value
        Sf = particulars.SF.value
        fa = rollingCondition.fa_r.value
        nu = const.nu.value
        rho = const.rho.value
        we = rollingCondition.we.value
        V = rollingCondition.V.value
        
        x = []
        y = []

        rf = np.round(((0.887+0.145*Cb)*(Sf/L)-2*OG)/pi,3)
        Rn = np.round((0.512*(rf*fa)**2*we)/nu,3)
        Cf = np.round(1.328*Rn**(-0.5)+0.014*Rn**(-0.114),3)

        B44F = np.round(0.5*rho*rf**3*Sf*Cf*(1+4.1*V/(we*L)),3) #non-linear fricntion damping component
        B44FEQ = np.round(B44F*((8/(3*pi))*fa*we),3) #linerized friction damping component
        
        for frame in sections.frames:
            self.sect_B44F.append([frame.xCoord,B44FEQ])
            x.append(frame.xCoord)
            y.append(B44FEQ)
        
        self.B44F.value = np.round(simps(y,x),3)
        self.B44.value = iadd(self.B44.value,self.B44F.value)
        
    def calculateEddyMakingDamping(self,sections,particulars,rollingCondition,const):
        OG = particulars.OG.value
        L = particulars.L.value
        fa = rollingCondition.fa_r.value
        rho = const.rho.value
        we = rollingCondition.we.value
        V = rollingCondition.V.value
        B44E = 0
        
        t = PrettyTable(['x-Coord','f1','f2','f3','rmax','R','l','psi','A0','B0','H','H0"','sig"','M','j','Cp','CR',"B44E'"])
        
        x = []
        y = []
        
        for frame in sections.frames:
            sig_s = frame.areaCoeff
            Bs = frame.Bs
            Ds = frame.Ds
            H0 = frame.H0
            a1 = frame.a1
            a3 = frame.a3
            R = frame.R
            l = frame.l
            
            f1 = np.round(0.5*(1+tanh(20*(sig_s-0.7))),3)
            f2 = np.round(0.5*(1-cos(pi*sig_s))-1.5*(1-exp(-5*(1-sig_s)))*((sin(pi*sig_s))**2),3)
            f3 = np.round(1+4*exp(-1.65E5*(1-sig_s)**2),3)
        
            psi1 = 0.0
            try:
                psi2 = np.round(0.5*acos((a1*(1+a3))/(4*a3)),3)
            except ValueError:
                psi2 = 0.0

            f = funct()
            rmax_psi1 = np.round(f.rmax(psi1,Bs,a1,a3),3)
            rmax_psi2 = np.round(f.rmax(psi2,Bs,a1,a3),3)

            if rmax_psi1 >= rmax_psi2:
                rmax = rmax_psi1
                psi = psi1
            else:
                rmax = rmax_psi2
                psi = psi2

            A0 = np.round(-2*a3*cos(5*psi)+a1*(1-a3)*cos(3*psi)+((6-3*a1)*a3**2+(a1**2-3*a1**2)*a3+a1**2)*cos(psi),3)
            B0 = np.round(-2*a3*sin(5*psi)+a1*(1-a3)*sin(3*psi)+((6+3*a1)*a3**2+(3*a1+a1**2)*a3+a1**2)*sin(psi),3)
            H = np.round(1+a1**2+9*a3**2+2*a1*(1-3*a3)*cos(2*psi)-6*a3*cos(4*psi),3)
            sig_s_star = np.round((sig_s-OG/Ds)/(1-OG/Ds),3)
            H0_star = np.round(H0/(1-OG/Ds),3)
            M = np.round(frame.Ms,3)
            gamma = np.round(((pi**(1/2)*f3)/(2*Ds*(H0*(sig_s+OG/Ds))**(1/2)))*(rmax+(2*M/H)*(A0**2+B0**2)**(1/2)),3)
            Cp = np.round(0.5*(0.87*exp(-gamma)-4*exp(-(0.187*gamma))+3),3)
            CR = np.round(((1-f1*R/Ds)*(1-OG/Ds)+f2*(H0-f1*R/Ds)**2)*Cp*(rmax/Ds)**2,3)
            K = V/(0.04*we*L)
            
            B44E = np.round((0.5*rho*Ds**4*CR)/(1+K**2),3) #non-linear eddy-making damping
            B44EEQ = np.round(B44E*((8/(3*pi))*fa*we),3) #linerized eddy-making damping
            
            t.add_row([frame.xCoord,f1,f2,f3,rmax,R,l,psi,A0,B0,H,H0_star,sig_s_star,M,gamma,Cp,CR,B44E])
            self.sect_B44E.append([frame.xCoord,B44EEQ])
            x.append(frame.xCoord)
            y.append(B44EEQ)
            
        self.B44E.value = np.round(simps(y,x),3)
        self.B44.value = iadd(self.B44.value,self.B44E.value)
        
        print t

    def calculateCriticalDamping(self,particulars,hydrostaticStiffness,addedMass):
        A44 = addedMass.matrix[4,4]
        C44 = hydrostaticStiffness.matrix[4,4]
        I44 = np.round((0.4*particilars.B.value)**2*particulars.Displ.value,3)
        M44 = I44 + A44
        
        return np.round(2*sqrt(M44*C44),3)

    def calculateBilgekeelDamping(self,sections,particulars,rollingCondition,const,bilgekeel):
        bbk = bilgekeel.B.value
        fa = rollingCondition.fa_r.value
        we = rollingCondition.we.value
        rho = const.rho.value
        OG = particulars.OG.value
        start = bilgekeel.x0.value
        end = bilgekeel.x1.value
        
        x = [start] #adding the start of the bilge-keel as an integration point
        y1 = []
        y2 = []
        
        for frame in sections.frames:
            if frame.xCoord >= start and frame.xCoord <= end:
                sig = frame.areaCoeff
                Ds = frame.Ds
                R = frame.R
                l = frame.l
                H0 = frame.H0
                f = np.round(1+0.3*exp(-160*(1-sig)),3)
                Cd = np.round(22.5*(bbk/(pi*l*fa*f))+2.4,3)
                
                B44BKN = np.round(rho*l**3*bbk*f**2*Cd,3) #non-linear bilge-keel normal damping
                B44BKNEQ = np.round(B44BKN*((8/(3*pi))*fa*we),3) #linerized bilge-keel normal damping
    
                S0 = np.round(0.3*pi*f*l*fa+1.95*bbk,3)
                m1 = np.round(R/Ds,3)
                m2 = np.round(OG/Ds,3)
                m3 = np.round(1-m1-m2,3)
                m4 = np.round(H0 - m1,3)
                m5 = np.round((0.414*H0+0.651*m1**2-(0.382*H0+0.0106)*m1)/((H0-0.215*m1)*(1-0.215*m1)),3)
                m6 = np.round((0.414*H0+0.651*m1**2-(0.382+0.0106*H0)*m1)/((H0-0.215*m1)*(1-0.215*m1)),3)
                m7 = 0
                m8 = np.round(m7 + 0.414*m1*(1-cos(S0/R)),3)
                if S0 > 0.25*pi*R:
                    m7 = np.round(S0/Ds-0.25*pi*m1,3)
                    m8 = np.round(m7 + 0.414*m1,3)
                B0 = np.round(m2**2/(3*(H0-0.215*m1))+((1-m1)**2*(2*m3-m2))/(6*(1-0.215*m1))+m1*(m3*m5+m4*m6),3)
                A0 = np.round((m3+m4)*m8-m7**2,3)
                CPP = 1.2
                CPM = np.round(-22.5*(bbk/(pi*l*fa*f))-1.2,3)
                integral = Ds**2*(-A0*CPM+B0*CPP)
    
                B44BKH = np.round(0.5*rho*l**2*f**2*integral,3) #non-linear bilge-keel pressure damping
                B44BKHEQ = np.round(B44BKH*((8/(3*pi))*fa*we),3) #linerized bilge-keel pressure damping
                
                x.append(frame.xCoord)
                y1.append(B44BKNEQ)
                y2.append(B44BKHEQ)
                self.sect_B44BKN.append([frame.xCoord,B44BKNEQ])
                self.sect_B44BKH.append([frame.xCoord,B44BKHEQ])
        
        x.append(end) #adding the end of the moonpool as an integration point
        y1.insert(0,y1[0])
        y2.insert(0,y2[0])
        y1.append(y1[-1])
        y2.append(y2[-1])
        self.B44BKN.value = np.round(simps(y1,x),3)
        self.B44BKH.value = np.round(simps(y2,x),3)
        
        self.B44.value = iadd(self.B44.value,self.B44BKN.value)
        self.B44.value = iadd(self.B44.value,self.B44BKH.value)
              
#container for ship particulars
class particulars:
    def __init__(self,L,B,D,d,Cb,Cm,KG,GM,V,const):
        self.L = variable(L,"Ship Length [m]")
        self.B = variable(B,"Ship Breadth [m]")
        self.D = variable(D,"Ship Depth [m]")
        self.d = variable(d,"Ship Draft [m]")
        self.Cb = variable(Cb,"Ship Block Coefficient [-]")
        self.Cm = variable(Cm,"Ship Midship Coefficient [-]")
        self.KG = variable(KG,"Ship Vertical Center of Gravity Coordinate [m]")
        self.GM = variable(GM,"Ship GM [m]")
        
        self.SF = variable(L*(1.7*d+Cb*B),"Hull Wetted Surface [m^2]")
        self.OG = variable(d-KG,"Vertical Distance from Mean Waterline to G [m]")
        self.Displ = variable(const.rho.value*L*B*d*Cb,"Ship Displacement [kg]")

    def getSF(self):
        return self.SF.value                   
                                                         
    def getDispl(self):
        return self.Displ.value
        
    def getOG(self):
        return self.OG.value
        
    def getDraft(self):
        return self.d.value

#container for 3D point    
class point3D:
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

#container for a ship frame with its points and particulars        
class frame:
    def __init__(self,xCoord):
        self.xCoord = xCoord
        self.area = 0.0
        self.perimeter = 0.0
        self.areaCoeff = 0.0
        self.points = []
        self.pointsUW = []
        self.qHullVertices = []
        
    def defineFramePosition(self,xCoord):
        self.xCoord = xCoord
    
    def printFrame(self):
        print '-------------------------/// NEW FRAME ///-------------------------'
        print 'Frame Position: {0} '.format(self.xCoord)
        print 'Number of Points: {0} '.format(len(self.points))
        print 'Bs: {0}'.format(self.Bs)
        print 'Ds: {0}'.format(self.Ds)
        print 'Area: {0} '.format(self.area)
        print 'sig: {0}'.format(self.areaCoeff)
        print 'H0: {0}'.format(self.H0)
        print 'Ms: {0}'.format(self.Ms)
        print 'a1: {0}'.format(self.a1)
        print 'a3: {0}'.format(self.a3)
        print 'R: {0}'.format(self.R)
        print 'l: {0}'.format(self.l)
         
    def LewisTransform(self,OG):
        self.H0 = np.round((self.Bs/2)/self.Ds,3)
        
        Bs = self.Bs
        Ds = self.Ds
        H0 = self.H0
        sig_s = self.areaCoeff
        
        self.c1 = np.round((3+(4*sig_s)/pi)+(1-(4*sig_s)/pi)*((H0-1)/(H0+1))**2,3)
        self.c2 = np.round(2*(self.c1-3),3)
        self.c3 = np.round(self.c1-4,3)
        
        self.a3 = np.round((-self.c1+3+(9-2*self.c1)**(1/2))/self.c1,3)
        self.a1 = np.round(((H0-1)/(H0+1))*(self.a3+1),3)
        
        self.Ms = np.round((Bs/2)/(1+self.a1+self.a3),3)
        
        self.lam_a = np.round(1+self.a1+self.a3,3)
        self.lam_b = np.round(1-self.a1+self.a3,3)
        self.b0 = np.round(2*self.Ms*self.lam_a,3)
        self.d0 = np.round(self.Ms*self.lam_b,3)

        R = np.round(2*Ds*sqrt((H0*(sig_s-1.0))/(pi-4)),3)
        
        if (H0>=1) and (R/Ds>1):
            R = Ds
        if (H0<=1) and (R/Ds>H0):
            R = np.round(Bs/2,3)
        
        l = np.round(Ds*sqrt((H0-(1-sqrt(2)/2)*(R/Ds))**2+(1-(OG/Ds)-(1-sqrt(2)/2)*(R/Ds))**2),3)       

        self.R = R
        self.l = l
           
    def calculateArea(self):
        y = []
        z = []
        for point in self.qHullVertices:
            y.append(np.round(point.y,4))
            z.append(np.round(point.z,4))
            
        y1 = list(reversed(y))
        z1 = list(reversed(z))

        self.Bs = np.amax(y1)*2
        self.Ds = np.amax(z1) - np.amin(z1)
                
        self.area = np.round(np.abs(np.dot(y1,np.roll(z1,1)) - np.dot(z1,np.roll(y1,1))))
        self.areaCoeff = np.round(self.area/(self.Bs*self.Ds),3)
        
    def slice(self,draft):
        self.points.append(point3D(self.xCoord,0.0,draft))
        self.draft = draft
        
        below = []
        above = []
        
        maxBelow = Point3D(0,0,0)
        minAbove = Point3D(0,0,0)
        
        for point in self.points:
            if point.z <= draft:
                below.append(point)
                if (point.z > maxBelow.z) or ((point.z == maxBelow.z) and (point.y > maxBelow.y)):
                    maxBelow = point
            elif point.z > draft:
                above.append(point)
                if (point.z < minAbove.z) or ((point.z == minAbove.z) and (point.y < minAbove.y)):
                    minAbove = point
        
        if (maxBelow.z - draft > 0.01) and (np.amax([point2.y for point2 in below]) - maxBelow.y > 0.001):
            dz = minAbove.z - maxBelow.z
            dy = minAbove.y - maxBelow.y
            slope = dz/dy
            below.append(point3D(self.xCoord,maxBelow.y + (draft-maxBelow.z)/slope,draft))
        else:
            below.append(point3D(self.xCoord,np.amax([point3.y for point3 in below]),draft))

        y = []
        z = []
        arr = []
        for pointA in below:
            y.append(pointA.y)
            z.append(pointA.z)

        for i in range(0,len(y)):
            arr.append([y[i],z[i]])
       
        hull = ConvexHull(arr) 
        
        for vertex in hull.vertices:
            x1 = self.xCoord
            y1 = hull.points[vertex][0]
            z1 = hull.points[vertex][1]
            self.qHullVertices.append(point3D(x1,y1,z1))
                               
#container for ship sections (input file location, points and frames)        
class sections:
    def __init__(self,filename,draft):
        self.filename = filename
        self.points = []
        self.frames = []
        self.draft = draft
                       
        with open(filename,'r') as textFile:
            for line in textFile:
                xyz = line.split(',')
                self.points.append(point3D(np.round(float(xyz[0]),6),np.round(float(xyz[1]),6),np.round(float(xyz[2]),6)))
        
    def showSections(self):
        for element in self.points:
            print element.x, element.y, element.z
            
    def formFrames(self,OG):
        tempArray = []
        for element in self.points:
            tempArray.append(element.x)
        
        uniqueXCoords = sorted(list(set(tempArray)))
        
        for xCoord in uniqueXCoords:
            self.frames.append(frame(xCoord))
        
        for i in range(len(self.frames)):
            for point in self.points:
                if point.x == self.frames[i].xCoord:
                    self.frames[i].points.append(point)
                    
        for frame1 in self.frames:
            frame1.slice(self.draft)
            frame1.calculateArea()
            frame1.LewisTransform(OG)
                    
    def plotPoints(self):       
        figure = pylab.figure()
        screen = Axes3D(figure)
        
        for frame in self.frames:
            x = []
            y = []
            z = []
            for point in frame.points:
                x.append(point.x)
                y.append(point.y)
                z.append(point.z)
            screen.plot(x,y,z)

        pyplot.show()
        
    def plotFrames(self):
        figure = pylab.figure()
        screen = Axes3D(figure)
        
        for frame in self.frames:
            x = []
            y = []
            z = []
            for point in frame.qHullVertices:
                x.append(point.x)
                y.append(point.y)
                z.append(point.z)
            screen.plot(x,y,z)

        pyplot.show()
        
    def plotSectionalArea(self):
        y = []
        x = []
        for frame in self.frames:
            x.append(frame.xCoord)
            y.append(frame.area)
        
        plt.figure(2)
        plt.plot(x,y)
        plt.show()
        
        
#container for all ship data                                    
class Ship:
    def __init__(self):
        self.bilgeKeel = []
        self.moonpool = []
        self.antiRollingTank = []
        
    def setParticulars(self,L,B,D,d,Cb,Cm,KG,GM,V,const):
        self.particulars = particulars(L,B,D,d,Cb,Cm,KG,GM,V,const)
              
    def addBilgeKeel(self,L,B,x0):
        self.bilgeKeel.append(bilgeKeel(L,B,x0))   
        
    def addMoonpool(self,L,B,h,rho,g):
        self.moonpool.append(moonpool(L,B,h,rho,g)) 
        
    def addAntiRollingTank(self,L,B,h,rho,g):
        self.antiRollingTank.append(antiRollingTank(L,B,h,rho,g))
        
    def addHydrostaticStiffnessMatrix(self,matrix):
        self.HydrostaticStiffnessMatrix = matrix
    
    def addSections(self,filename):
        self.sections = sections(filename,self.getDraft())
        
    def getDraft(self):
        return self.particulars.getDraft()
        
    def getOG(self):
        return self.particulars.getOG()
        
    def CalculateDamping(self,const,rollingCondition):
        self.damping = damping()
        L = self.particulars.L.value
        B = self.particulars.B.value
        d = self.particulars.d.value
        Cm = self.particulars.Cm.value
        OG = self.particulars.OG.value
        V = rollingCondition.V.value
        
        self.damping.calculateLiftDamping(const,L,B,d,Cm,OG,V)
        self.damping.calculateHullFrictionDamping(self.sections,self.particulars,rollingCondition,const)
        self.damping.calculateEddyMakingDamping(self.sections,self.particulars,rollingCondition,const)
        
        for bilgekeel in self.bilgeKeel:
            self.damping.calculateBilgekeelDamping(self.sections,self.particulars,rollingCondition,const,bilgekeel)
       
        for moonpool in self.moonpool:
            self.damping.calculateMoonpoolRollDamping(const,moonpool,rollingCondition,OG,B)

        for antiRollingTank in self.antiRollingTank:
            self.damping.calculateAntiRollingTankDamping(const,antiRollingTank,rollingCondition,OG,B)