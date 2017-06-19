from algorithms import *

#initializing global constants
const = constants()        

#initializing a ship instance
Asterious = Ship()

#initializing the Rolling Condition
Roll_Angle = 1 #degrees
Roll_Period = 5 #seconds
Wave_Direction = 90 #degrees
Vessel_Speed = 1 #m/s
RC = rollingCondition(const,Roll_Angle,Roll_Period,Wave_Direction,Vessel_Speed)

#defining the ship particulars
Asterious.setParticulars(170,32,12,6.7,0.846,0.998,13.32,3.80,1,const)

#adding objects to our ship instance
Asterious.addBilgeKeel(67.2,0.416,49)
Asterious.addMoonpool(4.2,3.75,6.7,const.rho.value,const.g.value)
Asterious.addAntiRollingTank(44.8,3,8.6,const.rho.value,const.g.value)
Asterious.addSections("Asterious.txt")
Asterious.sections.formFrames(Asterious.getOG())
#Asterious.sections.plotFrames()
#Asterious.sections.plotSectionalArea()
Asterious.CalculateDamping(const,RC)

print Asterious.damping.B44.value