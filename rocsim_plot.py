#
#
# CLEAR SCREEN, IMPORT MODULES SECTION
#
# clear the screen
import subprocess as sp
tmp =sp.call('clear',shell=True)

import plotly
import plotly.plotly as py
import plotly.graph_objs as go
import pandas as pd
import numpy as np
from decimal import Decimal
from math import pi
from ast import literal_eval

# CONSTANTS SECTION
#
# conversion factor for pounds from newtons
pounds_per_newton = 0.224809
# acceleration due to gravity in meters per second squared
Acc_Grav = 9.8101
# granularity of fit to thrust curve data in seconds
granularity = 0.02

# FUNCTIONS SECTION

# function returns true if string is a number
def is_numeric(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

# ensures numeric input
def prompt(clue):
	while True:
		y = raw_input(clue)
		if  is_numeric(y):
			return float(y)
		else:
			print y, "is not numeric"

# function to return engine type with capital letter
def cap_engine(type):
        if (ord(type[0]) <= ord('Z')):
                return type
        else:
                return chr(ord(type[0]) - 32) + type[1:]

# function for converting fahrenheit to celsius
def F_to_C(degF):
        return(degF - 32.0) * (5./9.)

# function for converting celsius to fahrenheit
def C_to_F(degC):
        return((9./5.) * degC + 32)
    
# function for calculating air density, temperature in centigrade
def air_density(temperature):
        kelvin = temperature + 273.15
        Pa = 101325 # standard sea level air pressure
        R = 287.058 #gas constant for dry air
        return(Pa/(R * kelvin))

# function for truncating a number to two digits after decimal
def two_digits(num):
        x=Decimal(num)
        return round(x, 2)

# function to display the current engine types
def display_engines():
        print "Engines: ",
        temp_type = engine_type
        num_engines = len(temp_type)

        for i in range (0, num_engines-1):
                for j in range (i+1, num_engines):
                        if temp_type[i][0] == temp_type[j][0]:
                                if (int(temp_type[j][1:]) < int(temp_type[i][1:])):
                                        temp_type[i], temp_type[j] = temp_type[j], temp_type[i]

        for engine in temp_type:
                print engine,

# user input for simulation
#
# rocket mass
rocket_mass_gms = prompt("Enter mass of empty rocket, no motor or payload, in grams: ")

# payload mass
payload_mass_gms = prompt("Enter mass of payload, in grams: ")

# rocket diameter
dia_mm = prompt("Enter rocket's maximum body tube diameter, in mm: ")

# input the drag coefficient (unitless number - usually 0.60 to 0.75)
Coefficient_drag = prompt("Enter drag coefficient: ")

# input the temperature , if Fahrenheit convert to centigrade
Temperature = int (prompt("Enter launch temperature: "))
while True:
    opt = raw_input("(F)ahrenheit or (C)entigrade: ")
    if (opt in ['F','C','f','c']):
        break
    else:
        print "Please enter 'F' or 'C'"
if (opt in ['F','f']):
    Temperature = F_to_C(Temperature)

# READ DATA SECTION
#
print 'Reading motor mass data'

# open the motor mass data file
data_file = open('./motor_mass.dat','r')

data=data_file.readline()
motor_mass_dictionary = literal_eval(data)
motor_mass_gms = pd.Series(motor_mass_dictionary)

# close the data file
data_file.close()

print 'Reading propellant mass data'

# open the propellant mass data file
data_file = open('./propellant_mass.dat','r')

data=data_file.readline()
propellant_mass_dictionary = literal_eval(data)
propellant_mass_gms = pd.Series(propellant_mass_dictionary)

# close the data file
data_file.close()

engine_type = list(motor_mass_gms.index)

print
display_engines()
print
print

engine = raw_input("Enter engine: ")
engine = cap_engine(engine)

print 'Reading thrust data for %s.'%(engine)
# read the excel thrust data for engine
while True:
        
        try:
                df_engine = pd.read_excel('engines/' + engine + '.xlsx')
                break
                
        except IOError:
                print ('No file engines/' + engine + '.xlsx')
                engine = cap_engine(raw_input("Enter engine: "))
        
# lists from excel thrust data               
seconds = list(df_engine.colA)
thrusts = list(df_engine.colB)

#

# enable charting for engine clusters
Engine_quantity = int(prompt("Enter number of engines: "))

# propellant and motor mass in kg for this engine, note must multiply by number of engines
Propellant_mass_kg = Engine_quantity * 0.001 * propellant_mass_gms[engine]
Motor_mass_kg = Engine_quantity * 0.001 * motor_mass_gms[engine]

print "Creating", engine, "thrust curve on", granularity, "second intervals."

new_seconds = []
new_thrusts = []

duration = seconds[len(seconds) - 1]
x = np.arange(0.0, duration, granularity)
new_seconds = list(x)

for i in range(0, len(new_seconds)):
    for j in range(0, len(seconds)):
        if (seconds[j] > new_seconds[i]):
            # here fit curve
            m = (thrusts[j] - thrusts[j-1])/(seconds[j] - seconds[j-1])
            b = thrusts[j] - m * seconds[j]
            new_thrust = m * new_seconds[i] + b
            x = Decimal(new_thrust)
            new_thrusts += [round (x, 2)]
            break
    continue

# true up lists by adding a last 0.0  value for new_thrusts, if needed 
if (len(new_thrusts) < len(new_seconds)):
    new_thrusts += [0.0]

    # SIMULATION SETUP
#
# launch mass and launch weight
empty_mass_gms = rocket_mass_gms + payload_mass_gms + Engine_quantity * (motor_mass_gms[engine])
Empty_mass_kg = 0.001 * empty_mass_gms # total mass, no propellant

launch_mass_gms = empty_mass_gms + Engine_quantity * (propellant_mass_gms[engine])
Launch_mass_kg = 0.001 * launch_mass_gms
Launch_weight = Acc_Grav * Launch_mass_kg

# measure diameter in meters
Dia_meters = dia_mm / 1000.0

print "Computing drag value function based on rocket diameter and air temperature."

# obtain the density of air at Temperature in Centigrade in kg/m^3
Rho = air_density(Temperature)

# compute drag value
Drag_Value = 0.5 * Rho * pi * Coefficient_drag*((Dia_meters/2.0)**2)

print "Computing thrust curve for", Engine_quantity, "engine(s)."

Thrust_curve = [(x * Engine_quantity) for x in new_thrusts]

print "Computing total impulse for rocket."

Thrust_total = sum(Thrust_curve)
Total_impulse = two_digits(granularity * Thrust_total)

print "Computing average thrust and average mass loss in", granularity, "second intervals."

# thrust duration in 0.02 seconds
Thrust_duration = len(Thrust_curve) - 2 # check why -2, probably because bracketed by 0.0

# average thrust in newtons in each granularity interval
Avg_thrust = float(Thrust_total / Thrust_duration)

# average incremental loss of mass, in kg, in granularity interval
Avg_mass_dec = float(Propellant_mass_kg / Thrust_duration)

print "Preparing simulation...."
print
#
# initialize simulation loop variables
#
Acceleration = 0.0
Velocity = 0.0
Altitude = 0.0
interval = 0
Max_velocity = 0.0
time_list = []
alt_list = []

# print header
#
print 'Time\tHeight\tSpeed\tdelta-V\tGees\tMass\tFuel\tfeet\tmph'
#

# SIMULATION MAIN LOOP
#
while((Velocity >= 0.0) or (interval < 10)):
    Gees=Acceleration/Acc_Grav      
    Time=interval * granularity
    Total_mass_kg = Empty_mass_kg + Propellant_mass_kg
    mph = 2.2369363 * Velocity
    feet = 3.2808399 * Altitude
    # lists for rendering graph
    time_list += [two_digits(Time)]
    alt_list += [two_digits(feet)]
    # print values for every 5 intervals (typically every 0.1 second)
    if ((interval % 5) == 0):
        print '%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f'%(Time, Altitude, Velocity, Acceleration, Gees, Total_mass_kg*1000.0, Propellant_mass_kg*1000.0, feet, mph)
    interval +=1
    if (interval<=Thrust_duration):
        Current_thrust = Thrust_curve[interval]
    else:
        Current_thrust = 0
    # check whether we have enough thrust to lift off
    if (Altitude <= 0.0):
        Altitude = 0.0
        liftoff = (Current_thrust >= (Acc_Grav*Total_mass_kg))
        if (liftoff):
            liftoff_time = Time
    if (Altitude > 0.0):
        liftoff = True
    # calculate net force acting on rocket
    NetForce = Current_thrust - (Acc_Grav*Total_mass_kg) - (Drag_Value*(Velocity**2))
    if (liftoff):
        Acceleration = NetForce/Total_mass_kg
        Velocity += Acceleration * granularity
    if (Velocity > Max_velocity):
        Max_velocity = Velocity
    Altitude += Velocity * granularity
    if (interval == Thrust_duration):
        print "Propellant exhaustion"
        # the following for rendering graph
        prop_out_time = two_digits(Time)
        prop_out_alt = two_digits(feet)
        prop_out_point = interval
        
    if (interval >= Thrust_duration):
        Propellant_mass_kg = 0.0
    if (interval < Thrust_duration):
        # propellant mass loss proportional to thrust in interval
        Mass_dec = Avg_mass_dec * (Current_thrust / Avg_thrust)
        Propellant_mass_kg -= Mass_dec
        if (Propellant_mass_kg < 0.0):
            Propellant_mass_kg = 0.0

delay = Time - granularity * Thrust_duration

#
# SIMULATION REPORT SUMMARY
#
print
print "SUMMARY"
print
print "Vehicle mass:", rocket_mass_gms, "grams"
print "Payload mass:", payload_mass_gms, "grams"
print "Launch weight:", two_digits(Launch_weight), "Newton-seconds."
print "Airframe diameter:", dia_mm, "millimeters"
print "Drag coefficient:", Coefficient_drag, "assumed"
print "Mass after propellant exhaustion", rocket_mass_gms + payload_mass_gms + motor_mass_gms[engine], "grams"
if (launch_mass_gms > 1500.0):
    print "Launch mass of", launch_mass_gms, "grams exceeds 1500 grams."
    print "Level 1 certification is required."
print
print 'Launch temperature: %.2f deg. C (%.2f deg. F),'%(Temperature, C_to_F(Temperature))
print
print "Motor used:", engine
print engine, "dry mass:", motor_mass_gms[engine], "gms."
print engine, "propellant mass:", propellant_mass_gms[engine], "gms."
if (Engine_quantity > 1):
    print "Motor cluster using", Engine_quantity, engine, "motors."
print "Total impulse:", Total_impulse, "Newton-seconds."
if (((Engine_quantity == 1) and (Total_impulse >= 160.01)) or (Total_impulse >= 320.01)\
    or (Avg_thrust > 80.0) or (propellant_mass_gms[engine] > 125.0)):
    print "Because of engine thrust/propellant, Level 1 certification is required."
else:
    print
print "Burn time:", two_digits(granularity * Thrust_duration), "seconds."
print
print "Launch weight to total impulse ratio check:",
if (Total_impulse > 5.0 * Launch_weight):
    print "YES"
else:
    print "NO"
if (liftoff_time > 0.0):
    print "liftoff occurs at", liftoff_time, "seconds after ignition."
print
print 'Maximum altitude = %.2f meters, = %.2f feet'%(Altitude, 3.2808399*Altitude)
print 'Maximum velocity = %.2f m/s, = %.2f mi/hr'%(Max_velocity, 2.237*Max_velocity)
print
print 'Recommended delay = %.2f seconds'%(delay)
print
    
while True:
    simcurvep = raw_input("Create simulation curve? (Y) or (N)")

    if (simcurvep in ['Y', 'y', 'N', 'n']):
        break
    else:
        print "Please enter 'Y' or 'N'"

if (simcurvep in ['N', 'n']):
    import sys
    sys.exit()

print "Preparing simulation curve"
trace0 = go.Scatter(
    x = time_list,
    y = alt_list,
    line = dict(
        color = ('rgb(205, 12, 24)'),
        width = 4)
    )

time_delay_list = []
alt_delay_list = []
time_from_out = prop_out_time
for i in range(prop_out_point, interval):
        time_delay_list += [time_from_out]
        time_from_out += granularity
        alt_delay_list += [prop_out_alt]

trace1 = go.Scatter(
    x = time_delay_list,
    y = alt_delay_list,
    line = dict(
        color = ('rgb(22, 96, 167)'),
        width = 4,),
    
    )


data = [trace0, trace1]

# Edit the layout
layout = dict(title = 'Flight: Mass = '+str(empty_mass_gms)+' g, Dia: ' +str(dia_mm) + ' mm, Coef. fr: '+ str(Coefficient_drag)\
              + '<br>using ' + str(Engine_quantity) + ' ' + engine \
              + ' engines(s) at a launch temperature of ' + str(two_digits(Temperature)) + ' deg.C (' + str(C_to_F(Temperature)) + ' deg. F)'\
              + '<br>Apogee: ' + str(two_digits(3.2808399*Altitude)) + ' feet, maximum velocity: ' + str(two_digits(2.237*Max_velocity)) + ' mph', 
                showlegend = False,
              xaxis = dict(title = 'Time in seconds'),
              yaxis = dict(title = 'Altitude in feet'),
              annotations=[
                dict(
                        x=prop_out_time,
                        y=prop_out_alt,
                        xref='x',
                        yref='y',
                        text='propellant exhaustion',
                        showarrow=True,
                        arrowhead=7,
                        ax=0,
                        ay=-40
                ),
                dict(
                        x=prop_out_time + (Time - prop_out_time)/2.0,
                        y=prop_out_alt,
                        xref='x',
                        yref='y',
                        text = 'recommended delay ' + str(delay) + ' seconds.',
                        ax=0,
                        ay = -12
                ),
                dict(
                        x=Time,
                        y=prop_out_alt,
                        xref='x',
                        yref='y',
                        text='ejection',
                        showarrow=True,
                        arrowhead=7,
                        ax=0,
                        ay=-40
                )
            ]
        )

fig = dict(data=data, layout=layout)

plotly.offline.plot(fig, filename = engine + 'rocsim flight curve.html')               

#


