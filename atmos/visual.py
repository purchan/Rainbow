import numpy  as np
import math   as m
from   vpython import *

# WIP
## in-angle out-angle
## h - E

# Ref
def hToTemperature(h):      # 0     < h < 85000
    if h <= 0:
        return 288.15
    elif h < 11000:
        return 288.15 - 6.5/1000*h
    elif h < 20000:
        return 216.65
    elif h < 32000:
        return 216.65 + 1/1000*h
    elif h < 47000:
        return 228.65 + 2.8/1000*h
    elif h < 51000:
        return 270.65
    elif h < 71000:
        return 270.65 - 2.8/1000*h
    elif h < 85000:
        return 214.65 - 2/1000*h
    else: return 200
def hToPressure(h):         # Let P(0)=1, unit is not important
    sum = 0
    for x in np.arange(0,h+100,100):
        sum += hToTemperature(x)
    sum /= (h/100)
    
    return m.exp(-9.8*h/287/sum)
def hToN(h):
    if h >= 85000:
        return 1
    else:
        return 1+2.98E-4*hToPressure(h)*hToTemperature(0)/hToTemperature(h)

# Config
## Simul Env
posE  = vec(0, 0, 0)
radE  = 6371000

scene = canvas(width=800, height=800, center=vec(0, radE, 0), background=vec(0.25, 0.23, 0.20))
dt    = 0.001
rt    = 1000
## Gen
genL  = [vec(-100, radE+85100, 0)]
lisL  = []

# lisR  = []
# def gatherResult(o, _str):
#     lisR.append({"PS":_str, "RL": o.rl, "WL": o.λ, "E ":o.E, "Ip":o.initpos, "Cp":o.pos})
# def printResult():
#     for res in lisR:
#         print('{')
#         for key in res:
#             print('    "',key,'" :',res[key])
#         print('},')

class LP(box):
    def __init__(self, pos, v=vec(100, -100, 0), E=1):
        box.__init__(self, pos=pos, axis=vec(0,1,0), size=vec(1, 1, 1), color=color.white, make_trail=True, trail_radius=1000)
        self.v  = v
        self.E  = E
    def nextpos(self):
        global dt, posE, radE
        # change pos->n->v
        self.nc   = hToN(mag(self.pos-posE)-radE)
        self.pos += self.v * dt
        
        hc= mag(self.pos-posE)-radE
        if hc < 0:
            return False
        self.n    = hToN(hc)
        
        self.fraction(self.nc, self.n)
        return True
    def fraction(self, n1, n2):
        global posE
        θi = diff_angle(posE-self.pos,self.v)
        if θi > m.pi/2: θi = m.pi - θi
        Lproj1 = proj(self.v, posE-self.pos)
        Pproj1 = self.v - Lproj1
        try:
            θ2 = asin(sin(θi)*n1/n2)
        except:
            self.v = Pproj1-Lproj1
        else:
            # Fresnel equation, sunlight is classify as unpolarized
            Rs = ((n1*m.cos(θi) - n2*m.cos(θ2))/(n1*m.cos(θi) + n2*m.cos(θ2)))**2
            Rp = ((n1*m.cos(θ2) - n2*m.cos(θi))/(n1*m.cos(θ2) + n2*m.cos(θi)))**2
            Erl = (Rs+Rp)/2
            
            Lproj2 = mag(self.v)*cos(θ2)*norm(Lproj1)
            Pproj2 = mag(self.v)*sin(θ2)*norm(Pproj1)
            
            self.v  = Pproj2+Lproj2
            self.E *= Erl
                
for pos in genL:       # generateLight
    lisL.append(LP(pos))

earth = sphere(radius=radE, color=vec(1, 0.7, 0.2))
atmos = sphere(radius=radE + 85000, color=color.white,opacity=0.05)
earth.r = earth.radius
earth.pos = posE
atmos.pos = posE

while len(lisL) > 0:  # main loop
    rate(rt)
    for ptcl in lisL:
        alive = ptcl.nextpos()
        if mag(ptcl.pos) > radE + 100000 or not alive:
            lisL.remove(ptcl)

print('bye')
# printResult()