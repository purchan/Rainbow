import numpy  as np
import math   as m
import functools
from vpython import *

# WIP
## camera rotation(?)
## energy instead of random

# Reference
def Z(T,p,xw): # compressibility
    t=T-273.15
    a0 = 1.58123e-6   #K·Pa^-1
    a1 = -2.9331e-8   #Pa^-1
    a2 = 1.1043e-10   #K^-1·Pa^-1
    b0 = 5.707e-6     #K·Pa^-1
    b1 = -2.051e-8    #Pa^-1
    c0 = 1.9898e-4    #K·Pa^-1
    c1 = -2.376e-6    #Pa^-1
    d  = 1.83e-11     #K^2·Pa^-2
    e  = -0.765e-8    #K^2·Pa^-2
    return 1-(p/T)*(a0+a1*t+a2*t**2+(b0+b1*t)*xw+(c0+c1*t)*xw**2) + (p/T)**2*(d+e*xw**2)
def n(λ,t,p,h,xc):
    # λ: wavelength, 0.3 to 1.69 μm 
    # t: temperature, -40 to +100 °C
    # p: pressure, 80000 to 120000 Pa
    # h: fractional humidity, 0 to 1
    # xc: CO2 concentration, 0 to 2000 ppm

    σ = 1/λ           #μm^-1
    
    T= t + 273.15     #Temperature °C -> K
    
    R = 8.314510      #gas constant, J/(mol·K)
    
    k0 = 238.0185     #μm^-2
    k1 = 5792105      #μm^-2
    k2 = 57.362       #μm^-2
    k3 = 167917       #μm^-2
 
    w0 = 295.235      #μm^-2
    w1 = 2.6422       #μm^-2
    w2 = -0.032380    #μm^-4
    w3 = 0.004028     #μm^-6
    
    A = 1.2378847e-5  #K^-2
    B = -1.9121316e-2 #K^-1
    C = 33.93711047
    D = -6.3431645e3  #K
    
    α = 1.00062
    β = 3.14e-8       #Pa^-1,
    γ = 5.6e-7        #°C^-2

    #saturation vapor pressure of water vapor in air at temperature T
    if(t>=0):
        svp = np.exp(A*T**2 + B*T + C + D/T) #Pa
    else:
        svp = 10**(-2663.5/T+12.537)
    
    #enhancement factor of water vapor in air
    f = α + β*p + γ*t**2
    
    #molar fraction of water vapor in moist air
    xw = f*h*svp/p
    
    #refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, 450 ppm CO2
    nas = 1 + (k1/(k0-σ**2)+k3/(k2-σ**2))*1e-8
    
    #refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, xc ppm CO2
    naxs = 1 + (nas-1) * (1+0.534e-6*(xc-450))
    
    #refractive index of water vapor at standard conditions (20 °C, 1333 Pa)
    nws = 1 + 1.022*(w0+w1*σ**2+w2*σ**4+w3*σ**6)*1e-8
    
    Ma = 1e-3*(28.9635 + 12.011e-6*(xc-400)) #molar mass of dry air, kg/mol
    Mw = 0.018015                            #molar mass of water vapor, kg/mol
    
    Za = Z(288.15, 101325, 0)                #compressibility of dry air
    Zw = Z(293.15, 1333, 1)                  #compressibility of pure water vapor
    
    #Eq.4 with (T,P,xw) = (288.15, 101325, 0)
    ρaxs = 101325*Ma/(Za*R*288.15)           #density of standard air
    
    #Eq 4 with (T,P,xw) = (293.15, 1333, 1)
    ρws  = 1333*Mw/(Zw*R*293.15)             #density of standard water vapor
    
    # two parts of Eq.4: ρ=ρa+ρw
    ρa   = p*Ma/(Z(T,p,xw)*R*T)*(1-xw)       #density of the dry component of the moist air    
    ρw   = p*Mw/(Z(T,p,xw)*R*T)*xw           #density of the water vapor component
    
    nprop = 1 + (ρa/ρaxs)*(naxs-1) + (ρw/ρws)*(nws-1)
    
    return nprop     # 0.3   < λ < 1.69
def lambdaToRgb(λ):
    λ *= 1000
    gamma = 0.80
    intensityMax = 1
    if   λ>=380 and λ<440:
        r = -(λ - 440) / 60
        g = 0.0
        b = 1.0
    elif λ>=440 and λ<490:
        r = 0.0
        g =  (λ - 440) / 50
        b = 1.0
    elif λ>=490 and λ<510:
        r = 0.0
        g = 1.0
        b = -(λ - 510) / 20
    elif λ>=510 and λ<580:
        r =  (λ - 510) / 70
        g = 1.0
        b = 0.0
    elif λ>=580 and λ<645:
        r = 1.0
        g = -(λ - 645)/ 65
        b = 0.0
    elif λ>=645 and λ<781:
        r = 1.0
        g = 0.0
        b = 0.0
    else:
        r = 0.0
        g = 0.0
        b = 0.0
        
    # Let the intensity fall off near the vision limits
    if   λ>=380 and λ<420:
        factor = 0.3 + 0.7*(λ - 380) / 40
    elif λ>=420 and λ<701:
        factor = 1.0
    elif λ>=701 and λ<781:
        factor = 0.3 + 0.7*(780 - λ) / 80
    else:
        factor = 0.0
    
    if not r==0:
        r = intensityMax * pow(r*factor, gamma)
    if not g==0:
        g = intensityMax * pow(g*factor, gamma)
    if not b==0:
        b = intensityMax * pow(b*factor, gamma)
    
    return vec(r, g, b)        # 0.38  < λ < 0.78
def lambdaToNW(λ):
    a  = 0.5672526103
    a1 = 0.1736581125
    a2 = 0.021211531502
    a3 = 0.1138493213
    
    b  = 0.005085550461
    b1 = 0.01814938654
    b2 = 0.02617260739
    b3 = 10.73888649
     
    return sqrt(1 + a/(1-b/λ**2) + a1/(1-b1/λ**2) + a2/(1-b2/λ**2) + a3/(1-b3/λ**2)) # 0.182 < λ < 1.129

lambdaToN   = functools.partial(n, t=19, p=101325, h=0, xc=400)  # check n for specific information
# Configurable Constants
## Simulation
background_color = vec(1, 0.92, 0.80)
dt     = 0.001   # time step
rt     = 1000    # max rate

light_r = 0.0001
light_v = vec(10, 0, 0)
light_maketr = True
light_trail  = 0.15

drop_r  = 50

## Generator
light_gen = [vec(-60, -40, x) for x in range(0, 1, 1)]
drop_gen  = [vec(0, 0, 0)]

light_list = [[]]
drop_list  = []

# functions
class LP(sphere):     # Light Particle
    def __init__(self, λ, pos, _rl, _v=light_v, _e=1):
        global light_r, light_v, light_maketr, light_trail
        self.λ = λ
        sphere.__init__(self, radius=light_r, color=lambdaToRgb(λ), make_trail=light_maketr, trail_radius=light_trail)
        
        self.n = lambdaToN(λ)
        self.nw = lambdaToNW(λ)
        
        self.pos = pos
        self.v = _v
        
        # assert _rl == "r" or "l"
        self.E   = _e         # Energy
        self._rl = _rl
    def nextpos(self):
        global dt
        return self.pos + self.v * dt
    def tryRefractReflect(self):
        global drop_list
        for drop in drop_list:
            dis_dl = drop.pos - self.pos
            # res[1] == hasRefract()
            if   mag(dis_dl)>drop.r and mag(drop.pos-self.nextpos())<drop.r:
                res = self.fraction("a2w", drop)
                self.v = res[0]
                while mag(dis_dl)>drop.r and res[1]:
                    self.pos = self.nextpos()
                    dis_dl = drop.pos - self.pos
            elif mag(dis_dl)<drop.r and mag(drop.pos-self.nextpos())>drop.r:
                res = self.fraction("w2a", drop)
                self.v = res[0]
                while mag(dis_dl)<drop.r and res[1]:
                    self.pos = self.nextpos()
                    dis_dl = drop.pos - self.pos
    def fraction(self, _str, drop):
        global light_list
        assert _str == "a2w" or _str == "w2a"
        if _str == "a2w":
            n1 = self.n
            n2 = self.nw
        elif _str == "w2a":
            n1 = self.nw
            n2 = self.n
        θi = diff_angle(drop.pos-self.pos,self.v)
        if θi > m.pi/2:
            θi = m.pi - θi
        linear_proj1 = proj(self.v,drop.pos-self.pos)
        planar_proj1 = self.v-linear_proj1
        try:
            θ2 = asin(sin(θi)*n1/n2)
        except:
            return (planar_proj1-linear_proj1, False)
        else:
            # Fresnel equation, sunlight is classify as unpolarized
            Rs = ((n1*m.cos(θi) - n2*m.cos(θ2))/(n1*m.cos(θi) + n2*m.cos(θ2)))**2
            Rp = ((n1*m.cos(θ2) - n2*m.cos(θi))/(n1*m.cos(θ2) + n2*m.cos(θi)))**2
            rl = (Rs+Rp)/2
            
            linear_proj2 = mag(self.v)*cos(θ2)*norm(linear_proj1)
            planar_proj2 = mag(self.v)*sin(θ2)*norm(planar_proj1)
            if self._rl=="r":     # Refraction
                self.E *= rl
                
                lp_l = LP(self.λ, self.pos, "l", planar_proj1-linear_proj1, self.E*(1-rl))
                light_list[0].append(lp_l)
                
                return (linear_proj2+planar_proj2, True)
            if self._rl=="l":     # Reflection
                self.E *= 1-rl
                
                lp_r = LP(self.λ, self.pos, "r", planar_proj2+linear_proj2, self.E*(rl))
                while mag(drop.pos - lp_r.pos)<drop.r:
                    lp_r.pos = lp_r.nextpos()
                
                light_list[0].append(lp_r)
                
                return (planar_proj1-linear_proj1, False)

# main code
## Set-up scene
scene = canvas(width=800, height=800, center=vec(0, 0, 0), background=background_color)

for pos in light_gen:       # generateLight
    sublist = []
    for λ in np.arange(0.38, 0.78, 0.01):
        tmpr = LP(λ, pos, "r")
        sublist.append(tmpr)
    light_list.append(sublist)

drop_furthest = mag(drop_gen[0])
for pos in drop_gen:        # generateDrop
    # assert every drop does not collide with each other
    tmp   = sphere(radius=drop_r, color=color.white,opacity=0.1)
    tmp.r = tmp.radius
    drop_list.append(tmp)
    
    if mag(pos) > drop_furthest:
        drop_furthest = mag(pos)

while len(light_list) > 0:  # main loop
    rate(rt)
    for sublist in light_list:
        if (len(sublist)==0 and light_list.index(sublist)!=0) or (len(light_list)==1 and len(light_list[0])==0):
            light_list.remove(sublist)
        for ptcl in sublist:
            ptcl.tryRefractReflect()
            ptcl.pos = ptcl.nextpos()
            ptcl.opacity = 0
            if mag(ptcl.pos) > (drop_furthest+2*drop_r):
                sublist.remove(ptcl)
