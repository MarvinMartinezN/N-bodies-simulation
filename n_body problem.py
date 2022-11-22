import numpy as np
import plotly.express as px
import pandas as pd
G=6.674e-11
#Set time of simulation and step size (seconds)
tf= 30e6
dt= 3600
Tdata=np.arange(0,tf,dt)
#body functions
BodyList=[]
class Body():
    def __init__(self,name,mass,x0,y0,z0,vx0,vy0,vz0):
        BodyList.append(self)
        self._name=name
        
        self._m=float(mass)
        self._r=np.array([x0,y0,z0],dtype=np.float64)
        self._v=np.array([vx0,vy0,vz0],dtype=np.float64)
        self._a=np.array([0,0,0],dtype=np.float64)

        self._rdata=self._r[np.newaxis,...]
        self._vdata=self._v[np.newaxis,...]
        self._adata=self._a[np.newaxis,...]
        self._Fdata=self._m*self._a
        self._a=np.array([0,0,0],dtype=np.float64)
    def NewCoords(self,ax,ay,az):
        self._a=np.array([float(ax),float(ay),float(az)],dtype=np.float64)
        self._r += self._v*dt
        self._v += self._a*dt

        self._rdata = np.concatenate((self._rdata,self._r[np.newaxis,...]),axis=0)
        self._vdata = np.concatenate((self._vdata,self._v[np.newaxis,...]),axis=0)
        self._adata = np.concatenate((self._adata,self._a[np.newaxis,...]),axis=0)
        self._Fdata = self._m*self._a
    def VelNorm(self,vx,vy,vz):
        self._v += np.array([-vx,-vy,-vz],dtype=np.float64)
        self._vdata += np.array([-vx,-vy,-vz],dtype=np.float64)
#Set Bodies:
Sola = Body("Sola",1.989e30,0,0,0,0,0,0)
Jorden = Body("Jorden",5.972e24,1.487850287972004E+11,-3.337359882077453E+09, 2.934286743328907E+07,2.149316965365740E+02,2.966087724794572E+04,-1.357445099479548E-00)
Merkur= Body("Merkur",3.285e23,5.231712151945569E+10,-5.542160322477117E+09,-5.358361976370921E+09,-4.281551989672854E+03,5.060465996645951E+04,4.529504366225421E+03)
Maned=Body("Maned",7.348e22,1.485104151036709E+11,-3.046923420478556E+09,6.522117052475084E+07,-4.759371607473837E+02, 2.896703732355526E+04,-2.859787360735666E+00)
Mars=Body('Mars',6.4171e23,1.783364318365985E+11,1.182603844695240E+11,-1.902090877974793E+09,-1.238392094305106E+04,2.230290950179835E+04,7.716582809797572E+02)
#Normalize v
ptot=np.array([0,0,0],dtype=np.float64)
for i in BodyList:
    ptot+= i._m*i._v
Mtot = 0
for i in BodyList:
    Mtot+= i._m
vtot=ptot /Mtot
for i in BodyList:
    i.VelNorm(vtot[0],vtot[1],vtot[2])
#Run simulation
for t in Tdata[1:]:
    for i in BodyList:
        i._a = np.array([0,0,0],dtype=np.float64)
        for j in BodyList:
            if i != j:
                dr = j._r-i._r
                ds = np.linalg.norm(dr)
                i._a += (G*j._m/(ds**3))*dr
    for i in BodyList:
        i.NewCoords(i._a[0],i._a[1],i._a[2])
print('Simulation done succesfully, waiting for graphics')
#plotting
range=2.5e11
def dataframe(x):
    body = list([x._name]*len(x._rdata[::100,0]))
    return pd.DataFrame({'x':x._rdata[::100,0],'y':x._rdata[::100,1],'z':x._rdata[::100,2],'body':body})
df=pd.concat([dataframe(i) for i in BodyList])

fig = px.line_3d(df, x='x', y='y', z='z',color='body',range_x=[-range,range],range_y=[-range,range],range_z=[-range,range])
fig.show()