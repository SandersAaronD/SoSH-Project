import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import sosh
import sys
import os

nargs=len(sys.argv)
#print(sys.argv)
if (nargs > 1):
  for i in range(1, nargs):
    exec(sys.argv[i])

try:
  freqscale
except:
  freqscale=1.0/3

try:
  bothm
except:
  bothm=False

try:
  dpi
except:
  dpi = 300

try:
  nframes
except:
  nframes = 64

try:
  colorshift
except:
  colorshift=1

try:
  animate
except:
  animate=1

try:
  show
except:
  show=1

try:
  caption
except:
  caption=0

try:
  rsurf
except:
  rsurf=1.0

try:
  pixels
except:
  pixels=1000

try:
  distobs
except:
  distobs=220.0

try:
  xoffset
except:
  xoffset=0.0

try:
  yoffset
except:
  yoffset=0.0

try:
  pangle
except:
  pangle=0.0

try:
  figsize
except:
  figsize=5

sosh.nframes=nframes
sosh.dpi=dpi
sosh.icolshift=colorshift
sosh.capflag=caption

sosh.loadmodel()

(r, theta, phi) = sosh.image2rtheta(xpixels=pixels,ypixels=pixels,distobs=distobs,pangle=pangle,xoffset=xoffset,yoffset=yoffset)
x=np.cos(theta)
(nx,ny)=x.shape
if (rsurf < 1.0):
  newind = (r <= rsurf)
  r.mask = np.logical_not(newind)
  x.mask = np.logical_not(newind)

#modeparms=np.loadtxt('mdi.average.modes')
# use model values instead of measured values
#modeparms=np.loadtxt('model.surface.modes')
#lmod=modeparms[:,0]
#nmod=modeparms[:,1]

arrlist=sosh.arrlist
lmaxlist=sosh.lmaxlist
varlist={}

rlist=[]
tlist=[]
plist=[]
freqlist=[]
amplist=[]
rclist=[]
hclist=[]

interval = 1.0/25
maxsave=0.0

#Animation function to call
def calcsumt(i):

  sumr=0.0
  sumt=0.0
  sump=0.0
  plotvar=sosh.plotvar

  if plotvar in ['Vr','Vmag','Vsq']:
    for k in range(nmodes):
      sumr += freqlist[k]*rclist[k]*rlist[k]*np.exp(-1.0j*2*np.pi*freqlist[k]*freqscale*float(i) / np.float(nframes))
  if plotvar in ['Vt', 'Vh', 'Vmag','Vsq']:
    for k in range(nmodes):
      sumt += freqlist[k]*hclist[k]*tlist[k]*np.exp(-1.0j*2*np.pi*freqlist[k]*freqscale*float(i) / np.float(nframes))
  if plotvar in ['Vt', 'Vh', 'Vmag','Vsq']:
    for k in range(nmodes):
      sump += freqlist[k]*hclist[k]*plist[k]*np.exp(-1.0j*2*np.pi*freqlist[k]*freqscale*float(i) / np.float(nframes))

  if plotvar in ['Vh', 'Vmag', 'Vsq']:
    vh2=np.square(sumt.real) + np.square(sump.real)
  if plotvar in ['Vmag', 'Vsq']:
    v2=np.square(sumr.real) + vh2

  if (plotvar == 'Vr'):
    d=sumr.real
  if (plotvar == 'Vt'):
    d=sumt.real
  if (plotvar == 'Vp'):
    d=sump.real
  if (plotvar == 'Vh'):
    d=np.sqrt(vh2)
  if (plotvar == 'Vmag'):
    d=np.sqrt(v2)
  if (plotvar == 'Vsq'):
    d=v2

  return d

def sumanimate(i):

  d=calcsumt(i)
  global maxsave
  mn,mx=d.min(),d.max()
  maxabs=np.abs([mn,mx]).max()
  if (maxabs > maxsave):
    maxsave=maxabs
#    print("new max=%f"%maxabs, end=" ", flush=True)

  im.set_data(d)
  fig.canvas.draw()
  return maxabs

sosh.animate=sumanimate
sosh.calcimage=calcsumt
sosh.isave=0


print("You may enter 'q' to quit or 'c' to change parameters.")
instr=sosh.catchc("Enter number of modes: ",2)
if (instr == 'q'):
  sys.exit()
nmodes=int(instr)
count=0
lstr=''
mstr=''
nstr=''
pm='\u00B1'
cap3=''
cap3b=r'$\ell$ = {0}, $m$ = {1}, $n$ = {2}' '\n'
for i in range(abs(nmodes)):
  print("MODE #%i" % (i+1))

  (l,m,n)=sosh.querylmn(nmodes)
  if l == -1:
    sys.exit()
  ntest = sosh.modeln[l==sosh.modell]
  while (n not in ntest):
    print("That n was not modelled for that l. Modelled values are in the range %i to %i. Try again." \
          % (ntest.min(),ntest.max()))
    (l,m,n)=sosh.querylmn(nmodes)
    ntest = sosh.modeln[l==sosh.modell]

  xir, xih = sosh.getradial(l,n)
  fr=np.interp(r,sosh.rmesh,xir*np.sqrt(sosh.rho))
  fh=np.interp(r,sosh.rmesh,xih*np.sqrt(sosh.rho))

  idx = ((sosh.modeln == n) & (sosh.modell == l))
  freq = sosh.modelnu[idx]*1000.0
  rclist.append(fr)
  hclist.append(fh)
  freqlist.append(freq)
  print("Scaled frequency = %f." % (freq*freqscale))

  signedm=m
  m=np.abs(m)
  if m not in lmaxlist.keys() or l > lmaxlist[m]:
    lmaxlist[m]=l
    plm=np.ma.array(np.zeros((nx,ny,l-m+1)))
    dplm=np.ma.array(np.zeros((nx,ny,l-m+1)))
    (y,dy)=sosh.setplm(l,m,x,plm,dplm)
    arrlist[m]=(plm,dplm)
    print("Spherical harmonic calculated and saved.")
  else:
    (plm,dplm)=arrlist[m]
    (y,dy)=(plm[...,l-m],dplm[...,l-m])
    print("Spherical harmonic retrieved.")

  if (i == 0):
    if (m > 0):
      phi+=(np.pi/4)/m
#      phi[:,int(nx/2):nx]=(np.pi/4)/m
#      phi[:,0:int(nx/2)]=(np.pi/4)/m + np.pi
#    else:
#      phi=0.0*theta
  ylm=y*np.exp(1.0j*signedm*phi)
  dylmt=-np.sin(theta)*dy*np.exp(1.0j*signedm*phi)
  dylmp=1.0j*signedm*y*np.exp(1.0j*signedm*phi)

  rlist.append(ylm)
  tlist.append(dylmt)
  plist.append(dylmp/np.sin(theta))

  count += 1

  if (bothm and m != 0):
    signedm *= -1
    ylm=y*np.exp(1.0j*signedm*phi)
    dylmt=-np.sin(theta)*dy*np.exp(1.0j*signedm*phi)
    dylmp=1.0j*signedm*y*np.exp(1.0j*signedm*phi)
    rlist.append(ylm)
    tlist.append(dylmt)
    plist.append(dylmp/np.sin(theta))
    rclist.append(rc)
    hclist.append(hc)
    freqlist.append(freq)
    count += 1
    lstr = lstr + str(l) + ', '
    mstr = mstr + pm + str(m) + ', '
    nstr = nstr + str(n) + ', '
    cap3 = cap3 + cap3b.format(str(l),pm+str(m),str(n))
  else:
    lstr = lstr + str(l) + ', '
    mstr = mstr + str(signedm) + ', '
    nstr = nstr + str(n) + ', '
    cap3 = cap3 + cap3b.format(l,signedm,n)


# if using measured mode parameters, this would 
# allow you to compute frequency as a function of m
# from the fitted a-coefficients
#  ai=np.append([0.0],modeparms[ind,12:18]/1000)
#  pols=sosh.apols(l,6)
#  fx=np.matmul(pols,ai)
#  freq=modeparms[ind,2]+fx
#  freqlist.append(freq[l+signedm]*freqscale)
#  amplist.append(modeparms[ind,3])

lstr=lstr.strip(', ')
mstr=mstr.strip(', ')
nstr=nstr.strip(', ')
cap3=cap3.strip('\n')

sumr=0.0
sumt=0.0
sump=0.0
for i in range(count):
  sumr += rclist[i]*rlist[i]*freqlist[i]
  sumt += hclist[i]*tlist[i]*freqlist[i]
  sump += hclist[i]*plist[i]*freqlist[i]

vh2=(np.square(sumt.real) + np.square(sump.real))
v2=(np.square(sumr.real) + vh2)

varlist['Vr']=sumr.real
varlist['Vt']=sumt.real
varlist['Vp']=sump.real
varlist['Vh']=np.sqrt(vh2)
varlist['Vmag']=np.sqrt(v2)
varlist['Vsq']=v2

while True:

  var=varlist[sosh.plotvar]
  maxsave=0.0

  if (sosh.capflag == 1):
    capblank=r'$\ell$ = {0}' '\n' r'$m$ = {1}' '\n' r'$n$ = {2}'
  elif (sosh.capflag == 2):
    capblank='\n' r'$\ell$ = {0};  $m$ = {1};  $n$ = {2}' '\n'
  elif (sosh.capflag == 3):
    capblank=cap3
  else:
    capblank=''
  caption=capblank.format(lstr,mstr,nstr)

  fig,im = sosh.drawfigure(var, fsize=figsize)
  if (sosh.isave != 0):
    sosh.savefigure(ianimate=animate, label='_interiorsum')
  if (show != 0):
    if (animate != 0):
      anim = animation.FuncAnimation(fig, sumanimate, frames=nframes, interval=interval)
    plt.show()

  print()
  c=sosh.queryplotparms()
  if (c == 'q'):
    break
  fs=sosh.catchc("Enter new freqscale: ", freqscale)
  if (fs == 'q'):
    break
  else:
    freqscale=float(fs)
  nf=sosh.catchc("Enter new nframes: ", nframes)
  if (nf == 'q'):
    break
  else:
    nframes=int(nf)
    sosh.nframes=nframes



