import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import sosh
import sys
from matplotlib import cm
from matplotlib.colors import ListedColormap

nargs=len(sys.argv)
if (nargs > 1):
  for i in range(1, nargs):
    exec(sys.argv[i])

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

arrlist=sosh.arrlist_int
lmaxlist=sosh.lmaxlist_int
varlist={}

interval = 1.0/25
maxsave=0.0

#Animation function to call
def calcylmt(i):

  plotvar=sosh.plotvar
  if plotvar in ['Vr','Vmag','Vsq']:
    ylmt=ylm*np.exp(-1.0j*2*np.pi*float(i) / np.float(nframes))
  if plotvar in ['Vt', 'Vh', 'Vmag','Vsq']:
    dylmtt=dylmt*np.exp(-1.0j*2*np.pi*float(i) / np.float(nframes))
  if plotvar in ['Vt', 'Vh', 'Vmag','Vsq']:
    dylmpt=dylmp*np.exp(-1.0j*2*np.pi*float(i) / np.float(nframes))

  if plotvar in ['Vh', 'Vmag', 'Vsq']:
    vh2=np.square(dylmtt.real) + np.square(dylmpt.real/np.sin(theta))
  if plotvar in ['Vmag', 'Vsq']:
    v2=np.square(fr*ylmt.real) + np.square(fh)*vh2

  if (plotvar == 'Vr'):
    d=fr*ylmt.real
  if (plotvar == 'Vt'):
    d=fh*dylmtt.real
  if (plotvar == 'Vp'):
    d=fh*dylmpt.real/np.sin(theta)
  if (plotvar == 'Vh'):
    d=fh*np.sqrt(vh2)
  if (plotvar == 'Vmag'):
    d=np.sqrt(v2)
  if (plotvar == 'Vsq'):
    d=v2

  return d

def ylmanimate(i):

  d=calcylmt(i)
  global maxsave
  mn,mx=d.min(),d.max()
  maxabs=np.abs([mn,mx]).max()
  if (maxabs > maxsave):
    maxsave=maxabs
#    print("new max=%f"%maxabs, end=" ", flush=True)

  im.set_data(d)
  fig.canvas.draw()
  return

sosh.animate=ylmanimate
sosh.calcimage=calcylmt
sosh.isave=0


print("You may enter 'q' to quit or 'c' to change parameters.")
print("Enter nothing to retain previous values.\n")
while True:

  (l,m,n)=sosh.querylmn(1)
  if l == -1:
    break
  ntest = sosh.modeln[l==sosh.modell]
  while (n not in ntest):
    print("That n was not modelled for that l. Modelled values are in the range %i to %i. Try again." \
          % (ntest.min(),ntest.max()))
    (l,m,n)=sosh.querylmn(1)
    ntest = sosh.modeln[l==sosh.modell]

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

  if (m > 0):
    phi+=(np.pi/4)/m
#    phi[:,int(nx/2):nx]=(np.pi/4)/m
#    phi[:,0:int(nx/2)]=(np.pi/4)/m + np.pi
#  else:
#    phi=0.0*theta
  ylm=y*np.exp(1.0j*signedm*phi)
  dylmt=-np.sin(theta)*dy*np.exp(1.0j*signedm*phi)
  dylmp=1.0j*signedm*y*np.exp(1.0j*signedm*phi)

  xir, xih = sosh.getradial(l,n)
  fr=np.interp(r,sosh.rmesh,xir*np.sqrt(sosh.rho))
  fh=np.interp(r,sosh.rmesh,xih*np.sqrt(sosh.rho))

  vh2=(np.square(dylmt.real) + np.square(dylmp.real/np.sin(theta)))
  v2=(np.square(fr*ylm.real) + np.square(fh)*vh2)

  varlist['Vr']=fr*ylm.real
  varlist['Vt']=fh*dylmt.real
  varlist['Vp']=fh*dylmp.real/np.sin(theta)
  varlist['Vh']=fh*np.sqrt(vh2)
  varlist['Vmag']=np.sqrt(v2)
  varlist['Vsq']=v2

  var=varlist[sosh.plotvar]
  maxsave=0.0

  if (sosh.capflag == 1):
    capblank=r'$\ell$ = {0}' '\n' r'$m$ = {1}' '\n' r'$n$ = {2}'
  elif (sosh.capflag == 2 or sosh.capflag == 3):
    capblank='\n' r'$\ell$ = {0}, $m$ = {1}, $n$ = {2}' '\n'
  else:
    capblank=''
  caption=capblank.format(l,signedm,n)

  fig,im = sosh.drawfigure(var, caption=caption, fsize=figsize)
  if (sosh.isave != 0):
    sosh.savefigure(ianimate=animate, label='_interior')
  if (show != 0):
    if (animate != 0):
      anim = animation.FuncAnimation(fig, ylmanimate, frames=nframes, interval=interval)
    plt.show()
# end of input loop
  print()

# stuff that didn't work
#writer = animation.writers['ffmpeg']
#writer = animation.FFMpegWriter(fps=15, bitrate=1800)
#  anim.save('test.mp4', writer=writer) #, dpi=600, fps = 1/interval, writer='ffmpeg')


