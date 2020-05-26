import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import sosh
import sys

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
  pixels
except:
  pixels=1000

try:
  bangle
except:
  bangle=30.0

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

#try:
#  rsurf
#except:
#  rsurf=1.0

sosh.nframes=nframes
sosh.dpi=dpi
sosh.icolshift=colorshift
sosh.capflag=caption

# not needed if using table of surface values
#sosh.loadmodel()
modeparms=np.loadtxt('model.surface.modes')
lmod=modeparms[:,0]
nmod=modeparms[:,1]

(phi,theta)=sosh.image2sphere(xpixels=pixels,ypixels=pixels,bangle=bangle,distobs=distobs,pangle=pangle,xoffset=xoffset,yoffset=yoffset)
x=np.cos(theta)
(nx,ny)=x.shape

arrlist=sosh.arrlist_surf
lmaxlist=sosh.lmaxlist_surf
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
    v2=np.square(rc*ylmt.real) + vh2*hc**2

  if (plotvar == 'Vr'):
    d=ylmt.real
  if (plotvar == 'Vt'):
    d=dylmtt.real
  if (plotvar == 'Vp'):
    d=dylmpt.real/np.sin(theta)
  if (plotvar == 'Vh'):
    d=np.sqrt(vh2)
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
  ntest = nmod[l==lmod]
  while (n not in ntest):
    print("That n was not modelled for that l. Modelled values are in the range %i to %i. Try again." \
          % (ntest.min(),ntest.max()))
    (l,m,n)=sosh.querylmn(1)
    ntest = nmod[l==lmod]

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

  ylm=y*np.exp(1.0j*signedm*phi)
  dylmt=-np.sin(theta)*dy*np.exp(1.0j*signedm*phi)
  dylmp=1.0j*signedm*y*np.exp(1.0j*signedm*phi)

# not needed if using table of surface values
#  xir, xih = sosh.getradial(l,n)
#  indsurf = np.abs(sosh.rmesh-rsurf).argmin()
#
#  rc=float(xir[indsurf])/1e8
#  hc=float(xih[indsurf])/1e8
#  rc=10.0
#  hc=10.0

  ind = ((l==lmod) & (n==nmod))
  freq = float(modeparms[ind,2])/1000.0
#  rc = modeparms[ind,5]/1e8
#  hc = modeparms[ind,6]/1e8
  rat1=modeparms[ind,5]/modeparms[ind,6]
  rc=100.0
  hc=rc/float(modeparms[ind,4])
#  print("Using theoretical vertical to horizontal ratio of %f, actual is %f." % (rc/hc, float(rat1)))
  print("Vertical to horizontal ratio = %f." % (rc/hc))

  vh2=(np.square(dylmt.real) + np.square(dylmp.real/np.sin(theta)))
  v2=(np.square(rc*ylm.real) + vh2*hc**2)

  varlist['Vr']=ylm.real
  varlist['Vt']=dylmt.real
  varlist['Vp']=dylmp.real/np.sin(theta)
  varlist['Vh']=np.sqrt(vh2)
  varlist['Vmag']=np.sqrt(v2)
  varlist['Vsq']=v2

  var=varlist[sosh.plotvar]
  maxsave=0.0

  if (sosh.capflag == 1):
    capblank=r'$\ell$ = {0}' '\n' r'$m$ = {1}' '\n'
  elif (sosh.capflag == 2 or sosh.capflag == 3):
    capblank='\n' r'$\ell$ = {0}, $m$ = {1}' '\n'
  else:
    capblank=''
  caption=capblank.format(l,signedm)

  fig,im = sosh.drawfigure(var, caption=caption, fsize=figsize)
#  plt.text(0.5,-0.1,sosh.caption, horizontalalignment='center', verticalalignment='bottom',transform=ax.transAxes)
  if (sosh.isave != 0):
    sosh.savefigure(ianimate=animate, label='_surface')
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


