import numpy as np
import os
import h5py
import matplotlib.pyplot as plt
#from matplotlib import animation
from matplotlib import cm
from matplotlib.colors import ListedColormap

modeldir='../mods'
modeln=np.array([])
modell=np.array([])
modelnu=np.array([])
rmesh=np.array([])
rho=np.array([])
R=0.0

arrlist_surf={}
lmaxlist_surf={}
arrlist_int={}
lmaxlist_int={}

def loadmodel():

#  in_dir = '/home/tplarson/solar/mods'
  in_dir = modeldir
  fname = 'mods_p3_eigfcn.h5'
  hf = h5py.File(os.path.join(in_dir, fname))

# load r, rho and R from "model"
  global rmesh, rho, R
  rmesh, rho, R = [hf['model'][k].value for k in ['r', 'rho', 'R']]
  rmesh /= R

# load l and n from "modes"
  global modeln, modell, modelnu
  modeln, modell, modelnu = [hf['modes'][k].value for k in ['n', 'l', 'nu']]


def getradial(l,n):

#  in_dir = '/home/tplarson/solar/mods'
  in_dir = modeldir
  fname = 'mods_p3_eigfcn.h5'
  hf = h5py.File(os.path.join(in_dir, fname))

# open y1 and y2 data sets from "modes" (without reading the data)
  y1ds = hf['modes/y1'] 
  y2ds = hf['modes/y2']

  idx = ((modeln == n) & (modell == l))
  y1, y2 = y1ds[idx,:], y2ds[idx,:]

# compute xi_r and xi_h
  L2 = l * (l + 1)
  xir = (y1 * R).flatten()
  if (l > 0):
    xih = (y2 * R / L2).flatten()
  else:
    xih = 0.0*xir

  return (xir, xih)


def writesurfacemodel(lmin=0, lmax=300, rsurf=1.0):

  outfile='model.surface.modes'
  file=open(outfile,'w')
  in_dir = modeldir
  fname = 'mods_p3_eigfcn.h5'
  hf = h5py.File(os.path.join(in_dir, fname))

# load r, rho and R from "model"
  rmesh, c, R = [hf['model'][k].value for k in ['r', 'c', 'R']]
  indsurf = np.abs(rmesh-rsurf*R).argmin()

# load l and n from "modes"
  modeln, modell, modelnu, modelsig2, modelE = [hf['modes'][k].value for k in ['n', 'l', 'nu', 'sigma2', 'E']]

  outfmt='{:d} {:d} {:f} {:f} {:f} {:f} {:f} {:f} {:f}'
  l=lmin
  while (l <= lmax):
    ind = (modell == l)
    nlist = modeln[ind]
    for n in nlist:
      xir, xih = getradial(l,n)
#      rc=float(xir[indsurf])
#      hc=float(xih[indsurf])
      rc=np.interp(rsurf*R,rmesh,xir)
      hc=np.interp(rsurf*R,rmesh,xih)
      ind2 = (modell == l) & (modeln == n)
      nu = float(modelnu[ind2])
      erg = float(modelE[ind2])*1e6
      sig2 = float(modelsig2[ind2])
      if (l > 0):
        L=np.sqrt(l*(l+1))
#        indturn=np.abs(rmesh/c - L/(2*np.pi*nu)).argmin()
#        rturn=float(rmesh[indturn]/R)
        rturn=np.interp(L/(2*np.pi*nu),rmesh/c,rmesh)/R
      else:
        rturn=0.0
      amp=nu*np.sqrt(np.square(rc)+l*(l+1)*np.square(hc))
      nu*=1e6
      amp=float(amp)/30e6
#      outstr=f'{l:d} {n:d} {nu:f} {amp:f} {sig2:f} {rc:f} {hc:f} {rturn:f} {erg:f}'
      outstr=outfmt.format(int(l),int(n),nu,amp,sig2,rc,hc,rturn,erg)
      file.write('%s\n' % (outstr))
    l+=1
  file.close()



def image2sphere(xpixels=1000, ypixels=1000, distobs=220.0, \
                 bangle=30.0, pangle=0.0, xoffset=0.0, yoffset=0.0):

  p = pangle*np.pi/180
  b0 = bangle*np.pi/180
  #distobs = 220.0  # solar radii
  #xpixels = 1000
  #ypixels = 1000 
  x0 = xpixels/2-0.5 + xoffset
  y0 = ypixels/2-0.5 + yoffset
  imscale = 1.97784*1024/xpixels
  scale = imscale/(180*60*60/np.pi)
  rsun=np.tan(np.arcsin(1/distobs))/scale

# Sun is sitting at the center of the main coordinate system and has radius 1.
# Observer is at robs=(xobs,yobs,zobs) moving with a velocity vobs.
# Start by setting robs from distobs and b0

  robs_x = distobs * np.cos(b0)
  robs_y = 0.0
  robs_z = distobs * np.sin(b0)

# image coordinates
  xx = np.linspace(0, xpixels-1, xpixels)
  yy = np.linspace(0, ypixels-1, ypixels)
  xx, yy = np.meshgrid(xx, yy)

  x2 = scale*(xx - x0)
  y2 = scale*(yy - y0)
# Rotate by the P-angle. New coordinate system has the y-axis pointing
# towards the solar north pole.
  x1 =  x2*np.cos(p) + y2*np.sin(p)
  y1 = -x2*np.sin(p) + y2*np.cos(p)

# Now transform to put the coordinates into the solar coordinate system.
# First find the directions (vecx and vecy) the x2 and y2 coordinate
# axis correspond to. vecsun points towards the Sun. Note that the
# (x2,y2,Sun) system is left handed. These vectors are unit vectors.

  vecx_x=0.0
  vecx_y=1.0
  vecx_z=0.0
  vecy_x=-np.sin(b0)
  vecy_y=0.0
  vecy_z=np.cos(b0)
  vecsun_x=-np.cos(b0)
  vecsun_y=0.0
  vecsun_z=-np.sin(b0)

# Now the proper direction can be found. These are not unit vectors.
  x = vecx_x*x1 + vecy_x*y1 + vecsun_x
  y = vecx_y*x1 + vecy_y*y1 + vecsun_y
  z = vecx_z*x1 + vecy_z*y1 + vecsun_z
  qq = 1/np.sqrt(x*x + y*y + z*z)
# Make them unit vectors.
  x*=qq
  y*=qq
  z*=qq

# Now find intersection with the Sun.
# Solve quadratic equation |robs+q*[x1,y1,z1]|=1 for q
# a, b and c are terms in a*x^2+bx+c=0. a==1 since [x1,y1,z1] is unit vector.
  c = robs_x*robs_x + robs_y*robs_y + robs_z*robs_z -1
  b = 2*(x*robs_x+y*robs_y+z*robs_z)
  d = b*b - 4*c
  index = (d >= 0)
  q = np.zeros([xpixels,ypixels])
  xsun = np.zeros([xpixels,ypixels])
  ysun = np.zeros([xpixels,ypixels])
  zsun = np.zeros([xpixels,ypixels])
  q[index]=(-b[index] - np.sqrt(d[index]))/2
  xsun[index]=robs_x + x[index]*q[index]
  ysun[index]=robs_y + y[index]*q[index]
  zsun[index]=robs_z + z[index]*q[index]

  phisun = np.arctan2(ysun,xsun)
  thetasun = np.pi/2 - np.arcsin(zsun)

  ph=np.ma.array(phisun, mask=np.logical_not(index))
  th=np.ma.array(thetasun, mask=np.logical_not(index))

  return (ph, th) 


def image2rtheta(xpixels=1000, ypixels=1000, distobs=220.0, xoffset=0.0, yoffset=0.0, \
                 bangle=0.0, pangle=0.0, gamma=0.0):

  p = pangle*np.pi/180
  b0 = bangle*np.pi/180
  g = gamma*np.pi/180
  imscale = 1.97784*1024/xpixels
  scale = imscale/(180*60*60/np.pi)
  rsun=np.tan(np.arcsin(1/distobs))/scale
  x0 = xpixels/2-0.5 + xoffset
  y0 = ypixels/2-0.5 + yoffset
  xx = (np.linspace(0, xpixels-1, xpixels)-x0)/rsun
  zz = (np.linspace(0, ypixels-1, ypixels)-y0)/rsun
  xx, zz = np.meshgrid(xx, zz)

  xx /= np.cos(g)
  zz /= np.cos(b0)
  x1 =  xx*np.cos(p) + zz*np.sin(p)
  z1 = -xx*np.sin(p) + zz*np.cos(p)

  rr = np.sqrt(x1*x1 + z1*z1)
  index = (rr <= 1.0)
  r = np.ma.array(rr, mask=np.logical_not(index))
  angle = np.arctan2(x1,z1)
  theta = np.ma.array(np.abs(angle), mask=r.mask)
  phi = 0.0*theta
  index = (angle < 0.0)
  phi[index]=np.pi

  return (r, theta, phi)


def setplm(l, m, x, plm, dplm): 

# adapted from setplm.pro by Jesper Schou
# Set plm(*,l)=P_l^m (x) for l=m,lmax
# optionally sets dplm(*,l)={{dP_l^m} \over dx} (x) for l=m,lmax
# P_l^m 's are normalized to have \int_{-1}^1 (P_l^m (x))^2 dx = 1
# Works up to l \approx 1800, see ~/invnew/plm.f for details

  eps = 1.0e-12
  x1 = np.maximum(np.minimum(x,1-eps),eps-1)
#  x1 = x
  x2 = 1.0/(x1*x1-1.0)
  x3 = x1*x2
  c = np.sqrt((2*m+1)/2.0)
  for i in range(1,m+1):
    c *= -np.sqrt(1.0-0.5/i)
  plm[...,0] = c*(np.sqrt(1.0-x1*x1))**m
  if (l > m):
    c = np.sqrt(2.0*m+3.0)
    plm[...,1] = c*x1*plm[...,0]
  i = m+2
  while (i <= l):
    c1 = np.sqrt((4.0*i*i-1.0)/(i*i-m*m))
    c2 = np.sqrt(((2.0*i+1.0)*(i+m-1.0)*(i-m-1.0))/((2.0*i-3.0)*(i*i-m*m)))
    plm[...,i-m] = c1*x1*plm[...,i-m-1] - c2*plm[...,i-m-2]
    i+=1

  dplm[...,0] = m*x3*plm[...,0]
  i = m+1
  while (i <= l):
    c = np.sqrt((2.0*i+1.0)*(i*i-m*m)/(2.0*i-1))
    dplm[...,i-m] = i*x3*plm[...,i-m] - c*x2*plm[...,i-m-1]
    i+=1

  return (plm[...,l-m],dplm[...,l-m])


lsave=[0]
msave=[0]
nsave=[1]

def querylmn(index):

  global lsave, msave, nsave
  if (index < 0):
    if (index < -len(lsave)):
      i=-len(lsave)
    else:
      i=index
  else:
    i=-1
  lstr=catchc("Enter spherical harmonic degree (l): ",lsave[i])
  if lstr == 'q':
    return (-1,-1,-1)
  else:
    while True:
      try:
        lval=int(lstr)
        if lval < 0:
          print("Degree must be greater than or equal to zero. Try again.")
          lstr=catchc("Enter spherical harmonic degree (l): ",lsave[i])
          if lstr == 'q':
            return (-1,-1,-1)
          continue
        break
      except:
        print("Invalid number, try again.")
        lstr=catchc("Enter spherical harmonic degree (l): ", lsave[i])
        if lstr == 'q':
          return (-1,-1,-1)

  mstr=catchc("Enter azimuthal order (m): ",msave[i])
  if mstr == 'q':
    return (-1,-1,-1)
  else:
    while True:
      try:
        mval=int(mstr)
        if np.abs(mval) > lval:
          print("Azimuthal order must be in the range -l <= m <= l. Try again.")
          mstr=catchc("Enter azimuthal order (m): ",msave[i])
          if mstr == 'q':
            return (-1,-1,-1)
          continue
        break
      except:
        print("Invalid number, try again.")
        mstr=catchc("Enter azimuthal order (m): ",msave[i])
        if mstr == 'q':
          return (-1,-1,-1)

  nstr=catchc("Enter radial order (n): ",nsave[i])
  if nstr == 'q':
    return (-1,-1,-1)
  else:
    while True:
      try:
        nval=int(nstr)
        if nval < 0:
          print("Radial order must be greater than or equal to zero. Try again.")
          nstr=catchc("Enter radial order (n): ",nsave[i])
          if nstr == 'q':
            return (-1,-1,-1)
          continue
        break
      except:
        print("Invalid number, try again.")
        nstr=catchc("Enter radial order (n): ",nsave[i])
        if nstr == 'q':
          return (-1,-1,-1)

  lsave.append(lval)
  msave.append(mval)
  nsave.append(nval)
  return (lval,mval,nval)


varlist=['Vr', 'Vt', 'Vp', 'Vh', 'Vmag', 'Vsq']
colormap="seismic"
plotvar="Vr"
isave=0
capflag=0

def queryplotparms():

  global colormap, plotvar, capflag, isave
  print("You may enter 'l' to list options.")
  cmap=input("Enter name of colormap: ")
  while True:
    if (cmap == 'q'):
      return cmap
    if (cmap == ''):
      print("Using saved value colormap = %s." % colormap)
      cmap=colormap
      break
    if (cmap == 'l'):
      printcmaps()
      print("Current colormap is %s." % colormap)
      cmap=input("Enter name of colormap: ")
      continue
    elif (cmap not in plt.colormaps()):
      print("That colormap not registered, try again.")
      cmap=input("Enter name of colormap: ")
      continue
    else:
      break
  colormap=cmap
    
  pvar=input("Enter variable to plot: ")
  while True:
    if (pvar == 'q'):
      return pvar
    if (pvar == ''):
      print("Using saved value plotvar = %s." % plotvar)
      pvar=plotvar
      break
    if (pvar == 'l'):
      print("Available plotting variables are: ", end="")
      for i in varlist[:len(varlist)-1]:
        print(i, end=", ")
      print("and %s" % varlist[len(varlist)-1], end=".\n")
      print("Currently plotting %s." % plotvar)
      pvar=input("Enter variable to plot: ")
      continue
    elif (pvar not in varlist):
      print("That variable not available, try again.")
      pvar=input("Enter variable to plot: ")
      continue
    else:
      break
  plotvar=pvar

  cfstr=input("Enter caption flag: ")
  while True:
    if (cfstr == 'q'):
      return cfstr
    if (cfstr == ''):
      print("Using saved value caption flag = %i." % capflag)
      cflag=capflag
      break
    if (cfstr == 'l'):
      print("Available caption flags are: ")
      print("0 - no caption")
      print("1 - stacked caption")
      print("2 - single line caption")
      print("3 - one line per mode")
      print("Current caption flag is %i." % capflag)
      cfstr=input("Enter caption flag: ")
      continue
    elif (cfstr not in ['0','1','2','3']):
      print("Invalid caption flag, try again.")
      cfstr=input("Enter caption flag: ")
      continue
    else:
      cflag=int(cfstr)
      break
  capflag=cflag

  ss=input("Save output to file? (y/n) ")
  if (ss == 'q'):
    return ss
  if (ss != 'y'):
    isave=0
    print("Save off.")
  else:
    isave=1
    print("Save on.")

  print()


def catchc(prompt, saveval):

  valstr=input(prompt)
  if (valstr == 'q'):
    return valstr
  while (valstr == 'c'):
    c=queryplotparms()
    if (c == 'q'):
      return c
    valstr=input(prompt)
  if (valstr == ''):
    print("Using saved value %s." % str(saveval))
    valstr=str(saveval)
  return valstr

def printcmaps():

  print("Options include the following nonexhaustive list.  You may append '_r' to any name to reverse the map. \
More information is available at https://matplotlib.org/tutorials/colors/colormaps.html")

  cmaps = [('Perceptually Uniform Sequential', [
            'viridis', 'plasma', 'inferno', 'magma', 'cividis']),
         ('Sequential', [
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']),
         ('Sequential (2)', [
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink',
            'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia',
            'hot', 'afmhot', 'gist_heat', 'copper']),
         ('Diverging', [
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu',
            'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic']),
#         ('Cyclic', ['twilight', 'twilight_shifted', 'hsv']),
         ('Qualitative', [
            'Pastel1', 'Pastel2', 'Paired', 'Accent',
            'Dark2', 'Set1', 'Set2', 'Set3',
            'tab10', 'tab20', 'tab20b', 'tab20c']),
         ('Miscellaneous', ['hsv', 
            'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar'])]

  for c in cmaps:
    print("%s:" % c[0])
    for m in c[1]:
      print(m,end=" ")
    print("")


def nothing():
  return None

nframes=64
dpi=300
icolshift=0
pngdirfmt = './png_out/{:s}{:s}.l{:d}m{:d}n{:d}'
animate=nothing
calcimage=nothing
bgcol='white'
txtcol='black'
capblank=r'$\mathbf{\ell} = 5$' '\n' r'${m = 3}$' '\n' r'${n = 1}$'

def drawfigure(var, caption=None, fsize=5):

#  below was written for single images. same functionality is now 
#  achieved by setting vmin and vmax.  probably exactly the same 
#  but new code is clearer. 
#  for animations.
#  if (icolshift != 0) and plotvar in ['Vr','Vt','Vp']:
#    mn,mx=var.min(),var.max()
#    absarr=np.abs([mn,mx])
#    frac=(mx-mn)/(2*absarr.max())
#    table=cm.get_cmap(colormap,256/frac)
#    if (absarr[0]>absarr[1]):
#      newcolors=table(np.linspace(0,frac,256))
#    else:
#      newcolors=table(np.linspace(1-frac,1,256))
#    newmap=ListedColormap(newcolors, name='soshtmp')
#    cm.register_cmap(cmap=newmap)
#    usemap='soshtmp'
#  else:
#    usemap=colormap

  fig, ax = plt.subplots(num=1, figsize=(fsize, fsize), facecolor=bgcol)
  ax.clear()
#  im=ax.imshow(var,cmap=usemap)
  im=ax.imshow(var, origin='lower', cmap=colormap)

  if (icolshift == 1) and plotvar in ['Vr','Vt','Vp']:
    mn,mx=var.min(),var.max()
    maxabs=np.abs([mn,mx]).max()
#    im=ax.imshow(var, cmap=colormap, vmin=-maxabs, vmax=maxabs)
    im.set_clim(vmin=-maxabs, vmax=maxabs)
    print("Scaling to maxval = %f"%maxabs)
  elif (icolshift == 2):
    maxval=0.0
    for i in range(nframes):
      d=calcimage(i)
      val=np.abs(d).max()
      if (val > maxval):
        maxval=val
    if plotvar in ['Vr','Vt','Vp']:
#      im=ax.imshow(var, cmap=colormap, vmin=-maxval, vmax=maxval)
      im.set_clim(vmin=-maxval, vmax=maxval)
    else:
      im.set_clim(vmin=0.0, vmax=maxval)
    print("Scaling to maxval = %f"%maxval)
#  else:
#    im=ax.imshow(var, cmap=colormap)

  ax.set_axis_off()
  fig.canvas.draw()
  if (caption != None):
    ax.text(0.5,-0.07,caption,size='large',multialignment='center',horizontalalignment='center',verticalalignment='center',transform=ax.transAxes,color=txtcol)

  return (fig,im)

def savefigure(ianimate=1, label=''):

  l=lsave[-1]
  m=msave[-1]
  n=nsave[-1]
  if (ianimate == 0):
    savestr='png_out/'+plotvar+label+'.l%im%in%i.png' % (l,m,n)
    plt.savefig(savestr, dpi=dpi, facecolor=bgcol)
    print("File saved.")
  else:
    pngdir = pngdirfmt.format(plotvar, label, l, m, n)
    if not os.path.exists(pngdir):
      os.makedirs(pngdir)
    print("Writing files...", end="")
    for i in range(nframes):
      print(i, end=" ", flush=True)
      animate(i)
      fpath = os.path.join(pngdir, '{:03d}.png'.format(i))
      plt.savefig(fpath, dpi=dpi, facecolor=bgcol)
    print("done.")


