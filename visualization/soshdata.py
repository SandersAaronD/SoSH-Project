import numpy as np
import os
import re
import drms
from astropy.io import fits
from urllib.request import urlretrieve
import soundfile as sf
import matplotlib.pyplot as plt

import warnings 
warnings.simplefilter(action='ignore', category=FutureWarning)

#datadir = '/home/tplarson/solar/data'
datadir = '../data'

def getmodeparms(instrument=None, day=None, lmin=0, lmax=300, makemsplit=True):

  if not os.path.exists(datadir):
    os.makedirs(datadir)

  if (instrument == None):
    inst=input("Enter name of instrument: ")
  else:
    inst=instrument

  if (day == None):
    daystr=input("Enter day number: ")
  else:
    daystr=str(day)

  mfile=getmfile(inst.lower(), daystr)
  if (mfile == None):
    return None

  if (makemsplit):
    msfile=mfile+'.msplit'

    if os.path.exists(msfile):
      print("Msplit file already exists.")
    else:
      msfile=mkmsplit(mfile,lmin=lmin,lmax=lmax)

    return mfile,msfile
  else:
    return mfile


def getmfile(inst, daystr):

  if (inst == 'hmi'):
    qblank = 'hmi.V_sht_modes[{:d}d][0][300][138240]'
    avgurl = 'http://sun.stanford.edu/~tplarson/audio/HMI/hmi.average.modes'
  elif (inst == 'mdi'):
    qblank = 'mdi.vw_V_sht_modes[{:d}d][0][300][103680]'
    avgurl = 'http://sun.stanford.edu/~tplarson/audio/MDI/mdi.average.modes'
  else:
    print("Invalid instrument, try again.")
    return None

  if (daystr == 'average'):
    mblank = '%s.average.modes'
    mfile=os.path.join(datadir,mblank % inst)
    if os.path.exists(mfile):
      print("Mode fits file already exists.")
    else:
      try:
        urlretrieve(avgurl, mfile)
      except:
        print("Failure of urlretrieve() for %s." % avgurl)
        return None
  else:
    mblank = '%s.%sd.modes'
    mfile=os.path.join(datadir,mblank % (inst, daystr))
    if os.path.exists(mfile):
      print("Mode fits file already exists.")
    else:
      try:
        day=int(daystr)
      except:
        print("Invalid number, try again.")
        return None
      qstr=qblank.format(day)
      c = drms.Client()
      m = c.query(qstr, seg='m6')
      if (len(m) == 0):
        print("No data found for that day number, try again.")
        return None
      url = 'http://jsoc.stanford.edu' + m.m6[0]
      try:
        urlretrieve(url, mfile)
      except:
        print("Failure of urlretrieve() for %s." % url)
        return None

  return mfile


def mkmsplit(mfile, lmin=0, lmax=300):

  modeparms=np.loadtxt(mfile)
  lmod=modeparms[:,0]
  nmod=modeparms[:,1]

  outfmt='{:d} {:d} {:d} {:f} {:f}'
  outfile = mfile + '.msplit'
  file=open(outfile,'w')

  l=lmin
  while (l <= lmax):
    pols=apols(l,6)
    ind = (lmod == l)
    nlist = nmod[ind]
    for n in nlist:
      ind2 = (lmod == l) & (nmod == n)
      ai=np.append([0.0],modeparms[ind2,12:18]/1000)
      ei=np.append([0.0],modeparms[ind2,18:24]/1000)
      fx=np.matmul(pols,ai)
      freq=modeparms[ind2,2]+fx
      var=np.matmul((pols*pols),(ei*ei))
      sig=modeparms[ind2,7]+np.sqrt(var)
      for m in range(-l,l+1):
        outstr=outfmt.format(int(l),int(n),int(m),freq[m+l],sig[m+l])
        file.write('%s\n' % (outstr))
    l+=1
  file.close()
  return outfile


def getwavfiles(instrument=None, day=None, l=None, m=None, delfits=False):

  if not os.path.exists(datadir):
    os.makedirs(datadir)

  if (instrument == None):
    inst=input("Enter name of instrument: ")
  else:
    inst=instrument

  if (inst.lower() == 'hmi'):
    hind=1
  elif (inst.lower() == 'mdi'):
    hind=0
  else:
    print("Invalid instrument, try again.")
    return None

  if (day == None):
    daystr=input("Enter day number: ")
    try:
      day=int(daystr)
    except:
      print("Invalid number, try again.")
      return None
  if (l == None):
    lstr=input("Enter spherical harmonic degree: ")
    try:
      l=int(lstr)
    except:
      print("Invalid number, try again.")
      return None
  if (m == None):
    mstr=input("Enter azimuthal order: ")
    try:
      m=int(mstr)
      m=np.abs(m)
      if (m > l):
        print("Abs(m) must be less than l, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None

  wblank = '%s.%id_72d.l=%i_m=%i_data%s.wav'
  wfile=os.path.join(datadir,wblank % (inst.lower(),day,l,m,'r'))
  wfile2=os.path.join(datadir,wblank % (inst.lower(),day,l,m,'i'))
  if os.path.exists(wfile) and (m==0 or os.path.exists(wfile2)):
    if (m == 0):
      print("Wav file already exists.")
      return wfile
    else:
      print("Wav files already exists.")
      return wfile,wfile2

  ffile=getfitsfile(inst.lower(),day,l)
  if (ffile == None):
    return None

  wf=convertfitstowav(ffile, hind, m, wfile)
  if (wf == None):
    return None

  if delfits:
    os.remove(ffile)
  return wf


def getfitsfile(inst, day, l):

  if (inst == 'hmi'):
    qblank = 'hmi.V_sht_gf_72d[{:d}d][{:d}][{:d}][138240]'
  elif (inst == 'mdi'):
    qblank = 'mdi.vw_V_sht_gf_72d[{:d}d][{:d}][{:d}][103680]'
  else:
    print("Invalid instrument, try again.")
    return None

  fblank='%s.%id.%i.fits'
  ffile=os.path.join(datadir,fblank % (inst,day,l))
  if os.path.exists(ffile):
    print("Fits file already exists.")
  else:
    qstr=qblank.format(day,l,l)
    c = drms.Client()
    f = c.query(qstr,seg='data')
    if (len(f) == 0):
      print("No data found for that day number and degree, try again.")
      return None
    url = 'http://jsoc.stanford.edu' + f.data[0]
    try:
      urlretrieve(url, ffile)
    except:
      print("Failure of urlretrieve() for %s." % url)
      return None

  return ffile


def convertfitstowav(ffile, hind, m, wfile):

  hdul = fits.open(ffile)
  h=hdul[hind]
  d=h.data
  x=d[m,0::2]/np.abs(d[m,:]).max()
  sf.write(wfile, x, 8000, subtype='FLOAT')

  if (m > 0):
    wfile2=wfile.replace('datar','datai')
    x=d[m,1::2]/np.abs(d[m,:]).max()
    sf.write(wfile2, x, 8000, subtype='FLOAT')

  hdul.close()
  if (m==0):
    return wfile
  else:
    return wfile,wfile2


def apols(lin, amaxin):

# adapted from IDL procedure apols.pro written by Jesper Schou

  l=np.long(lin)
  amax=np.minimum(amaxin,2*l)

  pols=np.zeros((2*l+1,amaxin+1))
# It is ESSENTIAL that x is set up such that x(-m)=-x(m) to full machine
# accuracy or that the re-orthogonalization is done with respect to all
# previous polynomials (second option below).
# x=(dindgen(2*l+1)-l)/l
  x=np.linspace(-1,1,2*l+1)

  pols[:,0]=1/np.sqrt(2*l+1.0)
  if (amax >= 1):
    pols[:,1]=x/np.sqrt((l+1.0)*(2*l+1.0)/(3.0*l))
# for n=2l,amax do begin
# Set up polynomials using exact recurrence relation.
  for n in range(2,amax+1):
    a0=2.0*l*(2*n-1.0)/(n*(2*l-n+1))
    b0=-(n-1.0)/n*(2*l+n)/(2*l-n+1)
    a=a0*np.sqrt((2*n+1.0)/(2*n-1.0)*(2*l-n+1.0)/(2*l+n+1.0))
    b=b0*np.sqrt((2*n+1.0)/(2*n-3.0)*(2*l-n+1.0)*(2*l-n+2.0)/(2*l+n+1.0)/(2*l+n))
    help=a*x*pols[:,n-1]+b*pols[:,n-2]
# Turns out that roundoff kills the algorithm, so we have to re-orthogonalize.
# Two choices here. First one is twice as fast and more accurate, but
# depends on the rounding being done in the standard IEEE way.
#  for j=n-2,0,-2 do begin 
    for j in range(n-2,-1,-2):    
# This choice is robust to the roundoff, but twice as slow and generally
# somewhat less accurate
# for j=n-1,0,-1 do begin
      help=help-pols[:,j]*np.sum(help*pols[:,j])
#  end
# Reset norm to 1.
    pols[:,n]=help/np.sqrt(np.sum(np.square(help)))

# Reset polynomials to have P_l(l)=l by setting overall norm.
# Note that this results in more accurate overall normalization, but
# that P_l(l) may be very far from l. This is the key to the algorithm.
  c=l**2*(2*l+1.0)
# for n=0l,amax do begin
  for n in range(0,amax+1):
    c=c*(2*l+n+1.0)/(2*l-n+1.0)
    pols[:,n]=pols[:,n]*np.sqrt(c/(2*n+1))

  return pols


def splitwavfile(wavfile=None, splitfactor=2, splitboth=True):

  wblank = '%s.%id_%id.l=%i_m=%i_data%s.wav'
  outlist=[]

  if (wavfile == None):
    inst=input("Enter name of instrument: ")
    inst=inst.lower()
    if (inst != 'mdi' and inst != 'hmi'):
      print("Invalid instrument, try again.")
      return None
    daystr=input("Enter day number: ")
    try:
      firstday=int(daystr)
    except:
      print("Invalid number, try again.")
      return None
    ndstr=input("Enter number of days: ")
    try:
      ndaysin=int(ndstr)
    except:
      print("Invalid number, try again.")
      return None
    lstr=input("Enter spherical harmonic degree: ")
    try:
      l=int(lstr)
      if (l < 0):
        print("Degree l must be nonnegative, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None
    mstr=input("Enter azimuthal order: ")
    try:
      m=int(mstr)
      m=np.abs(m)
      if (m > l):
        print("Abs(m) must be less than l, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None
    ri='r'
    wavfile=os.path.join(datadir,wblank % (inst,firstday,ndaysin,l,m,ri))
  else:
    regmatch = re.search(r"(\S+)\.(\d+)d_(\d+)d\.l=(\d+)_m=(\d+)_data(\S+).wav", wavfile)
    if (regmatch == None):
      print("Invalid file name, try again.")
      return None
    else:
      inst=regmatch.group(1)
      firstday=int(regmatch.group(2))
      ndaysin=int(regmatch.group(3))
      l=int(regmatch.group(4))
      m=int(regmatch.group(5))
      ri=regmatch.group(6)

  test= ndaysin % splitfactor
  if (test != 0):
    print("Timeseries not divisible by that factor.")
    return None

#  wfile = os.path.join(datadir,wavfile)
  try:
    x, sr = sf.read(wavfile)
  except:
    print("Invalid file, try again.")
    return None

  if (inst == 'mdi'):
    nperday = 1440
  elif (inst == 'hmi'):
    nperday = 1920

  ndaysout = ndaysin / splitfactor
  ndt=int(ndaysout*nperday)
  day = firstday
  index=0
  while (day < firstday+ndaysin):
    xout=x[index:index+ndt]
    wfileout = os.path.join(datadir,wblank % (inst,day,ndaysout,l,m,ri))
    sf.write(wfileout, xout, sr, subtype='FLOAT')
    day += ndaysout
    index += ndt
    outlist.append(wfileout)

  if (splitboth and m > 0):
    if (ri == 'i'):
      ir = 'r'
    elif (ri == 'r'):
      ir = 'i'
    wfile=wblank % (inst,firstday,ndaysin,l,m,ir)
    wavfile = os.path.join(datadir,wfile)
    try:
      x, sr = sf.read(wavfile)
    except:
      print("Invalid file, try again.")
      return None
    day=firstday
    index=0
    while (day < firstday+ndaysin):
      xout=x[index:index+ndt]
      wfileout = os.path.join(datadir,wblank % (inst,day,ndaysout,l,m,ir))
      sf.write(wfileout, xout, sr, subtype='FLOAT')
      day += ndaysout
      index += ndt
      outlist.append(wfileout)

  return outlist


def scanspectrum(wavfile=None,downshift=4,firstbin=None,lastbin=None,nbins=200,ndt=10000,binstep=None,ramp=50):

  if (binstep == None):
    binstep = int(nbins/2)

  wblank = '%s.%id_%id.l=%i_m=%i_data%s.wav'
  if (wavfile == None):
    inst=input("Enter name of instrument: ")
    inst=inst.lower()
    if (inst != 'mdi' and inst != 'hmi'):
      print("Invalid instrument, try again.")
      return None
    daystr=input("Enter day number: ")
    try:
      firstday=int(daystr)
    except:
      print("Invalid number, try again.")
      return None
    ndstr=input("Enter number of days: ")
    try:
      ndaysin=int(ndstr)
    except:
      print("Invalid number, try again.")
      return None
    lstr=input("Enter spherical harmonic degree: ")
    try:
      l=int(lstr)
      if (l < 0):
        print("Degree l must be nonnegative, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None
    mstr=input("Enter azimuthal order: ")
    try:
      m=int(mstr)
      signedm=m
      m=np.abs(m)
      if (m > l):
        print("Abs(m) must be less than l, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None
    ri='r'
    wavfile=os.path.join(datadir,wblank % (inst,firstday,ndaysin,l,m,ri))
    print("Wav file is %s" % wavfile)
  else:
    regmatch = re.search(r"(\S+)\.(\d+)d_(\d+)d\.l=(\d+)_m=(\d+)_data(\S+).wav", wavfile)
    if (regmatch == None):
      print("Invalid file name, try again.")
      return None
    else:
      inst=regmatch.group(1)
      firstday=int(regmatch.group(2))
      ndaysin=int(regmatch.group(3))
      l=int(regmatch.group(4))
      m=int(regmatch.group(5))
      signedm=m
      ri=regmatch.group(6)

  try:
    x, sr = sf.read(wavfile)
  except:
    print("Invalid file, try again.")
    return None

  z=np.zeros((len(x)),dtype=np.complex_)
  if (ri == 'r'):
    z.real=x
  else:
    z.imag=x

  if (m > 0):
    if (ri == 'i'):
      ir = 'r'
    elif (ri == 'r'):
      ir = 'i'
    wfile=wblank % (inst,firstday,ndaysin,l,m,ir)
    wavfile = os.path.join(datadir,wfile)
    try:
      x, sr = sf.read(wavfile)
    except:
      print("Invalid file, try again.")
      return None
  else:
    x=0.0*x

  if (ri == 'r'):
    z.imag=x
  else:
    z.real=x

  window=np.ones(ndt)
  if (ramp != 0):
    rampsamples=int(sr*ramp/1000)
    if (rampsamples > ndt/2):
      print("Ramp longer than ndt/2, try again.")
      return None
    window[0:rampsamples]=np.linspace(0.0,1.0,rampsamples,endpoint=False)
    window[ndt-1:ndt-rampsamples-1:-1]=np.linspace(0.0,1.0,rampsamples,endpoint=False)

  nx=len(x)
  nx2=int(nx/2)
  if (firstbin == None):
    firstbin = int(nx2/6)
  if (lastbin == None):
    lastbin = int(0.9*nx2)
  ftz = np.fft.fft(z)
  freq = np.fft.fftfreq(nx)
  if (signedm >= 0):
    mag = np.abs(ftz[0:nx2])
    phase = np.arctan2(ftz[0:nx2].imag,ftz[0:nx2].real)
    f = freq[0:nx2]/downshift
  else:
    mag = np.abs(ftz[nx:nx2:-1])
    phase = np.arctan2(-1*ftz[nx:nx2:-1].imag,ftz[nx:nx2:-1].real)
    f = -1*freq[nx:nx2:-1]/downshift

  nstep=int((lastbin-firstbin)/binstep)
  tlen=nstep*ndt
  t=np.linspace(0,ndt-1,ndt)
  xout=np.zeros((tlen))
  bindex=firstbin
  tindex=0
  t0in=0
  while (bindex + nbins <= lastbin):
    for i in range(nbins):
      xout[tindex:tindex+ndt]+=mag[bindex+i]*np.sin(2*np.pi*f[bindex+i]*(t+tindex)+phase[bindex+i])
    xout[tindex:tindex+ndt]*=window
    bindex+=binstep
    tindex+=ndt

  mx=np.abs(xout).max()
  xout /= mx
  sf.write('output.wav', xout, sr, subtype='FLOAT')

  return xout

plotnu=np.array([])

def plotspectrum(wavfile=None, transpose=None):

  wblank = '%s.%id_%id.l=%i_m=%i_data%s.wav'
  if (wavfile == None):
    inst=input("Enter name of instrument: ")
    inst=inst.lower()
    if (inst != 'mdi' and inst != 'hmi'):
      print("Invalid instrument, try again.")
      return None
    daystr=input("Enter day number: ")
    try:
      firstday=int(daystr)
    except:
      print("Invalid number, try again.")
      return None
    ndstr=input("Enter number of days: ")
    try:
      ndaysin=int(ndstr)
    except:
      print("Invalid number, try again.")
      return None
    lstr=input("Enter spherical harmonic degree: ")
    try:
      l=int(lstr)
      if (l < 0):
        print("Degree l must be nonnegative, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None
    mstr=input("Enter azimuthal order: ")
    try:
      m=int(mstr)
      signedm=m
      m=np.abs(m)
      if (m > l):
        print("Abs(m) must be less than l, try again.")
        return None
    except:
      print("Invalid number, try again.")
      return None
    ri='r'
    wavfile=os.path.join(datadir,wblank % (inst,firstday,ndaysin,l,m,ri))
    print("Wav file is %s" % wavfile)
  else:
    regmatch = re.search(r"(\S+)\.(\d+)d_(\d+)d\.l=(\d+)_m=(\d+)_data(\S+).wav", wavfile)
    if (regmatch == None):
      print("Invalid file name, try again.")
      return None
    else:
      inst=regmatch.group(1).split('/')[-1]
      firstday=int(regmatch.group(2))
      ndaysin=int(regmatch.group(3))
      l=int(regmatch.group(4))
      m=int(regmatch.group(5))
      signedm=m
      ri=regmatch.group(6)

  try:
    x, sr = sf.read(wavfile)
  except:
    print("Invalid file, try again.")
    return None

  z=np.zeros((len(x)),dtype=np.complex_)
  if (ri == 'r'):
    z.real=x
  else:
    z.imag=x

  if (m > 0):
    if (ri == 'i'):
      ir = 'r'
    elif (ri == 'r'):
      ir = 'i'
    wfile=wblank % (inst,firstday,ndaysin,l,m,ir)
    wavfile = os.path.join(datadir,wfile)
    try:
      x, sr = sf.read(wavfile)
    except:
      print("Invalid file, try again.")
      return None
  else:
    x=0.0*x

  if (ri == 'r'):
    z.imag=x
  else:
    z.real=x

  if (inst == 'mdi'):
    cadence=60.0
  elif (inst == 'hmi'):
    cadence=45.0

  if (transpose == None):
    fmult = 1000000/cadence
  else:
    fmult = transpose

  nx=len(x)
  nx2=int(nx/2)
  ftz = np.fft.fft(z)
  freq = np.fft.fftfreq(nx)
  if (signedm >= 0):
    pow = np.abs(ftz[0:nx2])**2
    f = freq[0:nx2] * fmult
  else:
    pow = np.abs(ftz[nx:nx2:-1])**2
    f = -1*freq[nx:nx2:-1] * fmult

  global plotnu
  plotnu = f

  return pow


def plotmodefit(l,n,m,mfile=None,transpose=None):

  if (mfile == None):
    inst=input("Enter name of instrument: ")
    inst=inst.lower()
    if (inst != 'mdi' and inst != 'hmi'):
      print("Invalid instrument, try again.")
      return None
    daystr=input("Enter day number: ")
    if (daystr == 'average'):
      mblank = '%s.average.modes'
      mfile=os.path.join(datadir,mblank % inst)
    else:
      try:
        day=int(daystr)
      except:
        print("Invalid number, try again.")
        return None
      mblank = '%s.%sd.modes'
      mfile=os.path.join(datadir,mblank % (inst, daystr))

  try:
    modeparms=np.loadtxt(mfile)
  except:
    print("File not found.")
    return None

  lmod=modeparms[:,0]
  nmod=modeparms[:,1]
  pols=apols(l,6)
  ind = (lmod == l) & (nmod == n)
  ai=np.append([0.0],modeparms[ind,12:18]/1000)
  fx=np.matmul(pols,ai)
  freq=modeparms[ind,2]+fx
  amp=modeparms[ind,3]
  wid=modeparms[ind,4]
  nu0=freq[m+l]
  if (transpose != None):
    nu0 *= transpose
    wid *= transpose
  p = 2*wid*amp*amp/(wid*wid + 4*(plotnu-nu0)*(plotnu-nu0))
  print("Peak frequency = %f." % nu0)

  return p


def mkfig(pow,p,interval=None,noise=50,ampmult=200):

  fig, ax = plt.subplots()
  ax.plot(plotnu,np.sqrt(pow))
  ax.plot(plotnu,ampmult*np.sqrt(p)+noise)
  if (interval != None):
    ax.set_xlim(interval)

  fig.show()

  return fig,ax


