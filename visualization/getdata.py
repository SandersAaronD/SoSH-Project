import soshdata
#import os
#import sys
import socket

def pdreceive(port):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.bind(('localhost', port))
    s.listen(1)
 
    while True:
        conn, addr = s.accept()
        for line in conn.makefile():
            yield line
        conn.close()
    s.close()

def pdsend(port, message):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.connect(('localhost', port))
    sent=s.send(message)
    s.close()
    return sent

#def sendtopd(message):
#  os.system("echo '" + message + ";' | pdsend 9119")

#while (True):
#  try:
#    recstr = sys.stdin.readline()
#  except:
#    break
#  recstr=recstr.strip(';\n')

print("Listening on port 1991.  You may now start the SoSH Tool in Pure Data.")
print("When done you may stop this script with ctrl-C.")

for message in pdreceive(1991):

  recstr=message.strip(';\n')
  print('"'+recstr+'" was received.')
  words=recstr.split()
  instr=words[0]
  day=int(words[1])
  l=int(words[2])
  m=int(words[3])
  stem=words[4]
  splitflag=int(words[5])
  makemsplit=(splitflag != 0)

  wf=soshdata.getwavfiles(instrument=instr, day=day, l=l, m=m, delfits=False)
  if (stem == 'average.modes'):
    pf=soshdata.getmodeparms(instrument=instr, day='average', makemsplit=makemsplit)
  else:
    pf=soshdata.getmodeparms(instrument=instr, day=day, makemsplit=makemsplit)

  if (wf != None and pf != None):
    print("Requested files:",wf,pf) 
#    sendtopd("go")
    sent=pdsend(9119, "go;".encode())
  else:
#    sendtopd("stop")
    sent=pdsend(9119, "stop;".encode())
  print()

