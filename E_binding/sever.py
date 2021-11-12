#!python
import sys

frag = int(sys.argv[1])

a=open("/home/xg23/CPEB3/SUMO/thread/cpeb3.seq", mode='r')
b=''

for i in a:
   b += i

c=b.strip()

for j in range(len(c)-frag+1):
   out=open('sim-'+str(j)+'.fasta', 'w')
   out.write('>SIM')
   out.write('\n')
   out.write(b[j:j+frag])
   out.write('\n')
   out.close()

