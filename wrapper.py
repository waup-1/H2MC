from h2randomwalk import h2randomwalk
import matplotlib.pyplot as plt
import time,os,sys

tstart = time.time()

d_list=[10,9,8,7,6,5,4,3,2,1,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01]
n=10000
psips=1000
delta_t=0.005

def func(*args):
    print ("processing",args, "...")
    d=args[0]
    i=args[1]
    f = open("protonpot_averages_delta_t=" +str(delta_t)+"_n="+str(n)+"_#psips="+str(psips)+"_d"+str(d)+"_it="+str(i)+".txt", "w")
    
    average_list=h2randomwalk(n,psips,d,delta_t,i*10000+int(d*100),i)
    print("repetition done: {:4.1f} {:2d}".format(d, i))

    # plotting the average energies
    ##############################
    fig = plt.figure()
    plt.plot(average_list)
    plt.title("delta_t=" +str(delta_t)+", n="+str(n)+ ", #psips="+str(psips)+", and d:"+str(d))
    plt.xlabel("steps")
    plt.ylabel("Energy [hartree]")
    #plt.xlim(0, n-1)
    fig.savefig("protonpot_averages_delta_t=" +str(delta_t)+"_n="+str(n)+"_#psips="+str(psips)+"_d"+str(d)+"_it="+str(i)+".png")
    plt.close(fig)
    j=0
    for average in average_list:
        print(j, average, file=f)
        j=j+1
    return "for {:4.1f} {:2d}  #averages: {:d}".format(d,i, j)


d = float (sys.argv[1])
for i in range(100):
  print (func (d, i))
print ("Elapsed: {:6.1f} seconds".format(time.time()-tstart))
