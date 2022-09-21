import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from scipy.optimize import curve_fit
import scipy
from numpy import exp

# parameters
d_list=[4.0,3.0,2.0,1.0,0.8,0.7,0.6,0.5,0.4,0.3,0.2]
n=10000
psips=1000
delta_t=0.005

# reading files
data_set=[]
data_per_d=[]
data_full={}
for d in d_list:
    for i in range(100):
        # INSERT YOUR FILEPATH HERE
        fn = "/media/veracrypt2/pdoc/uni/PA/output_protonpot/protonpot_averages_delta_t=" +str(delta_t)+"_n="+str(n)+"_#psips="+str(psips)+"_d"+str(d)+"_it="+str(i)+".txt"
        f = open(fn, "r")
        for line in f:
            readindex, readvalue = line.split()
            data_set.append(float(readvalue))
        f.close()
        data_per_d.append(data_set)
        data_set=[]
    # where all the data is kept; data_full has one list per entry in d_list. Each of these lists has 100 lists (these are the statistical attempts). And those lists are what is contained in the files that were just read (n or steplength entries of energies at each step)
    data_full[d]=data_per_d
    data_per_d=[]

# resorting the data so that for each d, there are n lists with each 100 entries fo the statistical attempts
new_dict= {}
inter_dict={}
resorted_data={}
for d in d_list:
    data_per_d=data_full[d]
    for i in range(n-1):
        for iteration in range(len(data_per_d)):
            new_dict[iteration]=data_per_d[iteration][i]
        inter_dict[i]=new_dict
        new_dict={}
    resorted_data[d]=inter_dict
    inter_dict={}

# calculating the average&standard deviation of the 100 attempts at each step n
average_n={}
average_d={}
standard_n={}
standard_d={}
sum_iter=0
for d in d_list:
    data_per_d=resorted_data[d]
    for i in range(n-1):
        data_set=data_per_d[i]
        for iteration in range(len(data_set)):
            sum_iter=sum_iter+data_set[iteration]
        average_n[i]=sum_iter/len(data_set)
        sum_sum_std=0
        for iteration in range(len(data_set)):
            sum_std=(average_n[i]-data_set[iteration])**2
            sum_sum_std+=sum_std
        standard_n[i]=math.sqrt((1/(n-1-1))*sum_sum_std)
        sum_iter=0
    average_d[d]=average_n
    standard_d[d]=standard_n
    average_n={}

# after a certain number of steps (chosen approximately by checking the graphs), the average energies don't change anymore. We calculate their average, gaining one final average energy per proton distance d
ave_d={}
stdv_d={}
sumofstandin=0
for d in d_list:
    average_n=average_d[d]
    for index in range (4000, n-1):
        if d not in ave_d:
            ave_d[d]=0   
        ave_d[d]=ave_d[d]+average_n[index]
    ave_d[d]=ave_d[d]/(n-1-4000)
    sumofstandin=0
    for index in range (4000, n-1):
        standin=(ave_d[d]-average_n[index])**2
        sumofstandin=sumofstandin+standin
    stdv_d[d]=math.sqrt((1/(n-1-4000-1))*sumofstandin)
    
# statistical error
i=0
sumofstandin=0
stat_error={}
for d in d_list:
    standin=0
    standard_n=standard_d[d]
    average_n=average_d[d]
    for index in range (4000, n-1):
        standin+=average_n[index]
    standin=((standin-average_n[i])*standard_n[i])**2
    sumofstandin+=standin
    stat_error[d]=(1/(n-1-4000))*math.sqrt(sumofstandin)
    
# full error
full_error={}
for d in d_list:
    full_error[d]=stdv_d[d]+stat_error[d]

# plotting needs lists, previously we had dicts
plot_ave=[]
for d in ave_d:
    plot_ave.append(ave_d[d])
plot_std=[]
for d in full_error:
    plot_std.append(full_error[d])

# fitting to the morse potential
# defining the morse potential (R, as the variable, has to be the first entry for curve_fit to work)
def morse(R,D_e,a,R_e,v):
    return -D_e*(1-np.exp(-a*(R-R_e)))**2+v
# starter value "guesses" (visual)    
p0=[-0.9,1,0.6,-1]
# original d_list was the distances of protons from zero; for sensible fitting, we need the full interproton distance
d_list_2=[]
for j in d_list:
    d_list_2.append(j*2)
# using scipy's curve_fit to compute parameters
params, pcov = curve_fit(morse, d_list_2, plot_ave, p0, maxfev=50000)
D_e, a, R_e, v = params[0], params[1], params[2], params[3]
lowestmorse=morse(R_e,D_e,a,R_e,v)
# computing the base frequency; 
omega_0=(a/(2*np.pi))*np.sqrt(-2*D_e/(1836))
# calculating the error of the fitted parameters
perr = np.sqrt(np.diag(pcov))
morseerr = lowestmorse*perr[0]/D_e+lowestmorse*2*a*perr[1]
omegaerr = (omega_0/a)*perr[1]+omega_0/(2*np.sqrt(2)*D_e)*perr[0]

# printing results!
print("D_e=",D_e,"+/-",perr[0],"a=",a,"+/-",perr[1],"R_e=",R_e,"+/-",perr[2])
print("best interproton distance=", R_e,"+/-",perr[2])
print("ground state energy=",lowestmorse,"+/-",morseerr)
print("base frequency=",omega_0,"+/-",omegaerr)

# creating SI values
R_e_SI=R_e*scipy.constants.physical_constants["Bohr radius"][0]
R_e_err_SI=perr[2]*scipy.constants.physical_constants["Bohr radius"][0]
lowestmorse_SI=lowestmorse*scipy.constants.physical_constants["Hartree energy"][0]
morseerr_SI=morseerr*lowestmorse*scipy.constants.physical_constants["Hartree energy"][0]
omega_0_SI=omega_0*(scipy.constants.physical_constants["Hartree energy"][0])/scipy.constants.hbar
omegaerr_SI=omegaerr*(scipy.constants.physical_constants["Hartree energy"][0])/scipy.constants.hbar

# printing SI results
print("In SI:")
print("best interproton distance=", "{:e}".format(R_e_SI),"+/-", "{:e}".format(R_e_err_SI))
print("ground state energy=","{:e}".format(lowestmorse_SI),"+/-","{:e}".format(morseerr_SI))
print("base frequency=","{:e}".format(omega_0_SI),"+/-","{:e}".format(omegaerr_SI))

print("base frequency in cm-1:", omega_0_SI/scipy.constants.c/100, "+/-", omegaerr_SI/scipy.constants.c/100)

# for plotting (no arithmetic operations with lists)
d_np=np.array(d_list_2)
morsefit = morse(d_np,D_e,a,R_e,v)

# plotting the fit curve!
plt.plot(d_list_2, morsefit,label="morse fit")
plt.errorbar(d_list_2, plot_ave, yerr=plot_std,fmt=".",ecolor="red",barsabove=True,label="averages")
plt.xlabel('Interproton distances [Bohr]')
plt.ylabel('Potential energy [Hartree]')
plt.legend(loc='best')
plt.grid(True)
plt.show() 

# plotting the averages per d from before!
colors = plt.cm.rainbow(np.linspace(0, 1, len(d_list)))
plot_av=[]
plot_st=[]
x_list=[x for x in range(n-1)]
j=0
for d in d_list:
    standin_d=average_d[d]
    std_standin_d=standard_d[d]
    for i in range(n-1):
        plot_av.append(standin_d[i])
        plot_st.append(std_standin_d[i])
    plt.plot(x_list, plot_av, color=colors[j], label=d)
    #plt.errorbar(x_list, plot_av, yerr=plot_st, color=colors[j], alpha=0.025)
    plot_av=[]
    plot_st=[]
    j+=1
plt.xlabel("steps")
plt.ylabel("Potential energy [Hartree]")
plt.xlim(0,n-1)    
plt.legend(loc="best")
plt.show()
