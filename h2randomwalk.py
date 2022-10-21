def h2randomwalk(n, psips, d, delta_t, thread_id, iteration):
    #n = 10000 is the number of steps(increase in the value of n increases the complexity of graph)
    #psips = 1000 number of psips
    #d = 0.7 placement of protons; along x-axis, symmetrical around center
    #delta_t = 0.005 time steps should be between 0.05 and 0.001

    import numpy as np 
    import math
    import matplotlib.pyplot as plt
    #import random 
    np.random.seed(thread_id)

    # basics
    ########
    D = 0.5 #diffusion constant
    mu = math.sqrt(2*delta_t*D) #for gauss
    sigma = mu #for gauss


    #initializing functions
    #######################
    #potential function
    def pot (coor_psips,d):#coordinates is a numpy array of vectors for each psip, the vectors have 6 values (3D coordinates for both electrons); d is the distance between 0 and the protons
        i = 0
        potential = [0 for x in range(len(coor_psips))] #I do not know if this is necessary; the code worked with and without it
        for vector in coor_psips:
            r_1a = math.sqrt((vector[0]-d)**2+vector[1]**2+vector[2]**2) #distance electron1 protona
            r_1b = math.sqrt((vector[0]+d)**2+vector[1]**2+vector[2]**2) #distance electron1 protonb
            r_2a = math.sqrt((vector[3]-d)**2+vector[4]**2+vector[5]**2) #distance electron2 protona
            r_2b = math.sqrt((vector[3]+d)**2+vector[4]**2+vector[5]**2) #distance electron2 protonb
            r_12 = math.sqrt((vector[0]-vector[3])**2+(vector[1]-vector[4])**2+(vector[2]-vector[5])**2) #distance electron1 electron2
            r_ab = 2*d #distance protona and protonb
            potential[i] = -1/r_1a-1/r_1b-1/r_2a-1/r_2b+1/r_12+1/r_ab #list of potential energies; each i=each psip
            i = i+1
        return potential

    #calculates the average potential of all psips for end result and pot_ref
    def ave_pot (pot_list):
        j = 0
        summe = 0
        for i in range(len(pot_list)):
            summe = summe + pot_list[i]
            j = j+1
        return summe/j
    
    #this walks the particles 
    def random_walk (coor_psips): #coordinates is a numpy array of vectors for each psip, the vectors have 6 values (3D coordinates for both electrons)
        i = 0
        for vector in coor_psips: #cycling through all psips
            for j in range(len(vector)): #cycling through all coordinates per psip
                coor_psips[i,j] = coor_psips[i,j] + np.random.normal(0, sigma)
            i = i+1


    # initializing values
    #####################
    coor_psips = np.random.normal(0, 1, (psips, 6)) #coordinates, initialized in random places
    initial_coor = coor_psips.copy() #to mark where the psips started
    births=0 #counter
    birth_dict = {} #record
    deaths=0 #counter
    death_dict = {} #record
    average_list = [0 for x in range(n-1)] #list of average potentials per step
    memorydanger = False #in case the births get out of hand, so the computer doesn't freeze


    # actual program
    ################
    for j in range(0, n-1): #each loop does a step
        potential = [0 for x in range(len(coor_psips))] #to fill potentials; have to initialize every time because the number of psips=size changes
        random_walk(coor_psips) #walking all coordinates of all psips
        potential = pot(coor_psips,d) #calculate potentials for each psip
        average_list[j] = ave_pot (potential) #calculate average of potentials
        index=0 
        birth_index=[] #to add birthplaces to coor_psips
        death_index=[] #to remove deaths from coor_psips
        while 1==1: #check births and deaths for each psip
            vector = coor_psips[index] #doing this complicated bit with the index because it's useful later
            k = (potential[index] - average_list[j])
            p_birth = (-1)*k*delta_t #likelihood for a birth
            p_death = k*delta_t #likelihood for a death
            n_rand = np.random.rand() #rolling the dice
            if k < 0:
                if p_birth > n_rand:                 
                     birth_dict[index] = vector.copy() #for plotting births and adding birthplaces
                     birth_index.append(index) #to find those birthplaces
                     births = births + 1 #counter
                     #print("birth", len(coor_psips), j)
            if k > 0:
                if p_death > n_rand:                 
                     death_dict[index] = vector.copy() #for plotting the deaths
                     death_index.append(index) #to find deathplaces
                     deaths = deaths + 1 #counter
                     #print("death", len(coor_psips), j)
            index=index+1
            if index > len(coor_psips)-1: #to get out of the while loop after all coordinates were examined
                break
            if len(coor_psips) > 10000: #to save the computer
                print("memory danger!")
                memorydanger = True
                break
        if memorydanger == True:
            break
        # deleting the dead psips
        dummy_deaths=0
        for index in death_index:
            coor_psips = np.delete(coor_psips, index-dummy_deaths, axis=0) #delete dead psips
            dummy_deaths=dummy_deaths+1
        # birthing new psips
        for index in birth_index:
            coor_psips=np.append(coor_psips, [birth_dict[index]], axis=0) #adding another psip at this location, but at the end of the array

    # sanity check                 
    print("d:",d, "it:", iteration, " births:", births, "deaths:", deaths, "Last energy:", average_list[j])

    # plotting prep
    ###############
    birth_array = np.empty((births, 6))
    i=0
    for index in birth_dict:
        birth_array[i]=birth_dict[index]
        i=i+1
    
    death_array = np.empty((deaths, 6))
    i=0
    for index in death_dict:
        death_array[i]=death_dict[index]
        i=i+1


    # plotting the r1 x r2 graph
    ############################
    plt.figure()
    #initial position of psips
    plt.scatter(initial_coor[:,0], initial_coor[:,3], c="gray", marker="+")
    #final position of psips
    plt.scatter(coor_psips[:,0],coor_psips[:,3], c="blue", alpha=0.15)
    #births and deaths
    plt.scatter(birth_array[:,0], birth_array[:,3], c="green", marker="*", alpha=0.5)
    plt.scatter(death_array[:,0], death_array[:,3], c="red", marker="+", alpha=0.5)
    #the place the psips should be - one e at one p, the other at the other
    plt.scatter(-d, d, color="black", marker="x")
    plt.scatter(d, -d, color="black", marker="x")
    #prettifying the graph
    plt.title("delta_t=" +str(delta_t)+", n="+str(n)+ " and #psips="+str(psips))
    plt.xlabel("x1")
    plt.ylabel("x2")
    fn = "protonpot_x1x2_delta_t=" +str(delta_t)+"_n="+str(n)+"_#psips="+str(psips)+"_d"+str(d)+"_it="+str(iteration)+".png"
    try:
        plt.savefig(fn)
    except:
        f=open(fn, "w")
        f.write ("Plot failed")
        f.close()
    plt.close()

    
    return(average_list)
