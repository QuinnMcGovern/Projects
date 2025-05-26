"""I want to run a statistical mechanics simulation. Spacificly, we will be simulating the Ising Model which, it's a model
which tell us how ferromagnetic (and antiferromagnetic) matirials work on a micro-scale. Ferromagnetic matirials being opposided to dielectric matirials;
copper, aliminium ect, which have a simple linear equation governing the magnetic feilds they output. We will build a 2D latice of random spins
then using the Metropolis–Hastings algorithm."""

#First, we import needed libaries
import numpy as np #need to run the complex mathamatics required, as well as being usefull for random numbers
import matplotlib as plt
import matplotlib.pyplot as plt #imports all needed elements of plt as plt so we can use them freely to present out results
import scipy.constants as C #imports all the constants we will need
from math import e#just gives us the value for euler's number for later
from matplotlib import animation    #this gives us a nessassary function for when we animate, artist animation

#these constants will be needed later in our probability calculation
T = 100 #defines the inital tempature in Kelvin
#If we had no inital tempature the system would tend to lower endergy, and therefore alighnment of spins, which is too simple
B = 1 / (T * C.k) #gives us the value of Beta we will use in our calculations   #I cheaked and it works ... mild rounding errors


J = 1 #defines the extange constant for our system
h = 0 #defines the external magnaetic feild
M = 1 #defines the magnaetic Moment

#to start we need to define an inital condition of spins
N = 150 #Number of elements we want in our model
Sigma_Sheet = np.random.binomial (1,0.5,[N,N])  #sets up and inital spread of spin data
Sigma_Sheet = 2 * Sigma_Sheet - 1 #this restructures our data to be -1 and 1 only with a 50/ 50 distribution rather than 0 to 1


"""With this basic information set up we can now graph our inital conditions"""
plt.imshow (Sigma_Sheet, cmap = "tab20b") #this shows our data 
plt.title (f"Monte Carlo Ising Model at T = {T}; Inital Condition") #gives a title
plt.xlabel ("Position") #lables our axis
plt.ylabel ("Position")
plt.colorbar (ticks = [-1,1], label = "Spin Direction") #this is used to indecate the sign of our iron 
plt.show ()


"""Now that we have defined our nessassary constants and inital conditions, we need to find the Hamiltonian of our system.
We can define the whole Hamiltonian as H = -J sum(S0 * Sij) - sum (M * h * S0)."""

#We now will use this equation to find the Hamiltonian of one of the points in our system

"""We use the Metropolis algorithm to simulate the evolution of the model. We will simulate with open boundaries, which end 
at the edges of our latince, opposed to looping boundaries which have no edges."""

def Hu (Spin_Array,Zero_Coord:int,One_Coord:int, J:float,h:float, m:float) -> float:
    """This function is used to find the energy in our system caused by one cell. We take in the whole array which it is stored within and the 
    x and y coordiantes, the joining constant J, the external magnetic field h and mangnetic dipole moment m.
    
    Inputs:
    Spin_Array: Array like
    Zero_Coord: int
    One_Coord: int
    J: float
    h: float
    m: float
    
    Return: float"""
    Ham = 0.0
    if Zero_Coord > 0:
        Ham += Spin_Array [Zero_Coord - 1,One_Coord]
    if Zero_Coord < len (Spin_Array[:,0]):
        Ham += Spin_Array [Zero_Coord + 1,One_Coord]
    if One_Coord > 0:
        Ham += Spin_Array [Zero_Coord,One_Coord - 1]
    if One_Coord < len (Spin_Array[0]):
        Ham += Spin_Array [Zero_Coord,One_Coord + 1]
    Ham *= -J * Spin_Array [Zero_Coord,One_Coord]
    Ham -= m * h * Spin_Array[Zero_Coord,One_Coord]
    return Ham

def Hv (Spin_Array,Zero_Coord:int,One_Coord:int, J:float,h:float, m:float) -> float:
    """This function is used to find the energy in our system caused by one cell. We take in the whole array which it is stored within and the 
    x and y coordiantes, the joining constant J, the external magnetic field h and mangnetic dipole moment m. Now switches the sign of the central cell.
    
    Inputs:
    Spin_Array: Array like
    Zero_Coord: int
    One_Coord: int
    J: float
    h: float
    m: float
    
    Return: float"""
    
    Ham = 0.0
    if Zero_Coord > 0:
        Ham += Spin_Array [Zero_Coord - 1,One_Coord]
    if Zero_Coord < len (Spin_Array[:,0]):
        Ham += Spin_Array [Zero_Coord + 1,One_Coord]
    if One_Coord > 0:
        Ham += Spin_Array [Zero_Coord,One_Coord - 1]
    if One_Coord < len (Spin_Array[0]):
        Ham += Spin_Array [Zero_Coord,One_Coord + 1]
    Ham *= J * Spin_Array [Zero_Coord,One_Coord]
    Ham += m * h * Spin_Array[Zero_Coord,One_Coord]
    return Ham



def Two_Random_Numbers (N:int)->int:
    """This function generates two random numbers we can use for indeies later."""
    return np.random.randint (0, N,1,dtype = 'int64')[0], np.random.randint (0, N,1,dtype = 'int64')[0]

iteration = 2000 #this is the number of times we want to iterate our simulation

Big_Results_Array = np.zeros ([iteration,N,N])

for i in range (iteration):#we loop through the number of iterations we have
    Number_0, Number_1 = Two_Random_Numbers (N - 1) #we first select an arbitary point in our range
    Hu_Value = Hu (Sigma_Sheet,Number_0,Number_1,J,h,M) #we calculate the local Hamiltonian
    Hv_Value = Hv (Sigma_Sheet,Number_0,Number_1,J,h,M) #we calculate the same but with a flipped sign of the cell's orrentation
    if Hv_Value < Hu_Value: #if the new energy is less than the old one we keep the switched sign
        Sigma_Sheet [Number_0,Number_1] *= -1
    if Hv_Value >= Hu_Value: #if not then we switch signs only if the probability is under a certain threshold given by e^-B (H2-H1)
        prob = e**(-B * (Hv_Value - Hu_Value))
        Threshold = np.random.rand ()
        if prob > Threshold:
            Sigma_Sheet [Number_0,Number_1] *= -1
    Big_Results_Array [i] += Sigma_Sheet #we save our results for each index



# plt.title (f"Monte Carlo Ising Model at T = {T}; Final Conditions")
# plt.xlabel ("Position")
# plt.ylabel ("Position")
# plt.imshow (Sigma_Sheet, cmap = "tab20b")
# plt.colorbar (ticks = [-1,1], label = "Spin Direction")
# plt.show ()




#now we can animate
fig, ax = plt.subplots() #set a figure, dpi is pixcel count
writervideo = animation.PillowWriter (fps = 8) #this starts the writer running at fps = 5 (for quicker load times)
ims = []

#ax.set (xlabel = "Position",ylabel = "Position")


ims = [] #holds each frame of our annimation
for i in range(iteration): #we save a frame for each iterations
    im = plt.imshow(Big_Results_Array[i], animated=True,cmap = "tab20b") #prints the results
    txt = ax.text (5,44,f"Number of iterations: {i}") #prints an update on the number of iterations 
    ims.append([im,txt])#adds through results
plt.colorbar (ticks = [-1,1], label = "Spin Direction") #adds detail to our graph
plt.title (f"Monte Carlo Ising Model at T = {T}")
plt.xlabel ("Position")
plt.ylabel ("Position")




ani = animation.ArtistAnimation (fig, ims, interval=50, blit=True, repeat_delay=1000) #collects our ainimation

ani.save('Ising_Model_Simulation.gif',writer = 'Pillow', fps = 8)  #saves the gif
plt.show ()