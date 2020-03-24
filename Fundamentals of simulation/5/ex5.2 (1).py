# 
"""
From Bhavya Joshi & Florian Jörg
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as mplt
import numpy as np
import pylab as plt
import math as m

# Testing the position: Only for smaller systems possible (for example: Ng=5)
def plotParticle(particle_x, particle_y, p, q, pp, qq, qstar, pstar, pf, qf):

    plt.close()
    
    plt.axis([0, 1, 0, 1])
    plt.plot(particle_x, particle_y ,'ro')
    plt.plot((p+0.5)*cellsize, (q+0.5)*cellsize, 'bo') 
    plt.plot((p+0.5)*cellsize, (qq+0.5)*cellsize, 'bo')      
    plt.plot((pp+0.5)*cellsize, (q+0.5)*cellsize, 'bo')  
    plt.plot((pp+0.5)*cellsize, (qq+0.5)*cellsize, 'bo') 
    plt.show()
    
    print 'Particle position: '+str(particle_x)+' '+str(particle_y)
    print 'p, q: '+str(p)+' '+str(q)    
    print 'pstar, qstar: '+str(pstar)+' '+str(qstar)    
    print 'pf, qf: '+str(pf)+' '+str(qf)

def threeDplot(Ng, field):
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    
    xpos = np.arange(0,Ng*Ng)
    ypos=[]
    for i in range(0,Ng):
        ypos.append(np.arange(0,Ng))
        
    xpos = np.repeat(np.arange(0,Ng), Ng)
    ypos=np.reshape(ypos, Ng*Ng)
    zpos = np.zeros(Ng*Ng)
    dx = np.ones(Ng*Ng)
    dy = np.ones(Ng*Ng)
    dz = np.arange(0,Ng*Ng)
    dz = np.reshape(field,Ng*Ng)
    #
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
    plt.show()
    

# Calculation of the force via CIC-Method

def add_particle_to_density_field(rho, Ng, cellsize, particle_x, particle_y, particle_mass):
    
    pf = particle_x / cellsize - 0.5
    qf = particle_y / cellsize - 0.5
    
    # Border conditions
    if (pf<0):
        pf=Ng+pf
    if (qf<0):
        qf=Ng+qf

    p=int(m.floor(pf))
    q=int(m.floor(qf))
    
    pstar = pf - p
    qstar = qf - q

    pp=p+1 
    qq=q+1
    
    # Border conditions
    if (pp>=Ng):
        pp=0
    if (qq>=Ng):
        qq=0
    
    rho[p][q] = (1-pstar)*(1-qstar)*particle_mass/cellsize/cellsize
    rho[pp][q] = pstar*(1-qstar)*particle_mass/cellsize/cellsize
    rho[p][qq] = (1-pstar)*qstar*particle_mass/cellsize/cellsize
    rho[pp][qq] = pstar*qstar*particle_mass/cellsize/cellsize
    
    return rho;
    

# Interpolate_force_field_to_particle_position via CIC-Interpolation
    
def ifftpp(force_field, Ng, cellcsize, particle_x, particl_y, particle_mass):
    pf = particle_x / cellsize - 0.5
    qf = particle_y / cellsize - 0.5

    # Border conditions
    if (pf<0):
        pf=Ng+pf
    if (qf<0):
        qf=Ng+qf

    p=int(m.floor(pf))
    q=int(m.floor(qf))
      
    pstar = pf - p
    qstar = qf - q
    
    pp=p+1 
    qq=q+1

    # Border conditions    
    if (pp>=Ng):
        pp=0
    if (qq>=Ng):
        qq=0
        
    a = (1-pstar)*(1-qstar)*particle_mass*force_field[p][q]
    a += pstar*(1-qstar)*particle_mass*force_field[pp][q]
    a += (1-pstar)*qstar*particle_mass*force_field[p][qq]
    a += pstar*qstar*particle_mass*force_field[pp][qq]
    
    return a


# Creating the vectors that are needed for the final plot
ax_vec=[]
ay_vec=[]
a_vec=[]
r_vec=[]
Ng = 256 # Dimension of field
L = 1.0 # Unit lenght
cellsize = L /Ng # cellsize
particle_mass = 1.0 # particle mass
G=1 # Green's function
    
# Repeat the calculation 10 times
for total in range (0,10):
    
    rho_real = [[0 for x in range(Ng)] for y in range(Ng)] # create empty density field
    
    # Calculate random particle positions
    particle_x = np.random.uniform(0,1)
    particle_y = np.random.uniform(0,1)

# a) 

    rho_real = add_particle_to_density_field(rho_real, Ng, cellsize, particle_x, particle_y, particle_mass)

    
# b) 

    rho_kspace = np.fft.fft2(rho_real)
    
    potential_kspace = [[0 for x in range(Ng)] for y in range(Ng)] # create density field
    
    for i in range(0,Ng):
        for j in range(0,Ng):
            if(i==0 and j==0):
                potential_kspace[0][0]=0
            else:
                potential_kspace[i][j]=-4.0*np.pi*G*L*L/(np.pi*np.pi*float(i*i+j*j))*rho_kspace[i][j]
    
    potential_rspace = np.fft.ifft2(potential_kspace)
    
# c) 

    
    ax = [[0 for x in range(Ng)] for y in range(Ng)] # create empty force field
    ay = [[0 for x in range(Ng)] for y in range(Ng)] # create empty force field
    
    for i in range(0,Ng):
        for j in range(0,Ng):
            if i==0:
                ax[i][j]=-(potential_rspace[i+1][j]-potential_rspace[Ng-1][j])/(2.*cellsize)
            elif (i==Ng-1):
                ax[i][j]=-(potential_rspace[0][j]-potential_rspace[i-1][j])/(2.*cellsize)
            else:
                ax[i][j]=-(potential_rspace[i+1][j]-potential_rspace[i-1][j])/(2.*cellsize)     
            
    
    for i in range(0,Ng):
        for j in range(0,Ng):
            if j==0:
                ay[i][j]=-(potential_rspace[i][j+1]-potential_rspace[i][Ng-1])/(2.*cellsize)
            elif (j==Ng-1):
                ay[i][j]=-(potential_rspace[i][0]-potential_rspace[i][j-1])/(2.*cellsize)
            else:
                ay[i][j]=-(potential_rspace[i][j+1]-potential_rspace[i][j-1])/(2.*cellsize)


# d) Calculate the forces inside the ring through interpolation

    rmin=0.3*L/Ng
    rmax=L/2.0
    
    for i in range(0,100):
        rp, rq = np.random.uniform(0,1), np.random.uniform(0,1)
        deltax=rmin*(rmax/rmin)**rp*m.cos(2*np.pi*rq)
        deltay=rmin*(rmax/rmin)**rp*m.sin(2*np.pi*rq)
        new_pos_x=particle_x+deltax
        new_pos_y=particle_y+deltay
    
        # Border conditions    
        if (new_pos_x>L):
            new_pos_x=new_pos_x-L
        elif (new_pos_x<0):
            new_pos_x=new_pos_x+L
        elif (new_pos_y>L):
            new_pos_y=new_pos_y-L
        elif (new_pos_y<0):
            new_pos_y=new_pos_y+L
        
        # Calculating the force-fields:
        r_vec.append(np.sqrt(deltax**2+deltay**2))
        ax_value=ifftpp(ax, Ng, cellsize, new_pos_x, new_pos_y, particle_mass)
        ay_value=ifftpp(ay, Ng, cellsize, new_pos_x, new_pos_y, particle_mass)
        ax_vec.append(ax_value)
        ay_vec.append(ay_value)
        a_vec.append(np.sqrt(ax_value**2+ay_value**2))


# e) Plot the result

a=[]
for i in range(0,1000):
    a.append(2./r_vec[i])

plt.close()
plt.plot(r_vec, a, label="a=2/r")
plt.axvline(1./Ng, color='k', label="r=1/Ng")
plt.loglog(r_vec, a_vec, 'ro')
plt.xlabel('r')
plt.ylabel('a')
plt.legend(loc=1)
plt.savefig('figures/ex5.2.eps', format='eps', dpi=1000)
plt.show()
