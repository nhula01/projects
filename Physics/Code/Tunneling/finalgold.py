# import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.animation import FuncAnimation

def doublewell(x):
    alpha=.01
    V0=16
    potential=-1/np.sqrt((x-10)**2/2000+alpha)-1/np.sqrt(x**2/2000+alpha)+V0 
    return potential

def expectedPosition(name): # plot where the expected positions are
    position = np.loadtxt('/Users/phihung/NumMethod/first/final/x_vs_t'+name+'.dat')
    plt.plot(position[:,0], position[:,1],label='expected position')
    plt.xlabel('Time t [a.u.]')
    plt.ylabel('Position x [a.u.]')
    plt.title('Expected Position of the Wave Packet')
    plt.grid(True)
    plt.legend()
    plt.show()
    return None
expectedPosition('doublewell01')
def expectedMomentum(name): #  plot where the expected momentum is 
    momentum = np.loadtxt('/Users/phihung/NumMethod/first/final/p_vs_t'+name+'.dat')
    plt.plot(momentum[:,0], momentum[:,1]*2,label='expected momentum')
    plt.xlabel('Time t [a.u.]')
    plt.ylabel('Momentum p [a.u.]')
    plt.title('Expected Momentum of the Wave Packet')
    plt.grid(True)
    plt.legend()
    plt.show()
    return None
expectedMomentum('doublewell01')
def potential(x,V0=20): #plot the well
    x0=10.0
    width=7.5   
    potential=0.0
    if x>=(x0-width/2.0) and x<=(x0+width/2.0):
        potential=V0
    else:
        potential=0.0
    return potential
def time_dependent_potential(x, t):
    return 20 * abs(np.sin(2 * np.pi * .3 * t) * np.exp(-.2*(x-10)**2))
def scatteringPlot(name): #plot the heat map of the wave function
    density = '/Users/phihung/NumMethod/first/final/density_psi_xt'+name+'.dat'
    Z = np.loadtxt(density)
    xshape, tshape= Z.shape
    # two boundaries
    a = -80
    b = 80
    # total time propation
    tp = 20
    # generate the meshgrid
    dx = (b-a)/xshape
    x = [(a+dx*i) for i in range(xshape)]
    dt = tp/tshape
    t = [(0+dt*i) for i in range(tshape)]
    X, T = np.meshgrid(x, t) 
    Z = np.transpose(Z)
    # Create the contour plot
    contours = plt.contour(X, T, Z, levels=75)  # You can adjust the number of levels as needed
    plt.xlabel('x [a.u.]')
    plt.ylabel('t [a.u.]')
    plt.title('Heatmap of Density Wave')
    plt.colorbar(label='Probability Density') 
    #plot the barrier
    x = np.linspace(-70,70,200)
    y  = [doublewell(i) for i in x]
    plt.plot(x,y)
    plt.show()
    return None
scatteringPlot('doublewell01')
def coordinate_parts(file_name, y_shift,V,s=0): #animate the real and imaginary parts of space wavee function
    im_part = '/Users/phihung/NumMethod/first/final/im_psi_xt'+file_name+'.dat'
    real_part = '/Users/phihung/NumMethod/first/final/real_psi_xt'+file_name+'.dat'
    psi_im = np.loadtxt(im_part)
    psi_real = np.loadtxt(real_part)
    xshape, tshape = psi_im.shape
    # create the grids
    a = -80
    b = 80
    tp = 20
    dx = (b-a)/xshape
    x = [(a+dx*i) for i in range(xshape)]
    dt = tp/tshape
    t = [(0+dt*i) for i in range(tshape)]
    # Initialize the line object at t=0
    fig, ax = plt.subplots()
    line1, = ax.plot(x, psi_real[:,0],label=r'Re($\psi(x)$)')
    line3, = ax.plot(x, psi_im[:,0],label=r'Im($\psi(x)$)')
    # Plot the barrier
    x1 = np.linspace(-80,80,200)
    y1 = [doublewell(i) for i in x1]
    line2, = ax.plot(x1,y1,label='Barrier')
    ax.legend()
    #ax.set_ylim(2, 3)
    title = ax.set_title('')
    # Function to update the plot for each frame of the animation
    def update(frame):
        frame = frame + 1
        t = frame*dt
        line1.set_ydata(psi_real[:,frame]+y_shift)
        line3.set_ydata(psi_im[:,frame]+y_shift)
        line2.set_ydata(y1)
        title.set_text('Time = {:.2f}'.format(t))
        return line1, line2, line3, title
    # Create the animation
    num_frames = tshape-1
    ani = FuncAnimation(fig, update, frames=num_frames, interval=80)
    plt.xlabel('x [a.u.]')
    plt.ylabel('Energy [a.u.]')
    plt.title('Animated Gaussian Wave Packet')
    #saving
    if s!=0: #save if not 0
        #ax.set_ylim(2, 3)
        ani.save('/Users/phihung/NumMethod/first/final/rix_potential_'+file_name+'.gif', writer='pillow')
    # Show the plot
    plt.show()
    return None
coordinate_parts('doublewell01',2.5,4,0)

p=[]
def momentum_parts(file_name,s=0):
    im_part = '/Users/phihung/NumMethod/first/final/im_psi_pt'+file_name+'.dat'
    real_part = '/Users/phihung/NumMethod/first/final/real_psi_pt'+file_name+'.dat'
    psi_im = np.loadtxt(im_part)
    psi_real = np.loadtxt(real_part)
    #create grid for momentum and time
    pshape, tshape= psi_im.shape
    N=2048
    a=-80
    b=80
    dx = (b-a)/pshape
    p_b = 2.0*np.pi*(N/2-1)/(dx*N) 
    p_a = -2.0*np.pi*(N+1-(N/2+1))/(dx*N)
    tp = 20 # time propagating
    dp = (p_b-p_a)/pshape
    for i in range(N):
        if (i<=N/2):
            p.append(2*np.pi*(i-1)/(dx*N))
        else:
            p.append(-2*np.pi*(N+1-i)/(dx*N))
    dt = tp/tshape
    t = [(0+dt*i) for i in range(tshape)]
    # Initialize the wave at t=0 in momentum space
    fig, ax = plt.subplots()
    line1, = ax.plot(p, psi_real[:,0],label=r'Re($\psi(p)$)')
    line2, = ax.plot(p, psi_im[:,0],label=r'Im($\psi(p)$)')
    ax.legend()
    title = ax.set_title('')
    # Function to update the plot for each frame of the animation
    def update(frame):
        frame = frame + 1
        t = frame*dt
        line1.set_ydata(psi_real[:,frame])
        line2.set_ydata(psi_im[:,frame])
        #line2.set_ydata(time_dependent_potential(x1,t))
        title.set_text('Time = {:.2f}'.format(t))
        return line1, line2, title
    # Create the animation
    num_frames = tshape-1
    ani = FuncAnimation(fig, update, frames=num_frames, interval=80)
    plt.xlabel('p [a.u.]')
    plt.ylabel(r'$\psi(p)$')
    plt.title('Animated Gaussian Wave Packet')
    # save the gif
    if s!=0:
        ani.save('/Users/phihung/NumMethod/first/final/rip_potential_'+file_name+'.gif', writer='pillow')
    plt.show()
    return None
momentum_parts('doublewell01',0)

def coordinate_density(file_name, y_shift,V,s=0): #plot the density of coordinate space
    density = '/Users/phihung/NumMethod/first/final/density_psi_xt'+file_name+'.dat'
    psi_density = np.loadtxt(density)
    xshape, tshape= psi_density.shape
    #creeate space and time grid
    a = -80
    b = 80
    tp = 20
    dx = (b-a)/xshape
    x = [(a+dx*i) for i in range(xshape)]
    dt = tp/tshape
    # Initialize the line object at t=0
    fig, ax = plt.subplots()
    line1, = ax.plot(x, psi_density[:,0],label=r'$\|\psi(x)\|^2$')
    #barrier
    x1 = np.linspace(-80,80,200)
    y1 = [doublewell(i) for i in x1]
    line2, = ax.plot(x1,y1,label='Barrier')
    ax.legend()
    title = ax.set_title('')
    # Function to update the plot for each frame of the animation
    def update(frame):
        frame = frame + 1
        t = frame*dt
        line1.set_ydata(psi_density[:,frame]+y_shift)
        line2.set_ydata(y1)
        title.set_text('Time = {:.2f}'.format(t))
        return line1, line2, title
    # Create the animation
    num_frames = tshape-1
    ani = FuncAnimation(fig, update, frames=num_frames, interval=80)
    plt.xlabel('x [a.u.]')
    plt.ylabel('Energy [a.u.]')
    plt.title('Coordinate Space')
    #ax.set_ylim(2, 3)
    # save the plot
    if s!=0:
       # ax.set_ylim(2, 3)
        ani.save('/Users/phihung/NumMethod/first/final/dx_potential_'+file_name+'.gif', writer='pillow')
    plt.show()
    return None
coordinate_density('doublewell01',2.5,4,0)

p=[]
def momentum_density(file_name,s=0): # plot the density of the momentum space
    density = '/Users/phihung/NumMethod/first/final/density_psi_pt'+file_name+'.dat'
    m_density = np.loadtxt(density)
    pshape, tshape= m_density.shape
    #grid for momentum and space
    N=2048
    a=-80
    b=80
    dx = (b-a)/pshape #given px
    tp = 20 # time propagating
    for i in range(N):
        if (i<=N/2):
            p.append(2*np.pi*(i-1)/(dx*N))
        else:
            p.append(-2*np.pi*(N+1-i)/(dx*N))
    dt = tp/tshape
    t = [(0+dt*i) for i in range(tshape)]
    # Initialize the line object 
    fig, ax = plt.subplots()
    line1, = ax.plot(p, m_density[:,0],label=r'$\|\psi(p)\|^2$')
    title = ax.set_title('')
    # Function to update the plot for each frame of the animation
    def update(frame):
        frame = frame + 1
        t = frame*dt
        line1.set_ydata(m_density[:,frame])
        title.set_text('Time = {:.2f}'.format(t))
        return line1, title
    # Create the animation
    num_frames = tshape-1
    ax.legend()
    ani = FuncAnimation(fig, update, frames=num_frames, interval=80)
    plt.xlabel('p [a.u.]')
    plt.ylabel(r'$\|\psi(p)\|^2$')
    plt.title('Momentum Density')
    #plt.legend(True)
    if s!=0:
        ani.save('/Users/phihung/NumMethod/first/final/dp_potential_'+file_name+'.gif', writer='pillow')
    # Show the plot
    plt.show()
    return None
momentum_density('doublewell01',0)
#plt.show()