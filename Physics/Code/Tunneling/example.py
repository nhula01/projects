import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Function to create a Gaussian wave packet
def gaussian_wave_packet(x, mu, sigma, k):
    return np.exp(-0.5 * ((x - mu) / sigma)**2) * np.exp(1j * k * x)

# Parameters
x_min, x_max = -10, 10
num_points = 500
mu1, sigma1, k1 = -2, 1, 1
mu2, sigma2, k2 = 2, 1, -1
dt = 0.05

# Create the x values
x_values = np.linspace(x_min, x_max, num_points)

# Create the figure and axis
fig, ax = plt.subplots()

# Initialize the line objects for the two functions
line1, = ax.plot(x_values, np.real(gaussian_wave_packet(x_values, mu1, sigma1, k1)), label='Function 1')
line2, = ax.plot(x_values, np.real(gaussian_wave_packet(x_values, mu2, sigma2, k2)), label='Function 2')

# Add legend
ax.legend()

# Function to update the plot for each frame of the animation
def update(frame):
    t = frame * dt
    wave_packet1 = gaussian_wave_packet(x_values, mu1, sigma1, k1) * np.exp(-1j * k1 * t)
    wave_packet2 = gaussian_wave_packet(x_values, mu2, sigma2, k2) * np.exp(-1j * k2 * t)

    line1.set_ydata(np.real(wave_packet1))
    line2.set_ydata(np.real(wave_packet2))

    return line1, line2

# Create the animation
num_frames = 200
ani = FuncAnimation(fig, update, frames=num_frames, interval=50)

# Show the plot
plt.xlabel('Position')
plt.ylabel('Amplitude')
plt.title('Animated Gaussian Wave Packets')
plt.show()
