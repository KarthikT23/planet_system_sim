# Import necessary libraries
import numpy as np                    # For numerical operations and random number generation
import matplotlib.pyplot as plt       # For plotting and visualization
import rebound                        # N-body gravitational simulation library
from astropy import units as u        # For unit conversions (Earth mass to Solar mass)
import time                           # For timing how long the simulation takes to run

# Initialize a new N-body simulation
sim = rebound.Simulation()

# Set simulation units: time in years, distance in AU, mass in solar masses
# This makes the gravitational constant G = 2π²/year² and simplifies calculations
sim.units = ('yr', 'AU', 'Msun')

# Add the central star (primary body)
# Mass = 0.867 solar masses (similar to a K-type main sequence star)
# This becomes particle 0 and the gravitational center of the system
sim.add(m=0.867)

# Add the first planet (inner planet)
sim.add(
    m=(1*u.Mearth).to(u.Msun).value,  # Mass: 1 Earth mass converted to solar masses (~3e-6 M_sun)
    a=0.7,                             # Semi-major axis: 0.7 AU (inside the habitable zone)
    e=0.05,                            # Eccentricity: 0.05 (nearly circular orbit, like Earth's)
    omega=np.random.uniform(0, 2*np.pi)  # Argument of periapsis: random angle (0 to 2π radians)
)

# Add the second planet (outer planet)
sim.add(
    m=(5*u.Mearth).to(u.Msun).value,  # Mass: 5 Earth masses (super-Earth)
    a=2.1,                             # Semi-major axis: 2.1 AU (asteroid belt region)
    e=0.1,                             # Eccentricity: 0.1 (slightly elliptical orbit)
    omega=2.85,                        # Argument of periapsis: 2.85 radians (~163 degrees)
    inc=(6.5*u.deg).to(u.rad).value,  # Inclination: 6.5 degrees converted to radians (tilted orbit)
    Omega=4.05                         # Longitude of ascending node: 4.05 radians (~232 degrees)
)

# Create and display an orbital plot showing the current configuration
# color=True: Use different colors for each orbit
# periastron=True: Mark the closest approach point to the star
op = rebound.OrbitPlot(sim, color=True, periastron=True)
plt.show()

# Create time array: 50 evenly spaced points from 0 to 10,000 years
# This simulates the system evolution over 10,000 years
times = np.linspace(0, 10000)

# Create empty array to store semi-major axes over time
# Shape: (2 planets, 50 time steps)
smas = np.empty((2, len(times)))

# Record the start time to measure simulation performance
time0 = time.time()

# Run simulation and record semi-major axis at each time step
for i, t in enumerate(times):
    sim.integrate(t)                  # Advance simulation to time t (years)
    p = sim.particles                 # Get all particles in the simulation
    smas[0, i] = p[1].a              # Record semi-major axis of planet 1 (inner)
    smas[1, i] = p[2].a              # Record semi-major axis of planet 2 (outer)

# Print the total time taken to run the simulation (in seconds)
print(time.time() - time0)

# Plot the evolution of semi-major axes over time
# .T transposes the array so each planet is plotted as a separate line
plt.plot(times, smas.T)
plt.xlabel('Time (years)')
plt.ylabel('Semi-major Axis (AU)')
plt.title('Orbital Evolution Over 10,000 Years')
plt.legend(['Inner Planet (1 M⊕)', 'Outer Planet (5 M⊕)'])
plt.grid(True, alpha=0.3)
plt.show()
