# Import necessary libraries
import numpy as np                    # For numerical operations and random number generation
import matplotlib.pyplot as plt       # For plotting and visualization
import rebound                        # N-body gravitational simulation library
from astropy import units as u        # For unit conversions (Earth mass to Solar mass)

# Initialize a new N-body simulation
sim = rebound.Simulation()

# Add the central star (primary body)
# Mass = 0.867 solar masses (similar to a K-type main sequence star)
# This becomes particle 0 and the gravitational center of the system
sim.add(m=0.867)

# Add a planet orbiting the star
sim.add(
    m=(1*u.Mearth).to(u.Msun).value,  # Mass: 1 Earth mass converted to solar masses (~3e-6 M_sun)
    a=0.7,                             # Semi-major axis: 0.7 AU (inside the habitable zone)
    e=0.05,                            # Eccentricity: 0.05 (nearly circular orbit, like Earth's)
    omega=np.random.uniform(0, 2*np.pi)  # Argument of periapsis: random angle (0 to 2π radians)
)

# Create and display an orbital plot showing the current configuration
# This visualizes the positions and orbital paths of all bodies in the simulation
op = rebound.OrbitPlot(sim)
plt.show()