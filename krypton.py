# Import necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import rebound
from astropy import units as u
import time
import sys

# Initialize a new N-body simulation
sim = rebound.Simulation()

# Set simulation units: time in years, distance in AU, mass in solar masses
sim.units = ('yr', 'AU', 'Msun')

# Add the central star (primary body)
# Mass = 0.867 solar masses (K-type main sequence star)
sim.add(m=0.867)

# Add Planet 1 (innermost)
sim.add(
    m=(0.6*u.Mearth).to(u.Msun).value,  # Mass: 0.6 Earth masses
    a=0.5,                               # Semi-major axis: 0.5 AU
    e=0.09,                              # Eccentricity: 0.09
    omega=5.37,                          
    inc=(-12*u.deg).to(u.rad).value,    # Inclination: -12 degrees
    Omega=2.14                           
)

# Add Planet 2 (inner habitable zone)
sim.add(
    m=(1*u.Mearth).to(u.Msun).value,    # Mass: 1 Earth mass
    a=0.7,                               # Semi-major axis: 0.7 AU
    e=0.05,                              # Eccentricity: 0.05 (nearly circular)
    omega=np.random.uniform(0, 2*np.pi) # Random argument of periapsis
)

# Add Planet 3 (middle planet)
sim.add(
    m=(3.2*u.Mearth).to(u.Msun).value,  # Mass: 3.2 Earth masses
    a=1.3,                               # Semi-major axis: 1.3 AU
    e=0.07,                              # Eccentricity: 0.07
    omega=3.92,                          
    inc=(8*u.deg).to(u.rad).value,      # Inclination: 8 degrees
    Omega=1.58                           
)

# Add Planet 4 (middle-outer)
sim.add(
    m=(5*u.Mearth).to(u.Msun).value,    # Mass: 5 Earth masses (super-Earth)
    a=2.1,                               # Semi-major axis: 2.1 AU
    e=0.1,                               # Eccentricity: 0.1
    omega=2.85,                          
    inc=(6.5*u.deg).to(u.rad).value,    # Inclination: 6.5 degrees
    Omega=4.05                           
)

# Add Planet 5 (outer)
sim.add(
    m=(95*u.Mearth).to(u.Msun).value,   # Mass: 95 Earth masses (Saturn-like)
    a=6.3,                               # Semi-major axis: 6.3 AU
    e=0.3,                               # Eccentricity: 0.3
    omega=5.86,                          
    inc=(21*u.deg).to(u.rad).value,     # Inclination: 21 degrees
    Omega=5.72                           
)

# Calculate Hill radius for stable moon placement
def hill_radius(a, m_planet, m_star):
    """Calculate Hill radius for a planet"""
    return a * (m_planet / (3 * m_star))**(1/3)

# Add Moon 1 orbiting Planet 2 (index 2 in particles list)
# For stability: place moon at ~5% of Hill radius (more conservative)
planet2_hill = hill_radius(0.7, (1*u.Mearth).to(u.Msun).value, 0.867)
moon1_distance = 0.05 * planet2_hill  # 5% of Hill radius
sim.add(
    primary=sim.particles[2],            # Orbits around Planet 2 (1 M⊕)
    m=(7.35e22*u.kg).to(u.Msun).value,  # Mass: ~Luna mass
    a=moon1_distance,                    # Safe distance at 5% Hill radius
    e=0.001,                             # Very low eccentricity
    omega=2.35,                          
    inc=(0.1*u.deg).to(u.rad).value,    # Very low inclination
    Omega=4.87 
)
print(f"Planet 2 Hill radius: {planet2_hill*1.496e8:.0f} km, Moon 1 distance: {moon1_distance*1.496e8:.0f} km ({100*moon1_distance/planet2_hill:.1f}% Hill)")

# Add Moon 2 orbiting Planet 3 (index 3 in particles list)
planet3_hill = hill_radius(1.3, (3.2*u.Mearth).to(u.Msun).value, 0.867)
moon2_distance = 0.04 * planet3_hill  # 4% of Hill radius
sim.add(
    primary=sim.particles[3],            # Orbits around Planet 3 (3.2 M⊕)
    m=(1.48e23*u.kg).to(u.Msun).value,  # Mass: ~2x Luna mass
    a=moon2_distance,                    # Safe distance
    e=0.001,                             # Very low eccentricity
    omega=1.28,                          
    inc=(0.1*u.deg).to(u.rad).value,    # Very low inclination
    Omega=5.13 
)
print(f"Planet 3 Hill radius: {planet3_hill*1.496e8:.0f} km, Moon 2 distance: {moon2_distance*1.496e8:.0f} km ({100*moon2_distance/planet3_hill:.1f}% Hill)")

# Add Moon 3 orbiting Planet 4 (index 4 in particles list) - First moon
planet4_hill = hill_radius(2.1, (5*u.Mearth).to(u.Msun).value, 0.867)
moon3_distance = 0.03 * planet4_hill  # 3% of Hill radius
sim.add(
    primary=sim.particles[4],            # Orbits around Planet 4 (5 M⊕ super-Earth)
    m=(1.48e23*u.kg).to(u.Msun).value,  # Mass: ~2x Luna mass
    a=moon3_distance,                    # Safe distance
    e=0.001,                             # Very low eccentricity
    omega=1.47,                          
    inc=(0.05*u.deg).to(u.rad).value,   # Very low inclination
    Omega=3.91 
)
print(f"Planet 4 Hill radius: {planet4_hill*1.496e8:.0f} km, Moon 3 distance: {moon3_distance*1.496e8:.0f} km ({100*moon3_distance/planet4_hill:.1f}% Hill)")

# Add Moon 4 orbiting Planet 4 (index 4 in particles list) - Second moon
moon4_distance = 0.06 * planet4_hill  # 6% of Hill radius (2x Moon 3's distance)
sim.add(
    primary=sim.particles[4],            # Orbits around Planet 4 (5 M⊕ super-Earth)
    m=(8.93e22*u.kg).to(u.Msun).value,  # Mass: ~1.2x Luna mass
    a=moon4_distance,                    # Safe distance, 2x Moon 3's distance
    e=0.001,                             # Very low eccentricity
    omega=4.71,                          
    inc=(0.08*u.deg).to(u.rad).value,   # Very low inclination
    Omega=2.19 
)
print(f"Planet 4 Hill radius: {planet4_hill*1.496e8:.0f} km, Moon 4 distance: {moon4_distance*1.496e8:.0f} km ({100*moon4_distance/planet4_hill:.1f}% Hill)")

# Move to center of mass frame (important for stability)
sim.move_to_com()

# CRITICAL: Use IAS15 for moon systems - it's slow but accurate
# WHFast doesn't handle planet-moon hierarchies well
sim.integrator = "ias15"
sim.dt = 0.001  # Small timestep for accuracy

print("\n⚠️  Using IAS15 integrator for moon stability (slower but accurate)")
print("    Expect ~30-60 seconds for 10,000 year simulation\n")

# Display initial orbital configuration
fig = plt.figure(figsize=(18, 12))

# Plot 1: Full planetary system with proper scaling
ax1 = plt.subplot(2, 3, (1, 3))
colors_planets = ['#FF6B6B', '#4ECDC4', '#6BCB77', '#45B7D1', '#FFA07A']
for i, (planet_idx, color) in enumerate(zip([1, 2, 3, 4, 5], colors_planets)):
    p_orb = sim.particles[planet_idx].orbit(primary=sim.particles[0])
    
    # Calculate orbit points
    theta = np.linspace(0, 2*np.pi, 200)
    r = p_orb.a * (1 - p_orb.e**2) / (1 + p_orb.e * np.cos(theta))
    x = r * np.cos(theta + p_orb.omega)
    y = r * np.sin(theta + p_orb.omega)
    
    ax1.plot(x, y, color=color, linewidth=2, label=f'Planet {i+1}')
    
    # Mark current position
    ax1.plot(sim.particles[planet_idx].x, sim.particles[planet_idx].y, 
             'o', color=color, markersize=8)

# Mark the star
ax1.plot(0, 0, '*', color='yellow', markersize=20, markeredgecolor='orange', markeredgewidth=1.5)

ax1.set_xlabel('x [AU]', fontsize=12, fontweight='bold')
ax1.set_ylabel('y [AU]', fontsize=12, fontweight='bold')
ax1.set_title('5-Planet System Configuration', fontsize=14, fontweight='bold')
ax1.set_aspect('equal')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right')
ax1.set_xlim(-8, 8)
ax1.set_ylim(-8, 8)

# Plot 2: Moon 1 around Planet 2
ax2 = plt.subplot(2, 3, 4)
moon1_orb = sim.particles[6].orbit(primary=sim.particles[2])
theta = np.linspace(0, 2*np.pi, 200)
r = moon1_orb.a * (1 - moon1_orb.e**2) / (1 + moon1_orb.e * np.cos(theta))
x = r * np.cos(theta + moon1_orb.omega) * 1.496e8  # Convert AU to km
y = r * np.sin(theta + moon1_orb.omega) * 1.496e8

ax2.plot(x, y, color='#95E1D3', linewidth=2)
ax2.plot(0, 0, 'o', color='#4ECDC4', markersize=10, label='Planet 2')
# Current position of moon relative to planet
moon_rel_x = (sim.particles[6].x - sim.particles[2].x) * 1.496e8
moon_rel_y = (sim.particles[6].y - sim.particles[2].y) * 1.496e8
ax2.plot(moon_rel_x, moon_rel_y, 'o', color='#95E1D3', markersize=6)
ax2.set_xlabel('x [km]', fontsize=10)
ax2.set_ylabel('y [km]', fontsize=10)
ax2.set_title('Moon 1 → Planet 2 (1 M⊕)', fontsize=11, fontweight='bold')
ax2.set_aspect('equal')
ax2.grid(True, alpha=0.3)

# Plot 3: Moon 2 around Planet 3
ax3 = plt.subplot(2, 3, 5)
moon2_orb = sim.particles[7].orbit(primary=sim.particles[3])
r = moon2_orb.a * (1 - moon2_orb.e**2) / (1 + moon2_orb.e * np.cos(theta))
x = r * np.cos(theta + moon2_orb.omega) * 1.496e8
y = r * np.sin(theta + moon2_orb.omega) * 1.496e8

ax3.plot(x, y, color='#FFD93D', linewidth=2)
ax3.plot(0, 0, 'o', color='#6BCB77', markersize=10, label='Planet 3')
moon_rel_x = (sim.particles[7].x - sim.particles[3].x) * 1.496e8
moon_rel_y = (sim.particles[7].y - sim.particles[3].y) * 1.496e8
ax3.plot(moon_rel_x, moon_rel_y, 'o', color='#FFD93D', markersize=6)
ax3.set_xlabel('x [km]', fontsize=10)
ax3.set_ylabel('y [km]', fontsize=10)
ax3.set_title('Moon 2 → Planet 3 (3.2 M⊕)', fontsize=11, fontweight='bold')
ax3.set_aspect('equal')
ax3.grid(True, alpha=0.3)

# Plot 4: Moons 3 & 4 around Planet 4
ax4 = plt.subplot(2, 3, 6)
moon3_orb = sim.particles[8].orbit(primary=sim.particles[4])
moon4_orb = sim.particles[9].orbit(primary=sim.particles[4])

# Moon 3 orbit
r = moon3_orb.a * (1 - moon3_orb.e**2) / (1 + moon3_orb.e * np.cos(theta))
x = r * np.cos(theta + moon3_orb.omega) * 1.496e8
y = r * np.sin(theta + moon3_orb.omega) * 1.496e8
ax4.plot(x, y, color='#C77DFF', linewidth=2, label='Moon 3')

# Moon 4 orbit
r = moon4_orb.a * (1 - moon4_orb.e**2) / (1 + moon4_orb.e * np.cos(theta))
x = r * np.cos(theta + moon4_orb.omega) * 1.496e8
y = r * np.sin(theta + moon4_orb.omega) * 1.496e8
ax4.plot(x, y, color='#A8DADC', linewidth=2, label='Moon 4')

ax4.plot(0, 0, 'o', color='#45B7D1', markersize=12, label='Planet 4')

# Current positions
moon_rel_x = (sim.particles[8].x - sim.particles[4].x) * 1.496e8
moon_rel_y = (sim.particles[8].y - sim.particles[4].y) * 1.496e8
ax4.plot(moon_rel_x, moon_rel_y, 'o', color='#C77DFF', markersize=6)

moon_rel_x = (sim.particles[9].x - sim.particles[4].x) * 1.496e8
moon_rel_y = (sim.particles[9].y - sim.particles[4].y) * 1.496e8
ax4.plot(moon_rel_x, moon_rel_y, 'o', color='#A8DADC', markersize=6)

ax4.set_xlabel('x [km]', fontsize=10)
ax4.set_ylabel('y [km]', fontsize=10)
ax4.set_title('Moons 3 & 4 → Planet 4 (5 M⊕)', fontsize=11, fontweight='bold')
ax4.set_aspect('equal')
ax4.grid(True, alpha=0.3)
ax4.legend(loc='upper right', fontsize=8)

plt.tight_layout()
plt.show()

# Setup for time evolution tracking
p = sim.particles
n_bodies = len(p) - 1  # Exclude the star
simulation_time = 10000  # years
n_steps = 500  # Number of snapshots

times = np.linspace(0, simulation_time, n_steps)

# Arrays to store orbital elements over time
smas = np.zeros((n_bodies, n_steps))
eccs = np.zeros((n_bodies, n_steps))

# Store initial energy for stability check
E0 = sim.energy()

print(f"Simulating {simulation_time} years with {n_steps} snapshots...")
print(f"Initial energy: {E0:.6e}")
print(f"Total bodies: {sim.N} (1 star + 5 planets + 4 moons)")
print(f"Integrator: IAS15 (adaptive timestep)")
time0 = time.time()

# Run simulation and record orbital elements
for i, t in enumerate(times):
    sim.integrate(t)
    
    # Record orbital elements for all bodies
    for j in range(n_bodies):
        try:
            orbit = p[j+1].orbit(primary=p[0])
            smas[j, i] = orbit.a
            eccs[j, i] = orbit.e
        except:
            # If orbit calculation fails (e.g., ejection), use NaN
            smas[j, i] = np.nan
            eccs[j, i] = np.nan
    
    # Progress indicator with better updates
    if i % 25 == 0:  # Less frequent updates for IAS15
        progress = 100 * t / simulation_time
        elapsed = time.time() - time0
        if i > 0:
            eta = elapsed / i * (n_steps - i)
            sys.stdout.write(f"\r  Progress: {progress:.1f}% | Elapsed: {elapsed:.1f}s | ETA: {eta:.1f}s")
        else:
            sys.stdout.write(f"\r  Progress: {progress:.1f}%")
        sys.stdout.flush()

# Print newline after progress indicator
print()

# Energy conservation check
E_final = sim.energy()
dE = abs((E_final - E0) / E0)
print(f"\nSimulation completed in {time.time() - time0:.2f} seconds")
print(f"Final energy: {E_final:.6e}")
print(f"Relative energy error: {dE:.2e}")

if dE > 1e-8:
    print("⚠️  Warning: Energy conservation may be poor.")
elif dE > 1e-10:
    print("✓ Energy reasonably conserved.")
else:
    print("✓✓ Energy excellently conserved!")

# Create comprehensive plots
fig, axes = plt.subplots(2, 1, figsize=(14, 10))

# Plot 1: Semi-major axis evolution
ax = axes[0]
planet_labels = ['Planet 1 (0.6 M⊕, 0.5 AU)', 
                 'Planet 2 (1 M⊕, 0.7 AU)',
                 'Planet 3 (3.2 M⊕, 1.3 AU)',
                 'Planet 4 (5 M⊕, 2.1 AU)', 
                 'Planet 5 (95 M⊕, 6.3 AU)',
                 'Moon 1 → P2',
                 'Moon 2 → P3',
                 'Moon 3 → P4 (inner)',
                 'Moon 4 → P4 (outer)']
colors = ['#FF6B6B', '#4ECDC4', '#6BCB77', '#45B7D1', '#FFA07A', 
          '#95E1D3', '#FFD93D', '#C77DFF', '#A8DADC']

# Plot planets with solid lines
for j in range(5):  # First 5 are planets
    ax.plot(times, smas[j], label=planet_labels[j], color=colors[j], linewidth=2)

# Plot moons with dashed lines (converted to AU for comparison)
for j in range(5, n_bodies):  # Remaining are moons
    ax.plot(times, smas[j], label=planet_labels[j], color=colors[j], 
            linewidth=1.5, linestyle='--', alpha=0.8)

ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
ax.set_ylabel('Semi-major Axis (AU)', fontsize=12, fontweight='bold')
ax.set_title('Orbital Evolution: Semi-major Axes', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=9, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, simulation_time)

# Plot 2: Eccentricity evolution
ax = axes[1]

# Plot planets with solid lines
for j in range(5):  # First 5 are planets
    ax.plot(times, eccs[j], label=planet_labels[j], color=colors[j], linewidth=2)

# Plot moons with dashed lines
for j in range(5, n_bodies):  # Remaining are moons
    ax.plot(times, eccs[j], label=planet_labels[j], color=colors[j], 
            linewidth=1.5, linestyle='--', alpha=0.8)

ax.set_xlabel('Time (years)', fontsize=12, fontweight='bold')
ax.set_ylabel('Eccentricity', fontsize=12, fontweight='bold')
ax.set_title('Orbital Evolution: Eccentricities', fontsize=14, fontweight='bold')
ax.legend(loc='upper right', fontsize=9, ncol=2)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, simulation_time)
ax.set_ylim(0, 0.5)  # Limit y-axis to see details better

plt.tight_layout()
plt.show()

# Summary statistics
print("\n" + "="*70)
print("ORBITAL STABILITY SUMMARY")
print("="*70)
for j in range(n_bodies):
    if not np.isnan(smas[j, -1]):
        a_change = abs(smas[j, -1] - smas[j, 0]) / smas[j, 0] * 100
        e_change = abs(eccs[j, -1] - eccs[j, 0])
        stability = "✓ STABLE" if a_change < 1.0 and e_change < 0.05 else "⚠ UNSTABLE"
        print(f"{planet_labels[j]}:")
        print(f"  Δa/a: {a_change:.3f}%  |  Δe: {e_change:.5f}  |  {stability}")
    else:
        print(f"{planet_labels[j]}:")
        print(f"  ✗ EJECTED or orbit calculation failed")
print("="*70)
