import numpy as np
import vedo
import matplotlib.cm as cm
from concurrent.futures import ThreadPoolExecutor

# --- Constants ---
G = 4.302e-3  # pc * (km/s)^2 / solar_mass
softening = 10.0
dt = 0.01
num_steps = 500

# --- Galaxy parameters ---
num_stars = 200
num_gas = 500
radius = 15
central_mass = 5e6
core_mass = 1e8
velocity_scale = 1.0  # Lowered for stability

# --- Galaxy Generators ---
def generate_exponential_disk(center, velocity, num, radius, direction=1):
    r = radius * np.sqrt(np.random.rand(num))  # exponential-like distribution
    theta = np.random.rand(num) * 2 * np.pi

    x = center[0] + r * np.cos(theta)
    y = center[1] + r * np.sin(theta)
    z = center[2] + np.random.randn(num) * (radius * 0.05)  # scale height restored

    halo_mass = 1e10
    r_core = 2.0
    v_max = 120
    v_circ = np.sqrt(G * (central_mass + halo_mass * r / (r + r_core)) / (r + softening))
    v_circ = np.minimum(v_circ, v_max)

    vx = velocity[0] - direction * np.sin(theta) * v_circ
    vy = velocity[1] + direction * np.cos(theta) * v_circ
    vz = velocity[2] + np.random.randn(num) * 0.01

    vx += np.random.randn(num) * 0.002 * v_circ
    vy += np.random.randn(num) * 0.002 * v_circ

    return np.stack([x, y, z, vx, vy, vz], axis=1)


def generate_core(center, velocity):
    return np.array([[*center, *velocity]], dtype=float)

# --- BH Tree ---
class BHNode:
    def __init__(self, center, half_size, indices):
        self.center = center
        self.half_size = half_size
        self.indices = np.array(indices, dtype=int)
        self.children = []
        self.mass = 0.0
        self.mass_center = np.zeros(3)
        self.is_leaf = True

    def subdivide(self, pos, masses, min_size=1e-2):
        if len(self.indices) <= 1 or self.half_size < min_size:
            self.mass = masses[self.indices].sum()
            if self.mass > 0:
                self.mass_center = np.average(pos[self.indices], axis=0, weights=masses[self.indices])
            else:
                self.mass_center = self.center
            return

        self.is_leaf = False
        offsets = np.array([[dx, dy, dz] for dx in (-0.5, 0.5)
                                       for dy in (-0.5, 0.5)
                                       for dz in (-0.5, 0.5)])
        hs = self.half_size / 2
        children_indices = [[] for _ in range(8)]

        for idx in self.indices:
            p = pos[idx]
            i = 0
            if p[0] > self.center[0]: i += 4
            if p[1] > self.center[1]: i += 2
            if p[2] > self.center[2]: i += 1
            children_indices[i].append(idx)

        for i in range(8):
            child = BHNode(self.center + offsets[i] * hs * 2, hs, children_indices[i])
            child.subdivide(pos, masses, min_size)
            self.children.append(child)

        self.mass = sum(c.mass for c in self.children)
        if self.mass > 0:
            self.mass_center = sum(c.mass * c.mass_center for c in self.children) / self.mass
        else:
            self.mass_center = self.center

def build_tree(pos, masses):
    center = np.mean(pos, axis=0)
    max_range = np.max(np.ptp(pos, axis=0)) / 2
    root = BHNode(center, max_range, np.arange(len(pos)))
    root.subdivide(pos, masses)
    return root

# Multithreaded acceleration calculation helper function
def acc_from_node(node, i, pos, masses, theta):
    if node.mass == 0 or (i in node.indices and node.is_leaf):
        return np.zeros(3)
    r = node.mass_center - pos[i]
    d = np.linalg.norm(r) + 1e-10
    if node.is_leaf or (node.half_size / d < theta):
        return G * node.mass * r / (d**3 + softening**3)
    return sum(acc_from_node(c, i, pos, masses, theta) for c in node.children)

def compute_acceleration_mt(pos, masses, theta=0.5, max_workers=8):
    acc = np.zeros_like(pos)
    root = build_tree(pos, masses)

    def compute_acc_i(i):
        return acc_from_node(root, i, pos, masses, theta)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(compute_acc_i, range(len(pos))))
    acc[:] = results
    return acc

# --- Initialization ---
g1_center = [-40, 0, 0]
g2_center = [40, 0, 0]
g1_vel = [5, 0, 0]
g2_vel = [-5, 0, 0]

stars1 = generate_exponential_disk(g1_center, g1_vel, num_stars, radius)
stars2 = generate_exponential_disk(g2_center, g2_vel, num_stars, radius, direction=-1)
gas1 = generate_exponential_disk(g1_center, g1_vel, num_gas, radius)
gas2 = generate_exponential_disk(g2_center, g2_vel, num_gas, radius, direction=-1)
core1 = generate_core(g1_center, g1_vel)
core2 = generate_core(g2_center, g2_vel)

# Scale velocities
for group in [stars1, stars2, gas1, gas2, core1, core2]:
    group[:, 3:6] *= velocity_scale

stars = np.concatenate([stars1, stars2, core1, core2])
gas = np.concatenate([gas1, gas2])
m_stars = np.ones(len(stars))
m_stars[-2:] = core_mass  # cores

m_gas = np.full(len(gas), 0.1)

bodies = np.concatenate([stars, gas])
masses = np.concatenate([m_stars, m_gas])
star_idx = np.arange(len(stars))
gas_idx = np.arange(len(stars), len(bodies))

# --- Simulation ---
def run_sim():
    global bodies
    plt = vedo.Plotter(bg='black', size=(1200, 900), title="Stable Galaxy Collision")

    star_pts = vedo.Points(bodies[star_idx, :3], r=3, c='white')
    gas_pts = vedo.Points(bodies[gas_idx, :3], r=2)

    plt.add(star_pts, gas_pts)

    jet = cm.get_cmap('jet')

    def update(event):
        global bodies
        pos = bodies[:, :3]
        vel = bodies[:, 3:6]

        acc = compute_acceleration_mt(pos, masses)
        vel += acc * dt
        pos += vel * dt

        # Damping for stabilization
        vel[star_idx] *= 0.995
        vel[gas_idx] *= 0.99

        # Angular momentum-based spiral encouragement
        L = np.cross(pos, vel)
        L_mean = np.mean(L[star_idx], axis=0)
        norm_L = np.linalg.norm(L_mean)
        if norm_L > 0:
            L_dir = L_mean / norm_L
            vel_correction = np.cross(L_dir, vel)
            vel += 0.01 * vel_correction

        bodies[:, :3] = pos
        bodies[:, 3:6] = vel

        star_pts.points = pos[star_idx]
        gas_pts.points = pos[gas_idx]

        vmag = np.linalg.norm(vel[gas_idx], axis=1)
        vmag_norm = (vmag - vmag.min()) / (np.ptp(vmag) + 1e-10)
        colors = (jet(vmag_norm)[:, :3] * 255).astype(np.uint8)
        gas_pts.pointcolors = colors

        plt.render()

    plt.add_callback('timer', update)
    plt.timer_callback('start', 10)
    plt.show(interactive=True)

if __name__ == "__main__":
    run_sim()
