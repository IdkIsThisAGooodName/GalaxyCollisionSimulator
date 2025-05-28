
# 🌌 Stable Galaxy Collision Simulation

This project simulates the collision of two galaxies using a Barnes-Hut tree algorithm for efficient gravitational force calculations. Stars and gas clouds interact under gravity, and the simulation is visualized in real-time 3D using the `vedo` library.

---

## 🚀 Features

* ✅ **Barnes-Hut tree** implementation for fast N-body gravitational force calculations (O(n log n))
* 🌍 Simulation of **two colliding galaxies** with:

  * Stars
  * Gas clouds
  * Central massive cores
* ⚙️ **Velocity damping** for numerical stability and realistic dynamics
* 🌀 **Angular momentum-based adjustments** to promote spiral structure formation
* 🎨 Real-time 3D visualization with color-coded gas velocities for enhanced insights

---

## 🛠️ Installation

### 1. Clone the repository

```bash
git clone https://github.com/IdkIsThisAGooodName/GalaxyCollisionSimulator.git
cd GalaxyCollisionSimulator
```

### 2. Python Dependencies

Install the required Python packages using pip:

```bash
pip install vedo numpy matplotlib
```

### 3. System Dependencies for C++ Backend

The project includes a C++ acceleration backend for improved performance. To build this, you need:

* Python development headers (e.g., `python3-devel`)
* `pybind11` headers (e.g., `pybind11-devel`)
* A C++17 compatible compiler (e.g., `g++`)

On Fedora, install them via:

```bash
sudo dnf install python3-devel pybind11-devel gcc-c++
```

### 4. Build the C++ Extension

Compile the Barnes-Hut tree C++ module with:

```bash
g++ -O3 -Wall -shared -std=c++17 -fPIC \
    $(python3 -m pybind11 --includes) \
    bh_tree.cpp -o bh_tree$(python3-config --extension-suffix)
```

---

## 🚀 Usage

Run the simulation script:

```bash
python GalaxyCollisionSimV1.py
```

This will launch a 3D visualization window showing the colliding galaxies with color-coded gas cloud velocities.

---

## 🧩 How it works

* The Barnes-Hut algorithm approximates gravitational forces efficiently using a tree data structure, reducing the computational complexity compared to naive all-pairs calculations.
* Two galaxies are initialized with stars, gas clouds, and massive central cores.
* Velocity damping and angular momentum adjustments help maintain numerical stability and realistic spiral arms.
* The `vedo` library provides smooth, interactive 3D visualization.

---

## 📂 Project Structure

```
GalaxyCollisionSimulator/
│
├── GalaxyCollisionSimV1.py
├── LICENSE
├── README.MD
├── bh_tree.cpp
└── requirements.txt
```

---

## 📝 License

This project is licensed under the MIT License. See the LICENSE file for details.

---

## 🙏 Contributions & Feedback

Contributions, bug reports, and feature requests are welcome! Feel free to open an issue or submit a pull request.

---

## 📬 Contact

If you have questions or want to collaborate, reach out via GitHub or email.

---

Thank you for checking out Stable Galaxy Collision Simulation! 🚀🌌

---
