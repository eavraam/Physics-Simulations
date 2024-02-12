<div align="center">

  <h3 align="center">Physics Simulations</h3>

  <div align="center">
    An introduction to physics simulations in C++/Qt, from numerical integrators to collisions, cloth and fluid simulation.
  </div>
  
  <br />

  <div>
    <img src="https://img.shields.io/badge/OpenGL-black?style=flat" alt="opengl" />
    <img src="https://img.shields.io/badge/C%2B%2B-blue?style=flat" alt="cpp" />
    <img src="https://img.shields.io/badge/QtCreator%205-green?style=flat" alt="qtcreator5" />
    <img src="https://img.shields.io/badge/Ubuntu-orange?style=flat" alt="ubuntu" />
  </div>
</div>



## 📋 <a name="table">Table of Contents</a>
1. 🤖 [Introduction](#introduction)
2. ⚙️  [How to run](#how-to-run)
3. 📐 [Integrators](#integrators)
4. 💥 [Collisions](#collisions)
5. 👗 [Cloth Simulation](#cloth) 
6. 🌊 [Fluid Simulation](#fluid)


## <a name="introduction">🤖 Introduction</a>
Developed using C++ inside QtCreator 5 and OpenGL as the rendering backend, as part of the Computer Animation (CA) course. This project's goal is to teach students the core concepts in physics simulations, like particle systems, numerical integrators and collisions, as well as their application in more complex topics, like cloth and fluid simulation.


## <a name="how-to-run">⚙️ How to run</a>

**Prerequisites**

Make sure you have the following installed on your machine:
- [Ubuntu Linux](https://ubuntu.com/) (Operating System)
- [QtCreator 5](https://doc.qt.io/qt-5/gettingstarted.html)

**Running the project**
- Navigate to the `/Lab-code/` directory.
- Open Simulations.pro
- Compile & run


## <a name="integrators">📐 Integrators</a>

In the first step of the assignment, I implemented various different numerical solvers in order to properly compute the movement of the particles. At this stage, some basic (but important) forces were added as options for the whole particle system, namely: constant acceleration force (gravity), linear and quadratic drag force. On a later stage, I also added a black hole (gravitational attraction force).

**Provided Numerical Solvers**

- Analytical motion
- Euler

**Implemented Numerical Solvers**

- Symnplectic Euler
- Verlet
- Midpoint
- Runge-Kutta 2nd order (RK-2)
- Runge-Kutta 4th order (RK-4)
- Velocity Verlet

Projectiles - Integrators Simulation|
:-------------------------:|
![](repo_images/sim_projectiles.webm)|


## <a name="collisions">💥 Collisions</a>

The next goal is to render thousands of particles and properly test/resolve collisions
with primitive shapes. Here, the particle simulation is a fountain which originates from a single point
and spreads upwards.

The tasks for this exercise are:

- Collisions: Test and Response (Plane, AABB, Sphere)
- Gravity + another force (Drag, Gravitational attraction - black hole, etc.)
- Mouse interaction: e.g. modify / move colliders
- Particle-particle collision (Spatial Hash Tables)

To play around with the application, you could use the following key combinations:

- CTRL + Right click: Move the fountain's origin.
- Shift + Right click (on a primitive shape): Move the shape around

**<ins>Disclaimer</ins>:** The implementation of Spatial Hash Tables is based on Matthias Muller's approach, found [here](https://matthias-research.github.io/pages/tenMinutePhysics/15-selfCollision.pdf).

**<ins>Note</ins>:** For the best particle-AABB collision testing and resolution, I tried implementing `Continuous-Collision Detection (CCD)`. However, there appears to be some issue that I could not resolve at the time, meaning that the particle-AABB collisions are not totally accurate.

Collisions - Fountain Simulation|
:-------------------------:|
![](repo_images/sim_fountain.webm)|


## <a name="cloth">👗 Cloth Simulation</a>

**Provot's Spring Model**

The goal of the current exercise is to implement `Provot's spring model` for cloth simulation. The tasks for this exercise are:

- Create the cloth as a grid of particles, which consider and apply forces based on Provot's spring
model. (**DONE**)
- Particle-collider test/resolution (**DONE**)
- Bonus: Cloth surface - Collider test/resolution (**NOT-DONE**)
- Bonus++: Prevent self-intersections (**APPROXIMATED**)

For this simulation, the Verlet integrator was used, since this provided the more stable and consistent result. Prevention of self-intersection with the use of spatial hashing was attempted, and seems to
provide some improvement, but it doesn't look to work exactly as expected. An additional checkbox
was added in the UI to enable/disable this spatial hashing functionality while the application is
running.

**<ins>Note</ins>:** To move the Sphere prop object that appears in the scene, you can press CTRL + Right click
while targeting it (via ray-casting) to move it around. By selecting a particle of the cloth using
Right click, you can "pinch and hold" the cloth from this particle and move the cloth around. While
holding the particle with Right click, if you press the button "F" you can fix this particle as an anchor
to make it static.

**<ins>Note 2</ins>:** To run the simulation, you also need to click the `Update` button.

Cloth Simulation|
:-------------------------:|
![](repo_images/sim_cloth.webm)|



## <a name="fluid">🌊 Fluid Simulation</a>