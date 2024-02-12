#include "integrators.h"


void IntegratorEuler::step(ParticleSystem &system, double dt) {
    Vecd x0 = system.getState();
    Vecd dx = system.getDerivative();

    Vecd x1 = x0 + dt*dx;

    // Save state
    system.setState(x1);
}


void IntegratorSymplecticEuler::step(ParticleSystem &system, double dt) {
    Vecd x0 = system.getPositions();
    Vecd v0 = system.getVelocities();
    Vecd a0 = system.getAccelerations();

    Vecd v1 = v0 + a0*dt;
    Vecd x1 = x0 + v1*dt;

    // Since we handle everything separately, and not the state as a whole,
    // we need to manually set the new positions, velocities and updateForces
    // (which is otherwise called via setState if applyForces==true).
    system.setPositions(x1);
    system.setVelocities(v1);
    system.updateForces();
}


void IntegratorMidpoint::step(ParticleSystem &system, double dt) {
    Vecd x0  = system.getState();
    Vecd dx0 = system.getDerivative();

    // Use the midpoint method to estimate the state at the next time step
    Vecd x_mid = x0 + (dt / 2.0) * dx0;
    system.setState(x_mid);

    // Calculate derivative at the midpoint, and compute the new state with it
    Vecd dx_mid = system.getDerivative();
    Vecd x1 = x0 + dt * dx_mid;

    // Save state
    system.setState(x1);
}


void IntegratorVerlet::step(ParticleSystem &system, double dt) {
    Vecd x0     = system.getPositions();
    Vecd x_prev = system.getPreviousPositions();
    Vecd a0     = system.getAccelerations();

    // Initialize x_prev if it hasn't been initialized
//    static bool initialized = false;
//    if (!initialized) {
//        x_prev = x0 - system.getVelocities() * dt;
//        initialized = true;
//    }

    Vecd x1 = 2 * x0 - x_prev + a0 * dt * dt;
    Vecd v1 = (x1 - x0) / dt;

    // Since we handle everything separately, and not the state as a whole,
    // we need to manually set the new positions, velocities, and updateForces
    // (which is otherwise called via setState if applyForces == true).
    system.setPositions(x1);
    system.setPreviousPositions(x0);
    system.setVelocities(v1);
    system.updateForces();
}
// --------- ???? ----------


void IntegratorVelocityVerlet::step(ParticleSystem &system, double dt) {
    Vecd x0 = system.getPositions();
    Vecd v0 = system.getVelocities();
    Vecd a0 = system.getAccelerations(); // Save the initial accelerations for now

    // I handle the positions for now, since
    Vecd x1 = x0 + v0*dt + a0*dt*dt/2;
    system.setPositions(x1);
    // I update the forces, so I can get the next step's acceleration afterwards
    system.updateForces();

    // Get the next step's acceleration, and compute the Velocity Verlet formula
    Vecd a1 = system.getAccelerations();
    Vecd v1 = v0 + (a1 + a0)*dt/2;
    system.setVelocities(v1);
    system.updateForces();
}


void IntegratorRK2::step(ParticleSystem &system, double dt) {
    Vecd x0 = system.getState();
    Vecd k1 = system.getDerivative();

    // Compute a temporary x_temp, in order to store the state using k1
    Vecd x_temp = x0 + k1*dt;
    system.setState(x_temp);

    // Get the derivative of the above expression as k2,
    // and compute the formula of RK2: Xn+1 = Xn + (k1+k2)dt/2
    Vecd k2 = system.getDerivative();
    Vecd x1 = x0 + (k1+k2)*dt/2;

    // Save state
    system.setState(x1);
}


void IntegratorRK4::step(ParticleSystem &system, double dt) {
    Vecd x0 = system.getState();
    Vecd k1 = system.getDerivative();

    // First step: k1
    Vecd x_temp = x0 + k1*dt;
    system.setState(x_temp);

    // Second step: k2
    Vecd k2 = system.getDerivative();
    x_temp = x0 + k2*dt;
    system.setState(x_temp);

    // Third step: k3
    Vecd k3 = system.getDerivative();
    x_temp = x0 + k3*dt;
    system.setState(x_temp);

    // Fourth step: k4
    Vecd k4 = system.getDerivative();

    // Compute RK4 formula
    Vecd x1 = x0 + (k1 + 2*k2 + 2*k3 + k4) * dt/6;

    // Save state
    system.setState(x1);
}





