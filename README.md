# Monte Carlo Event Chain Algorithm
This method is used for the time evolution. It is called Event Chain Monte Carlo (ECMC) method and was developed from the Markov Chain Monte Carlo aiming at faster speed and better performance. In this method many disks are moved in each iteration and the total distance dislocated by all disks is the parameter `D`. For better performance, the directions of movement are `+x` and `+y`.


![alt text](/figures/fig2.gif "fig2.png")

## Installation
First clone this repository
```
git clone https://github.com/AndreyGFranca/mcec.git
```

Go to the folder and type
```
pip install ./mcec/
```

wait for the installation to finish. And you can run your simulation.


## Usage
```python
# Import mcec module
import mcec

# Make an instance of the system by choosing the number of particles and the density
sys = mcec.System(num_particles=8**2, density=5.0)

# Initialize the system
sys.init()

# Run the system n iterations, here we run for n = 100
sys.run(100)

```

![alt text](/figures/fig1.gif "Simulation 1")
![alt text](/figures/fig3.gif "Simulation 2")
