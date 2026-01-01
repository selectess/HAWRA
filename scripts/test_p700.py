
import sys
import os
import json
import numpy as np

# Add project root to path
sys.path.append('/Users/mehdiwhb/Desktop/HAWRA')
sys.path.append('/Users/mehdiwhb/Desktop/HAWRA/02_bioos_engine')

from simulations.multiphysics_simulator.simulator import MultiphysicsSimulator

config = {
    "max_time": 500,
    "dt": 5,
    "env": {
        "light_schedule": [
            {
                "start_time": 0,
                "end_time": 500,
                "intensity": 800
            }
        ]
    },
    "quantum": {
        "p700_threshold": 0.8,
        "decoherence_rate": 0.01
    },
    "bio": {
        "p700_synthesis_rate": 1.0,
        "p700_degradation_rate": 0.1
    }
}

sim = MultiphysicsSimulator(config)
log = sim.run()

print(f"Final P700 concentration: {log[-1]['p700_concentration']:.4f}")
