# Environment Engine for the Multiphysics Simulator

class EnvironmentEngine:
    def __init__(self, config):
        print("Initializing Environment Engine...")
        self.config = config
        self.pulse_configs = config.get('pulse_configs', [])
        self.light_intensity = 0

    def update(self, time, dt):
        # Default to off
        current_intensity = 0
        # Check for active pulses
        for pulse in self.pulse_configs:
            if pulse['start'] <= time < pulse['end']:
                current_intensity = pulse['intensity']
                break # First active pulse wins
        
        if current_intensity != self.light_intensity:
            print(f"EVENT: Light intensity changed to {current_intensity} at t={time}")
            self.light_intensity = current_intensity

        return {
            'light_intensity': self.light_intensity
        }
