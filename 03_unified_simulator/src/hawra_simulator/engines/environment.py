
class Environment:
    """
    Engine that simulates external environmental conditions.
    Manages temperature and a lighting schedule.
    """
    def __init__(self, config):
        """
        Initializes the environment with a configuration.

        Args:
            config (dict): Configuration dictionary.
        """
        self.config = config if config else {}
        self.temperature = self.config.get('temperature', 35.0)
        self.light_intensity = self.config.get('light_intensity', 0.0)
        self.light_wavelength = self.config.get('light_wavelength', None)
        schedule = self.config.get('light_schedule', [])
        normalized = []
        for ev in schedule:
            if isinstance(ev, dict) and 'time' in ev and 'intensity' in ev:
                item = {'time': float(ev['time']), 'intensity': float(ev['intensity'])}
                if 'wavelength' in ev:
                    item['wavelength'] = float(ev['wavelength'])
                normalized.append(item)
            elif isinstance(ev, (list, tuple)) and len(ev) == 2:
                normalized.append({'time': float(ev[0]), 'intensity': float(ev[1])})
        self.light_schedule = sorted(normalized, key=lambda x: x['time'])
        self._cursor = 0
        if self.light_schedule:
            self.update_from_schedule(0)

    def get_state(self):
        """
        Returns the current state of the environment.

        Returns:
            dict: A dictionary containing the environmental parameters.
        """
        return {
            'temperature': self.temperature,
            'light_intensity': self.light_intensity,
            'light_wavelength': self.light_wavelength
        }

    def set_temperature(self, value):
        """Sets the environment temperature."""
        self.temperature = value

    def set_light_intensity(self, value):
        """Sets the light intensity of the environment."""
        self.light_intensity = value

    def set_light_wavelength(self, value):
        self.light_wavelength = value

    def compact_schedule(self):
        if not self.light_schedule:
            return
        sorted_ev = sorted(self.light_schedule, key=lambda x: x['time'])
        compact = []
        last_intensity = None
        last_wavelength = None
        for ev in sorted_ev:
            intensity = ev['intensity']
            wavelength = ev.get('wavelength')
            if last_intensity is None or intensity != last_intensity or wavelength != last_wavelength:
                item = {'time': ev['time'], 'intensity': intensity}
                if wavelength is not None:
                    item['wavelength'] = wavelength
                compact.append(item)
                last_intensity = intensity
                last_wavelength = wavelength
        self.light_schedule = compact

    def update_from_schedule(self, time):
        if not self.light_schedule:
            return
        while self._cursor + 1 < len(self.light_schedule) and self.light_schedule[self._cursor + 1]['time'] <= time:
            self._cursor += 1
        current = self.light_schedule[self._cursor]
        self.light_intensity = current['intensity']
        self.light_wavelength = current.get('wavelength', self.light_wavelength)
