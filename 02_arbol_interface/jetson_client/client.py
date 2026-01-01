from flask import Flask, request, jsonify
import time
import logging

app = Flask(__name__)
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("JetsonClient")

# État simulé du matériel
hardware_state = {
    "light": {"status": "ready", "channels": {}},
    "em_field": {"status": "ready", "current": None},
    "sensors": {"camera": "ready", "electrodes": "ready"}
}

@app.route('/api/status', methods=['GET'])
def get_status():
    return jsonify({
        "status": "ok",
        "devices": [
            {"id": "light_controller", "status": hardware_state["light"]["status"]},
            {"id": "em_field_controller", "status": hardware_state["em_field"]["status"]},
            {"id": "readout_cam", "status": hardware_state["sensors"]["camera"]}
        ]
    })

@app.route('/api/actuators/light/set', methods=['POST'])
def set_light():
    data = request.json
    channel_id = data.get('channel_id')
    intensity = data.get('intensity')
    duration_ms = data.get('duration_ms', 0)
    
    logger.info(f"Setting light {channel_id} to {intensity} for {duration_ms}ms")
    # Simulation de l'action matérielle (GPIO/PWM sur Jetson)
    # TODO: Implémenter Jetson.GPIO ici
    
    return jsonify({
        "status": "ok",
        "message": f"Channel {channel_id} activated at {intensity*100}% for {duration_ms}ms."
    })

@app.route('/api/actuators/em_field/set', methods=['POST'])
def set_em_field():
    data = request.json
    freq = data.get('frequency_hz')
    amp = data.get('amplitude_v')
    
    logger.info(f"Applying EM field: {freq}Hz, {amp}V")
    # Simulation de l'action matérielle (DAC sur Jetson)
    
    return jsonify({
        "status": "ok",
        "message": f"EM field ({freq}Hz, {amp}V) applied."
    })

@app.route('/api/sensors/electrode/read', methods=['GET'])
def read_electrode():
    channel_id = request.args.get('channel_id')
    # Simulation de lecture ADC
    value = -70.0 + (time.time() % 1.0) # Valeur simulée
    
    return jsonify({
        "status": "ok",
        "channel_id": channel_id,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "data": {"unit": "millivolt", "value": value}
    })

if __name__ == '__main__':
    logger.info("HAWRA Jetson Client starting on port 5001...")
    app.run(host='0.0.0.0', port=5001)
