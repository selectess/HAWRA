# Schéma d'API pour les Drivers Matériels HAWRA

**Version:** 1.0
**Statut:** Proposition

## 1. Introduction

Ce document définit le contrat d'interface programmatique (API) entre l'**Orchestrateur BioOS** et la couche de **drivers matériels**. L'Orchestrateur, en interprétant un contrat `.bsim`, émet des appels à cette API pour manipuler les actionneurs et lire les données des capteurs.

L'API est conçue pour être simple, robuste et agnostique quant à l'implémentation sous-jacente des drivers. Elle peut être exposée via un serveur local sur le Jetson (ex: gRPC, REST sur un socket Unix).

## 2. Principes Généraux

- **Format des Données**: Toutes les requêtes et réponses utilisent le format JSON.
- **Gestion des Erreurs**: Les réponses d'API incluent systématiquement un champ `status` (`"ok"` ou `"error"`) et un champ `message` en cas d'erreur.
- **Adressage**: Chaque périphérique est identifié par un `id` unique (ex: `"light_red_660nm"`, `"electrode_channel_1"`).

---

## 3. API des Actionneurs (BioOS → Matériel)

Ces points d'accès permettent à BioOS de contrôler l'environnement physique de la PQPE.

### 3.1. Contrôleur d'Éclairage

Permet de contrôler les sources lumineuses (LEDs, lasers).

- **Endpoint**: `/api/actuators/light/set`
- **Méthode**: `POST`
- **Requête**:
  ```json
  {
    "channel_id": "RED_660nm",
    "intensity": 0.85,
    "duration_ms": 2000
  }
  ```
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "message": "Channel RED_660nm activated at 85% for 2000ms."
  }
  ```

### 3.2. Contrôleur de Champ Électromagnétique

Permet de générer des champs EM pour les opérations de porte quantique.

- **Endpoint**: `/api/actuators/em_field/set`
- **Méthode**: `POST`
- **Requête**:
  ```json
  {
    "frequency_hz": 9800,
    "amplitude_v": 5.0,
    "waveform": "sine",
    "duration_ms": 50
  }
  ```
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "message": "EM field (9.8kHz, 5V sine) applied for 50ms."
  }
  ```

### 3.3. Contrôleur Microfluidique

Permet de gérer l'injection de nutriments ou de réactifs.

- **Endpoint**: `/api/actuators/fluidics/pump`
- **Méthode**: `POST`
- **Requête**:
  ```json
  {
    "pump_id": "nutrient_pump_A",
    "volume_ul": 100,
    "flow_rate_ul_s": 10
  }
  ```
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "message": "Pump nutrient_pump_A dispensed 100uL at 10uL/s."
  }
  ```

---

## 4. API des Capteurs (Matériel → BioOS)

Ces points d'accès permettent à BioOS de lire l'état de la PQPE.

### 4.1. Capteur d'Imagerie (Bioluminescence)

Permet de capturer les signaux de bioluminescence pour la mesure des qubits.

- **Endpoint**: `/api/sensors/camera/read`
- **Méthode**: `GET`
- **Requête (paramètres d'URL)**: `?camera_id=readout_cam&exposure_ms=1000`
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "camera_id": "readout_cam",
    "timestamp": "2025-11-09T14:30:01Z",
    "data": {
      "format": "base64_png",
      "value": "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNkYAAAAAYAAjCB0C8AAAAASUVORK5CYII="
    }
  }
  ```

### 4.2. Capteur d'Électrode

Permet de mesurer les potentiels de surface ou d'autres signaux électriques.

- **Endpoint**: `/api/sensors/electrode/read`
- **Méthode**: `GET`
- **Requête (paramètres d'URL)**: `?channel_id=electrode_A1`
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "channel_id": "electrode_A1",
    "timestamp": "2025-11-09T14:30:02Z",
    "data": {
      "unit": "millivolt",
      "value": -70.5
    }
  }
  ```

### 4.3. Spectromètre

Permet d'analyser le spectre d'émission ou d'absorption.

- **Endpoint**: `/api/sensors/spectrometer/read`
- **Méthode**: `GET`
- **Requête (paramètres d'URL)**: `?integration_time_ms=500`
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "timestamp": "2025-11-09T14:30:03Z",
    "data": {
      "unit": "nm_vs_intensity",
      "value": [
        [400, 0.1],
        [401, 0.12],
        [...],
        [750, 0.05]
      ]
    }
  }
  ```

---

## 5. Point d'Accès de Statut

Permet de surveiller l'état de santé de l'ensemble du matériel.

- **Endpoint**: `/api/status`
- **Méthode**: `GET`
- **Réponse (`status: "ok"`)**:
  ```json
  {
    "status": "ok",
    "devices": [
      {"id": "light_red_660nm", "status": "ready"},
      {"id": "readout_cam", "status": "ready"},
      {"id": "em_field_controller", "status": "error", "message": "Calibration needed."}
    ]
  }
  ```
