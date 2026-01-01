#!/usr/bin/env python3
"""
Isolation Environnementale HAWRA - Priorité 1
Calcule blindage Faraday + refroidissement modéré (système Peltier)
"""

import math
import json
from dataclasses import dataclass, asdict
from typing import Dict, Tuple

@dataclass
class FaradayShieldConfig:
    material: str
    conductivity: float  # S/m
    permeability: float  # H/m
    thickness_mm: float
    attenuation_db: float
    weight_kg: float

@dataclass
class CoolingConfig:
    target_temp_c: float
    ambient_temp_c: float
    heat_load_w: float
    peltier_cop: float
    required_power_w: float
    fan_airflow_cfm: float

@dataclass
class IsolationPlan:
    faraday_shield: FaradayShieldConfig
    cooling: CoolingConfig
    cost_estimate_usd: float
    feasibility: float

COPER_PERFORMANCE = {
    0.5: 1.5,
    0.75: 1.2,
    1.0: 0.9
}

MATERIAL_PROPERTIES = {
    "copper": {"conductivity": 5.96e7, "density": 8.96},
    "aluminum": {"conductivity": 3.5e7, "density": 2.70},
    "steel": {"conductivity": 6.0e6, "density": 7.85}
}

MU_0 = 4 * math.pi * 1e-7

class IsolationDesigner:
    def __init__(self, area_m2: float = 0.8, target_attenuation_db: float = 60):
        self.area_m2 = area_m2
        self.target_attenuation_db = target_attenuation_db

    def compute_faraday_shield(self, material: str = "copper", frequency_khz: float = 1.0) -> FaradayShieldConfig:
        props = MATERIAL_PROPERTIES[material]
        conductivity = props["conductivity"]
        permeability = MU_0
        skin_depth = math.sqrt(2 / (2 * math.pi * frequency_khz * 1e3 * conductivity * permeability))
        thickness = max(3 * skin_depth, 0.0005)
        attenuation = 8.686 * thickness / skin_depth
        if attenuation < self.target_attenuation_db:
            factor = self.target_attenuation_db / attenuation
            thickness *= factor
            attenuation *= factor
        volume = self.area_m2 * thickness
        weight = volume * props["density"]
        return FaradayShieldConfig(
            material=material,
            conductivity=conductivity,
            permeability=permeability,
            thickness_mm=thickness * 1000,
            attenuation_db=attenuation,
            weight_kg=weight
        )

    def compute_cooling(self, ambient_temp_c: float = 25.0, target_temp_c: float = 18.0, heat_load_w: float = 50.0) -> CoolingConfig:
        delta_t = ambient_temp_c - target_temp_c
        if delta_t <= 0:
            raise ValueError("Target temperature must be below ambient temperature")
        if delta_t <= 5:
            cop = 1.5
        elif delta_t <= 10:
            cop = 1.2
        else:
            cop = 0.9
        required_power = heat_load_w / cop
        airflow_cfm = (heat_load_w * 3.413) / (1.08 * delta_t)
        return CoolingConfig(
            target_temp_c=target_temp_c,
            ambient_temp_c=ambient_temp_c,
            heat_load_w=heat_load_w,
            peltier_cop=cop,
            required_power_w=required_power,
            fan_airflow_cfm=airflow_cfm
        )

    def estimate_cost(self, shield: FaradayShieldConfig, cooling: CoolingConfig) -> float:
        material_cost = {
            "copper": 20.0,
            "aluminum": 8.0,
            "steel": 5.0
        }
        area_mm2 = self.area_m2 * 1e6
        thickness_mm = shield.thickness_mm
        cost_material = material_cost[shield.material] * (area_mm2 * thickness_mm) / 1e6
        peltier_cost = 300.0
        fan_cost = 80.0
        control_cost = 150.0
        total_cost = cost_material + peltier_cost + fan_cost + control_cost
        return round(total_cost, 2)

    def create_plan(self) -> IsolationPlan:
        shield = self.compute_faraday_shield()
        cooling = self.compute_cooling()
        cost = self.estimate_cost(shield, cooling)
        feasibility = 0.9
        return IsolationPlan(
            faraday_shield=shield,
            cooling=cooling,
            cost_estimate_usd=cost,
            feasibility=feasibility
        )

    def to_json(self, plan: IsolationPlan) -> str:
        return json.dumps(asdict(plan), indent=2)


def main():
    designer = IsolationDesigner()
    plan = designer.create_plan()
    print("🛡️  PLAN ISOLATION HAWRA")
    print(designer.to_json(plan))


if __name__ == "__main__":
    main()
