#!/usr/bin/env python3
"""Unified BSIM Schema Adapter for HAWRA project."""

import json
from typing import Dict, List, Any

class BSIMSchemaAdapter:
    """Unified adapter for BSIM schema conversion and validation."""
    
    def __init__(self):
        self.version = "0.3"
        
    def validate_bsim_schema(self, data: Dict[str, Any]) -> bool:
        """Validate that the data conforms to BSIM schema."""
        required_fields = ["metadata", "instructions"]
        if not all(field in data for field in required_fields):
            return False
            
        if "source_arbol" not in data["metadata"]:
            return False
            
        for instruction in data["instructions"]:
            if not all(field in instruction for field in ["instruction_id", "command", "params"]):
                return False
                
        return True
    
    def generate_demo_schema(self, demo_type: str = "bio_quantum") -> Dict[str, Any]:
        """Generate a demo BSIM schema for testing and demonstration."""
        if demo_type == "bio_quantum":
            return {
                "metadata": {
                    "source_arbol": "demo_bio_quantum.arbol",
                    "version": self.version,
                    "description": "Bio-quantum hybrid circuit with light stimulus"
                },
                "instructions": [
                    {
                        "instruction_id": 0,
                        "command": "INITIALIZE",
                        "params": {
                            "max_time": 500,
                            "dt": 5,
                            "env": {"light_schedule": []},
                            "quantum": {"p700_threshold": 0.8, "decoherence_rate": 0.05},
                            "bio": {
                                "p700_synthesis_rate": 1.0,
                                "p700_degradation_rate": 0.1,
                                "optimal_temp": 35.0,
                                "temp_sensitivity": 0.1,
                                "genes": {},
                                "grn": {}
                            }
                        }
                    },
                    {
                        "instruction_id": 1,
                        "command": "STIMULUS_APPLY",
                        "params": {
                            "stimulus": "light",
                            "target": "q1",
                            "arguments": {"wavelength": "660nm", "duration": "10ms"}
                        }
                    },
                    {
                        "instruction_id": 2,
                        "command": "QUANTUM_OP",
                        "params": {"gate": "H", "qubits": ["q1"]}
                    },
                    {
                        "instruction_id": 3,
                        "command": "MEASURE",
                        "params": {"qubit": "q1", "classical_bit": "c1"}
                    }
                ]
            }
        else:
            # Basic quantum circuit demo
            return {
                "metadata": {
                    "source_arbol": "demo_basic.arbol",
                    "version": self.version,
                    "description": "Basic quantum circuit with H gate and measurement"
                },
                "instructions": [
                    {
                        "instruction_id": 0,
                        "command": "INITIALIZE",
                        "params": {"qubits": 1, "classical_bits": 1}
                    },
                    {
                        "instruction_id": 1,
                        "command": "QUANTUM_OP",
                        "params": {"gate": "H", "qubits": ["q0"]}
                    },
                    {
                        "instruction_id": 2,
                        "command": "MEASURE",
                        "params": {"qubit": "q0", "classical_bit": "c0"}
                    }
                ]
            }

def main():
    """Demonstrate the BSIM schema adapter functionality."""
    adapter = BSIMSchemaAdapter()
    
    print("=== BSIM Schema Adapter Demo ===\n")
    
    # Generate a bio-quantum demo schema
    print("1. Generating bio-quantum demo schema...")
    bio_schema = adapter.generate_demo_schema("bio_quantum")
    print(f"   Generated {len(bio_schema['instructions'])} instructions")
    
    # Validate the generated schema
    print("2. Validating generated schema...")
    is_valid = adapter.validate_bsim_schema(bio_schema)
    print(f"   Schema is valid: {is_valid}")
    
    # Save to file
    print("3. Saving schema to file...")
    with open("demo_bio_quantum.bsim.json", "w") as f:
        json.dump(bio_schema, f, indent=2)
    print("   Saved to demo_bio_quantum.bsim.json")
    
    print("\n=== Demo completed successfully! ===")

if __name__ == "__main__":
    main()
