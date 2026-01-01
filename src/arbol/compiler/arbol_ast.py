from dataclasses import dataclass, field
from typing import List, Optional, Any

class Node:
    pass

@dataclass(frozen=True)
class Identifier(Node):
    name: str

@dataclass
class Parameter(Node):
    name: str
    type: str
    value: Optional[Any] = None

@dataclass
class Argument(Node):
    value: Any

@dataclass
class Program(Node):
    statements: List[Node]

@dataclass
class CircuitDefinition(Node):
    name: Identifier
    parameters: List[Parameter]
    body: List[Node]

@dataclass
class GateDefinition(Node):
    name: Identifier
    qubit_refs: List[Identifier]
    body: List[Node]

@dataclass
class QuantumOperation(Node):
    name: Identifier
    qubit_refs: List[Identifier]

@dataclass
class MeasureOperation(Node):
    qubit: Identifier
    classical_bit: Identifier

@dataclass
class ClassicalBitDefinition(Node):
    name: Identifier
    line: int
    column: int

@dataclass
class QubitDeclaration(Node):
    name: Identifier
    line: int = 0
    column: int = 0

@dataclass
class GateApplication(Node):
    gate_name: Identifier
    arguments: List[Argument]
    qubits: List[Identifier]
    line: int
    column: int

@dataclass
class StimulusDefinition(Node):
    name: Identifier
    promoter: 'Promoter'
    is_clause: 'IsClause'
    line: int
    column: int
    parameters: List[Parameter] = field(default_factory=list)

@dataclass
class Promoter(Node):
    type: str
    name: str

@dataclass
class IsClause(Node):
    promoters: List[Promoter]

@dataclass
class LogicalQubitDefinition(Node):
    name: Identifier
    size: int
    is_clause: 'IsClause'
    line: int
    column: int

@dataclass
class ConfigParameter(Node):
    name: Identifier
    value: Any

@dataclass
class StimulusApplication(Node):
    stimulus_name: Identifier
    target: Identifier
    arguments: List[ConfigParameter]

@dataclass
class MeasureStatement(Node):
    qubit: Identifier
    classical_bit: Optional[Identifier]
    line: int
    column: int

@dataclass
class ConfigDefinition(Node):
    sections: List['ConfigSection']

@dataclass
class ConfigSection(Node):
    name: Identifier
    parameters: List['ConfigParameter']

@dataclass
class RunCircuit(Node):
    circuit_name: Identifier
    arguments: List[Argument]
    line: int = 0
    column: int = 0

@dataclass
class GeneDefinition(Node):
    gene_name: Identifier
    properties: List[ConfigParameter]

@dataclass
class GenesBlock(Node):
    genes: List[GeneDefinition]

@dataclass
class GRNInteraction(Node):
    regulator: Identifier
    interaction_type: str
    target: Identifier
    params: List[ConfigParameter]

@dataclass
class GRNDefinition(Node):
    interactions: List[GRNInteraction]

@dataclass
class StepInstruction(Node):
    duration: Argument

@dataclass
class RunBlock(Node):
    instructions: List[Node]
