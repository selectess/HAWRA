class PhotosynthesisManager:
    def check(self):
        return {"status": "ok", "p700_sites": 1024}

class IonChannelManager:
    def measure_ca_wave(self):
        return {"freq_hz": 0.5, "amplitude": 1.0}

class EMInterface:
    def __init__(self):
        self.supported_freq_range_khz = (1, 100)  # Gamme de fréquences supportée

    def check_frequency(self, target_khz: float) -> bool:
        """Vérifie si la fréquence est dans la gamme biologiquement supportée."""
        min_freq, max_freq = self.supported_freq_range_khz
        return min_freq <= target_khz <= max_freq

class VesicleManager:
    def transport(self, cargo: str) -> dict:
        return {"transported": cargo, "status": "ok"}

class CircadianClock:
    def now_phase(self) -> str:
        return "day"  # simple stub

class MetabolicNetwork:
    def status(self) -> dict:
        return {"CAM_active": True, "ATP_level": 0.8}

from Bio import SeqIO

class BioOS:
    """
    ARBOL Bio-Operating System
    Interfaces avec systèmes biologiques natifs
    """
    def __init__(self):
        self.photosynthesis = PhotosynthesisManager()  # P700 qubits
        self.ion_channels = IonChannelManager()        # Ca²⁺ signaling
        self.em_interface = EMInterface()              # CRY2 interface
        self.vesicle_transport = VesicleManager()      # Traffic cellulaire
        self.circadian = CircadianClock()             # 24h cycles
        self.metabolic = MetabolicNetwork()           # CAM/C3 switch
        self.genome = None  # Pour stocker les données du génome

    def load_genome(self, genome_file: str):
        """Charge et analyse le génome à partir d'un fichier GenBank."""
        print(f"   -> Chargement du génome depuis: {genome_file}")
        try:
            self.genome = list(SeqIO.parse(genome_file, "genbank"))
            if not self.genome:
                print("   ⚠️  Avertissement: Le fichier du génome est vide ou n'a pas pu être analysé.")
                return

            total_len = sum(len(rec) for rec in self.genome)
            print(f"   ✅ Génome chargé avec succès: {len(self.genome)} enregistrements, {total_len:,} bp au total.")

        except FileNotFoundError:
            print(f"   ❌ Erreur: Fichier du génome non trouvé à l'emplacement: {genome_file}")
        except Exception as e:
            print(f"   ❌ Erreur inattendue lors du chargement du génome: {e}")

    def find_gene(self, gene_name: str):
        """Recherche un gène par son nom dans le génome chargé."""
        if not self.genome:
            return None

        for record in self.genome:
            for feature in record.features:
                if feature.type == "gene" and "gene" in feature.qualifiers:
                    if feature.qualifiers["gene"][0] == gene_name:
                        return feature
        return None

    def compile_to_bio_signals(self, compiled_instructions: dict) -> dict:
        """Traduit les instructions ARBOL compilées en signaux biologiques plus réalistes."""
        bio_signals = {"actions": []}

        # 1. CIBLAGE DE GÈNE
        if "gene_targets" in compiled_instructions:
            for target in compiled_instructions["gene_targets"]:
                gene_name = target["gene_name"]
                gene_feature = self.find_gene(gene_name)
                if gene_feature:
                    action = f"TARGET_GENE_SUCCESS_{gene_name}_AT_{gene_feature.location}"
                else:
                    action = f"TARGET_GENE_FAILED_{gene_name}_NOT_FOUND"
                bio_signals["actions"].append(action)

        # 2. Gérer les portes quantiques -> transport de vésicules
        if compiled_instructions.get("quantum_gates"):
            transport_action = self.vesicle_transport.transport(cargo="[P700_cofactors]")
            if transport_action["status"] == "ok":
                bio_signals["actions"].append(f"VESICLE_TRANSPORT_FOR_{transport_action['transported']}")

        # 3. Gérer les impulsions EM -> horloge biologique
        if "em_events" in compiled_instructions:
            for em_event in compiled_instructions["em_events"]:
                if self.em_interface.check_frequency(em_event["freq_khz"]):
                    bio_signals["actions"].append(f"ACTIVATE_CRY2_CLOCK_WITH_FREQ_{em_event['freq_khz']}KHZ")
                else:
                    bio_signals["actions"].append(f"WARN_UNSUPPORTED_EM_FREQ_{em_event['freq_khz']}KHZ")

        # 4. Gérer les flashs IR -> lecture
        if "ir_events" in compiled_instructions:
            for ir_event in compiled_instructions["ir_events"]:
                bio_signals["actions"].append(f"PERFORM_LUC_READOUT_AT_{ir_event['wavelength_nm']}nm")

        return bio_signals