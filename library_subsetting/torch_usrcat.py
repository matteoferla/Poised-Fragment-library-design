import torch
from typing import Dict, List
from library_classifier import Classifier, InchiType, RoboDecomposer
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def calc_usrscores(d1: torch.Tensor, d2_matrix: torch.Tensor) -> torch.Tensor:
    """
    This copies what the C++ code ``calcUSRScore`` did, but in vectorised form.
    The variables names d1, d2 come from the C++ code and stand for ...?
    Re the chunk of 12, I need to read the paper.

    .. code-block:: cpp

        double calcUSRScore(const std::vector<double> &d1,
                            const std::vector<double> &d2) {
          // This code is stolen from RDKit, but stripped of weights by MF
          unsigned int num = 12;  // length of each subset
          double score = 1.0;
          for (unsigned int w = 0; w < (d1.size() / num); ++w) {
            double tmpScore = 0.0;
            unsigned int offset = num * w;
            for (unsigned int i = 0; i < num; ++i) {
              tmpScore += fabs(d1[i + offset] - d2[i + offset]);
            }
            tmpScore /= num;
            score += tmpScore;
          }
          return 1.0 / score;
        }
    """
    num = 12  # length of each subset
    # Reshape d1 and d2 into chunks of size 12, because that is what `calcUSRScore` did
    d1_chunks = d1.view(-1, num)  # Shape: (num_chunks, 12)
    d2_chunks = d2_matrix.view(d2_matrix.size(0), -1, num)  # Shape: (num_rows, num_chunks, 12)
    # d1_chunks will be broadcast to match d2 to make the absolute differences
    # for each chunk between d1 and each row in d2_matrix work, this averaged across the extra dim (num_rows)
    # dim=1 is chunks. dim=0 (num_rows) should kept!
    protoscores = torch.abs(d2_chunks - d1_chunks).mean(dim=2).sum(dim=1)  # Shape: (num_rows)
    # Final scores (one for each row in d2_matrix)
    return 1.0 / (1.0 + protoscores)  # Shape: (num_rows)

def calc_summed_scores(d1: torch.Tensor, d2_matrix: torch.Tensor, d2_weights: torch.Tensor, cutoff=0.7) -> torch.Tensor:
    scores = calc_usrscores(d1, d2_matrix)  # Shape: (num_rows)
    return d2_weights[scores > cutoff].sum()  # Shape: 1

class GPUClassifier(Classifier):
    def __init__(self,
                 common_synthons_tally: Dict[InchiType, int],
                 common_synthons_usrcats: Dict[InchiType, list]):
        self.common_synthons_tally = torch.tensor(common_synthons_tally, device='cuda')
        self.common_synthons_usrcats = torch.tensor(common_synthons_usrcats, device='cuda')
        self.dejavu_synthons: Dict[InchiType, int] = {}
        self.robodecomposer = RoboDecomposer()

    def calc_sociability(self, synthon: Chem.Mol) -> float:
        synthon_inchi = Chem.MolToInchi(synthon)
        if synthon_inchi in self.dejavu_synthons:
            return self.dejavu_synthons[synthon_inchi]
        if synthon is None:
            return -1
        AllChem.EmbedMolecule(synthon)
        if Chem.Mol.GetNumHeavyAtoms(synthon) < 3 or Chem.Mol.GetNumConformers(synthon) == 0:
            return -1
        synthon_usrcat = torch.tensor(rdMolDescriptors.GetUSRCAT(synthon), device='cuda')
        sociability = calc_summed_scores(synthon_usrcat, self.common_synthons_usrcats, self.common_synthons_tally).tolist()
        self.dejavu_synthons[synthon_inchi] = sociability
        return sociability
    def calc_synthon_info(self, mol, verdict):
        synthons: List[Chem.Mol] = self.robodecomposer.decompose(mol)
        verdict['N_synthons'] = len(synthons)
        verdict['synthon_sociability'] = sum(
            [self.calc_sociability(synthon) for synthon in synthons])
