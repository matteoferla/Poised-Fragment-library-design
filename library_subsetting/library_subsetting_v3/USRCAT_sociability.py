"""
This submodule contains the functions to calculate the USR and CAT scores for the sociability of a molecule,
using Torch tensors.
As a result there are import safeguards.
"""

import torch

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