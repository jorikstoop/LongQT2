"""
Helper functions for postprocessing scripts

Author: Jorik Stoop
"""
from pathlib import Path
import re

# helper functions to extract info from filename
def get_cl(filename):
    """Extract cycle length from filename"""
    match = re.search(r"(?:^|_)CL_(\d+(?:\.\d+)?)(?=_|$)",
                      Path(filename).stem,
                     flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Missing CL value in filename: {filename}")
    return float(match.group(1))


def get_label(filename):
    """Extract model and treatment(s) from filename"""
    tokens = Path(filename).stem.upper().split("_")
    model = next((token for token in tokens
                  if token in {"WT", "T474I", "R56Q"}),
                  "")
    treatment_names = {
        "ISO": "Isoproterenol",
        "KSUPP": "K+_supplement",
        "VERA": "Verapamil",
        "NODRUG": "No_drug",
    }
    treatments = [treatment_names[token]
                  for token in tokens
                  if token in treatment_names]
    if not treatments:
        raise ValueError(f"Missing treatment tag in filename: {filename}")
    treatment = "+".join(treatments)
    return f"{model}—{treatment}" if model else treatment