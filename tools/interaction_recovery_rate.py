import warnings

warnings.filterwarnings("ignore", category=UserWarning, module=r"prolif")
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")
warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"MDAnalysis(\.|$)")

from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

import argparse
import pandas as pd
from rdkit import Chem
import prolif as plf
from pathlib import Path

def parse_arguments():
    parser = argparse.ArgumentParser(description="Generate ProLIF fingerprints for crystal and ligand structures.")
    parser.add_argument("-p", "--protein", required=True, help="Path to protein PDB file")
    parser.add_argument("-x", "--xtal", required=True, help="Path to crystal ligand file (.sdf)")
    parser.add_argument("-l", "--ligand", required=True, help="Path to docked ligand file (.sdf)")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory for CSVs")
    return parser.parse_args()


def create_df(df_fingerprint):

    records = []
    df_t = df_fingerprint.T
    for i in range(df_t.shape[0]):
        count = int(df_t.iloc[i, 0])
        res = df_t.index[i][1] 
        interaction = df_t.index[i][2]  
        records.append({
            "count": count,
            "res": res,
            "type": interaction
        })
    return pd.DataFrame(records)


def load_ligand(path):
    path = Path(path)
    if path.suffix != ".sdf":
        raise ValueError(f"Unsupported file type: {path.suffix}")

    supplier = plf.sdf_supplier(str(path))
    for mol in supplier:
        if mol is not None:
            return mol

    raise ValueError(f"No valid molecule found in {path}")

    
def interaction_recovery(ref, lig):
    TS = []
    MS = []

    if ref.empty or lig.empty:
        return 0
    
    for i in range(ref.shape[0]):
        count_ref = ref.loc[i, 'count']
        res_ref = ref.loc[i, 'res']
        type_ref = ref.loc[i, 'type']

        match = lig[(lig['res'] == res_ref) & (lig['type'] == type_ref)]
        if match.empty or pd.isnull(match['count'].values[0]):
            TS.append(0)
        else:
            count_lig = match['count'].values[0]
            TS.append(min(count_lig, count_ref))
        MS.append(count_ref)
    if sum(MS) == 0:
        return 0
    else:
        return sum(TS)/sum(MS)


def main():
    args = parse_arguments()
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    protein_path = str(Path(args.protein).resolve())
    rdkit_prot = Chem.MolFromPDBFile(str(protein_path), removeHs=False, sanitize=False)
    protein_mol = plf.Molecule(rdkit_prot)

    xtal_mol = load_ligand(args.xtal)
    ligand_mol = load_ligand(args.ligand)

    fp_xtal = plf.Fingerprint(["HBAcceptor", "HBDonor", "Cationic", "Anionic", "PiStacking"], count=True)
    fp_lig = plf.Fingerprint(["HBAcceptor", "HBDonor", "Cationic", "Anionic", "PiStacking"], count=True)

    fp_xtal.run_from_iterable([xtal_mol], protein_mol)
    fp_lig.run_from_iterable([ligand_mol], protein_mol)

    records_xtal = create_df(fp_xtal.to_dataframe())
    records_xtal.to_csv(out_dir / "xtal_prolif.csv", index=False)

    records_lig = create_df(fp_lig.to_dataframe())
    records_lig.to_csv(out_dir / "ligand_prolif.csv", index=False)

    interaction_recovery_rate = interaction_recovery(records_xtal, records_lig)
    
    with open(out_dir / "result.txt", "w") as f:
        f.write(f"interaction_recovery_rate: {interaction_recovery_rate:}\n")

if __name__ == "__main__":
    main()
