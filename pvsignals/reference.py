"""Reference sets of positive and negative controls, outcome definitions, and drug-name matching.

Controls come from two published reference sets, as distributed in the OHDSI ``MethodEvaluation`` R package
(Apache-2.0):

- OMOP (Ryan et al., Drug Safety 2013; 36 Suppl 1:S33-47): 165 positive and 234 negative controls for acute liver
  injury, acute kidney injury, acute myocardial infarction and upper GI bleeding.
- EU-ADR (Coloma et al., Drug Safety 2013; 36:13-23): 43 positive and 50 negative controls for 10 outcomes (the
  OHDSI copy drops nimesulide, which is not in RxNorm).

The outcomes are clinical concepts, so each is mapped here to a narrow list of MedDRA preferred terms (PTs). These
lists were written for this project from the outcome descriptions; they are not MedDRA SMQs.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

RDA_URL = "https://raw.githubusercontent.com/OHDSI/MethodEvaluation/main/data/{name}.rda"
REFERENCE_SETS = {"omop": "omopReferenceSet", "euadr": "euadrReferenceSet"}

# Reference-set outcome name -> outcome key used in this project.
OUTCOME_KEYS = {
    "OMOP Acute Liver Failure 1": "acute_liver_injury",
    "OMOP Acute Renal Failure 1": "acute_kidney_injury",
    "OMOP Acute myocardial Infarction  1": "acute_myocardial_infarction",
    "HOI Upper GI #3": "upper_gi_bleeding",
    "Anaphylaxis #1": "anaphylaxis",
    "Cardiac Valve Fibrosis #1": "cardiac_valve_fibrosis",
    "Leukopenia Including Neutropenia and Agranulocytosis": "leukopenia_agranulocytosis",
    "OMOP Aplastic Anemia 1": "aplastic_anaemia",
    "Rhabdomyolysis #1": "rhabdomyolysis",
    "Stevens-Johnson Syndrome #1": "sjs_ten",
}

# Outcome key -> MedDRA preferred terms (as spelled in FAERS: upper case, British spelling).
OUTCOME_PTS: dict[str, list[str]] = {
    "acute_liver_injury": [
        "DRUG-INDUCED LIVER INJURY", "HEPATOTOXICITY", "LIVER INJURY", "HEPATOCELLULAR INJURY",
        "CHOLESTATIC LIVER INJURY", "MIXED LIVER INJURY", "ACUTE HEPATIC FAILURE", "HEPATIC FAILURE",
        "SUBACUTE HEPATIC FAILURE", "HEPATITIS", "HEPATITIS ACUTE", "HEPATITIS CHOLESTATIC", "HEPATITIS FULMINANT",
        "HEPATIC NECROSIS", "JAUNDICE", "JAUNDICE CHOLESTATIC",
    ],
    "acute_kidney_injury": [
        "ACUTE KIDNEY INJURY", "RENAL FAILURE ACUTE", "RENAL FAILURE", "RENAL INJURY", "NEPHROPATHY TOXIC",
        "RENAL TUBULAR NECROSIS", "PRERENAL FAILURE", "ANURIA",
    ],
    "acute_myocardial_infarction": ["ACUTE MYOCARDIAL INFARCTION", "MYOCARDIAL INFARCTION",
                                    "SILENT MYOCARDIAL INFARCTION"],
    "upper_gi_bleeding": [
        "UPPER GASTROINTESTINAL HAEMORRHAGE", "GASTROINTESTINAL HAEMORRHAGE", "GASTRIC HAEMORRHAGE",
        "GASTRIC ULCER HAEMORRHAGE", "DUODENAL ULCER HAEMORRHAGE", "PEPTIC ULCER HAEMORRHAGE",
        "OESOPHAGEAL HAEMORRHAGE", "DUODENAL HAEMORRHAGE", "HAEMATEMESIS", "MELAENA", "GASTRITIS HAEMORRHAGIC",
    ],
    "anaphylaxis": ["ANAPHYLACTIC REACTION", "ANAPHYLACTIC SHOCK", "ANAPHYLACTOID REACTION", "ANAPHYLACTOID SHOCK"],
    "cardiac_valve_fibrosis": [
        "CARDIAC VALVE FIBROSIS", "CARDIAC VALVE DISEASE", "MITRAL VALVE INCOMPETENCE", "AORTIC VALVE INCOMPETENCE",
        "TRICUSPID VALVE INCOMPETENCE",
    ],
    "leukopenia_agranulocytosis": [
        "AGRANULOCYTOSIS", "NEUTROPENIA", "FEBRILE NEUTROPENIA", "LEUKOPENIA", "GRANULOCYTOPENIA",
        "NEUTROPHIL COUNT DECREASED", "WHITE BLOOD CELL COUNT DECREASED",
    ],
    "aplastic_anaemia": ["APLASTIC ANAEMIA", "BONE MARROW FAILURE"],
    "rhabdomyolysis": ["RHABDOMYOLYSIS"],
    "sjs_ten": ["STEVENS-JOHNSON SYNDROME", "TOXIC EPIDERMAL NECROLYSIS", "SJS-TEN OVERLAP"],
}
OUTCOME_LABEL = {
    "acute_liver_injury": "Acute liver injury", "acute_kidney_injury": "Acute kidney injury",
    "acute_myocardial_infarction": "Acute myocardial infarction", "upper_gi_bleeding": "Upper GI bleeding",
    "anaphylaxis": "Anaphylaxis", "cardiac_valve_fibrosis": "Cardiac valve fibrosis",
    "leukopenia_agranulocytosis": "Leukopenia / agranulocytosis", "aplastic_anaemia": "Aplastic anaemia",
    "rhabdomyolysis": "Rhabdomyolysis", "sjs_ten": "SJS / TEN",
}

# Single-PT events for the textbook signals named in the project brief (rhabdomyolysis reuses the outcome above).
EVENT_PTS: dict[str, list[str]] = {
    "tendon_rupture": ["TENDON RUPTURE"],
    "depression": ["DEPRESSION"],
}


@dataclass(frozen=True)
class DrugSpec:
    """How one ingredient is found in FAERS active-substance names.

    ``pattern`` is matched against each component of a lower-cased active-substance term (combination products are
    written ``A\\B``); ``search`` is a token or phrase whose openFDA search returns every report that could contain
    a matching term (a superset; several items are OR-ed), used to restrict count queries to the reports of a
    batch of drugs.
    """

    name: str
    pattern: str
    search: tuple[str, ...]


# Reference-set names that differ from FAERS active-substance names, or whose default pattern would be ambiguous.
SPECIAL: dict[str, tuple[str, tuple[str, ...]]] = {
    "amylases": (r"^pancrelipase\b", ("pancrelipase",)),
    "lipase": (r"^pancrelipase\b", ("pancrelipase",)),
    "endopeptidases": (r"^pancrelipase\b", ("pancrelipase",)),
    "chlorazepate": (r"^clorazepate\b", ("clorazepate",)),
    "estrogens, conjugated (usp)": (r"^estrogens,? conjugated\b", ("conjugated",)),
    "factor viia": (r"factor viia\b", ("viia",)),
    "interferon beta-1a": (r"^interferon beta-1a\b", ("interferon",)),
    "isosorbide": (r"^isosorbide (mono|di)nitrate\b|^isosorbide$", ("isosorbide",)),
    "lithium carbonate": (r"^lithium carbonate\b", ("lithium",)),
    "lithium citrate": (r"^lithium citrate\b", ("lithium",)),
    "ferrous gluconate": (r"^ferrous gluconate\b", ("ferrous",)),
    "ferrous sulfate": (r"^ferrous sulfate\b", ("ferrous",)),
    "mefenamate": (r"^mefenamic acid\b", ("mefenamic",)),
    "olmesartan medoxomil": (r"^olmesartan\b", ("olmesartan",)),
    "penicillin v": (r"^penicillin v\b", ("penicillin",)),
    "potassium chloride": (r"^potassium chloride\b", ("potassium chloride",)),
    "regular insulin, human": (r"^insulin human( regular)?$", ("insulin human",)),
    "sodium phosphate, monobasic": (r"^sodium phosphate", ("sodium phosphate",)),
    "tetrahydrocannabinol": (r"^(dronabinol|tetrahydrocannabinol)\b", ("dronabinol",)),
    "thyroxine": (r"^levothyroxine\b", ("levothyroxine",)),
    "valproate": (r"^(valproic acid|valproate|divalproex)\b", ("valproic", "valproate", "divalproex")),
    "vitamin a": (r"^vitamin a\b|^retinol\b", ("vitamin", "retinol")),
    "darbepoetin alfa": (r"^darbepoetin alfa\b", ("darbepoetin",)),
    "epoetin alfa": (r"^epoetin alfa\b", ("epoetin",)),
    "cephalexin": (r"^cephalexin\b|^cefalexin\b", ("cephalexin", "cefalexin")),
}


def drug_spec(name: str) -> DrugSpec:
    """The matching rule for a reference-set drug name (default: the name at the start of a component).

    The spec's ``name`` is the canonical drug (see ``canonical_drug``), which keys every count table.
    """
    key = name.strip().lower()
    if key in SPECIAL:
        pattern, search = SPECIAL[key]
        return DrugSpec(canonical_drug(key), pattern, search)
    return DrugSpec(key, rf"^{re.escape(key)}\b", (key,))


def canonical_drug(name: str) -> str:
    """Drug identity after matching (the three OMOP pancreatic-enzyme entries are all pancrelipase)."""
    key = name.strip().lower()
    return "pancrelipase" if key in {"amylases", "lipase", "endopeptidases"} else key


def components(term: str) -> list[str]:
    """Lower-cased components of a FAERS active-substance term (``A\\B`` -> ``["a", "b"]``)."""
    return [c.strip() for c in term.lower().split("\\") if c.strip()]


def term_matches(term: str, spec: DrugSpec) -> bool:
    """True if any component of the active-substance term matches the drug's pattern."""
    rx = re.compile(spec.pattern)
    return any(rx.search(c) for c in components(term))


def load_reference(raw_dir: Path, fetch: bool = True) -> pd.DataFrame:
    """Both reference sets as one table: set, drug, canonical drug, outcome, ground truth, indication.

    The ``.rda`` files are downloaded once into ``raw_dir`` (reading them needs ``pyreadr``). Pairs that become
    identical after mapping (the OMOP pancreatic-enzyme entries) are kept once; the 7 pairs present in both sets are
    kept in each set, with ``in_both`` marking them.
    """
    import pyreadr  # optional dependency of the data step only
    import requests

    raw_dir = Path(raw_dir)
    raw_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    for key, name in REFERENCE_SETS.items():
        path = raw_dir / f"{name}.rda"
        if not path.exists():
            if not fetch:
                raise FileNotFoundError(path)
            r = requests.get(RDA_URL.format(name=name), timeout=60)
            r.raise_for_status()
            path.write_bytes(r.content)
        df = list(pyreadr.read_r(str(path)).values())[0]
        frames.append(pd.DataFrame({
            "set": key,
            "drug_name": df["exposureName"].astype(str).str.strip(),
            "outcome_name": df["outcomeName"].astype(str),
            "ground_truth": df["groundTruth"].astype(int),
            "indication": df["indicationName"].astype(str),
        }))
    ref = pd.concat(frames, ignore_index=True)
    ref["drug"] = ref["drug_name"].map(canonical_drug)
    ref["outcome"] = ref["outcome_name"].map(OUTCOME_KEYS)
    if ref["outcome"].isna().any():
        raise ValueError(f"unmapped outcomes: {sorted(ref.loc[ref.outcome.isna(), 'outcome_name'].unique())}")
    ref = ref.drop_duplicates(["set", "drug", "outcome"]).reset_index(drop=True)
    both = ref.groupby(["drug", "outcome"])["set"].transform("nunique") > 1
    ref["in_both"] = both
    conflict = ref[both].groupby(["drug", "outcome"])["ground_truth"].nunique()
    if (conflict > 1).any():
        raise ValueError(f"conflicting ground truth: {conflict[conflict > 1].index.tolist()}")
    return ref
