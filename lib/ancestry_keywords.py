"""
Ancestry keyword lists for classification.

Extracted from 01_create_mapping.py to reduce clutter and improve maintainability.
"""

# European subregions
KW_EUR_NW = [
    "british", "english", "irish", "welsh", "scottish", "scotland",
    "northern european", "north european", "northwestern european",
    "northwest european", "nw european", "scandinavian", "scandinavia",
    "scandavia", "swedish", "norwegian", "norweig", "danish", "dutch",
    "anglo", "celtic", "viking", "huguenot", "frisian",
    "french", "orcadian",
]

KW_EUR_S = [
    "south european", "southern european", "mediterranean",
    "italian", "iberian", "iberia", "maltese", "balkan",
    "sardinian", "portugese", "portuguese", "spanish", "greek", "romanian",
]

KW_EUR_E = [
    "eastern european", "east european", "north eastern european",
    "ukrainian", "ukranian", "polish", "latvian", "baltic",
    "slavic", "slovak", "sorbian", "czech", "hungarian",
    "finnish", "finn", "russian",
]

KW_EUR_C = [
    "german", "swiss", "austrian", "central european", "frankish", "bavaria",
]

KW_EUR_GENERAL = ["european", "europe"]

# Non-European populations
KW_AFRICAN = ["african", "sub-saharan", "west african"]

KW_EAST_ASIAN = [
    "east asian", "chinese", "korean", "japanese", "vietnamese",
    "siberian", "southeast asian",
    "asia",
]

KW_SOUTH_ASIAN = ["south asian", "sri lankan", "pashtun"]

KW_MIDDLE_EASTERN = [
    "middle eastern", "m. eastern", "arab", "iranian",
    "algerian", "western asian", "asia minor", "western asia",
]

KW_AMERICAS = [
    "native american", "amerindian", "leni lenape",
    "southern american", "panamanian",
]

KW_OCEANIAN = ["melanesia", "melanesian", "austronesian"]

# Founder populations
KW_JEWISH = ["ashkenazi", "sephardic", "jewish", "levite"]

# Special labels
SOCIAL_LABELS = {"caucasian"}
VAGUE_LABELS = {"mixed", "mixed ancestry"}
