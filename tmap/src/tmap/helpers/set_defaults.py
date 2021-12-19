from typing import Dict, Any


def set_defaults(kws: Dict[str, Any], settings: Dict[str, Any]) -> Dict[str, Any]:
    for key, value in settings.items():
        if key not in kws:
            kws[key] = value
    return kws
