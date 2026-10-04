from __future__ import annotations

from dataclasses import dataclass
from typing import Any


class ChemistryError(RuntimeError):
    pass


@dataclass(frozen=True)
class ComponentCapability:
    name: str
    aliases: tuple[str, ...]
    insane_name: str | None
    martini_topology: bool
    mstool_mapping: bool
    charmm_topology: bool
    optional_plugin: str | None = None
    headgroup_beads: tuple[str, ...] = ()


COMPONENTS = {
    "POPC": ComponentCapability("POPC", (), "POPC", True, True, True, headgroup_beads=("PO4",)),
    "POPE": ComponentCapability("POPE", (), "POPE", True, True, True, headgroup_beads=("PO4",)),
    "POPS": ComponentCapability("POPS", (), "POPS", True, True, True, headgroup_beads=("PO4",)),
    "POPG": ComponentCapability("POPG", (), "POPG", True, True, True, headgroup_beads=("PO4",)),
    "DOPC": ComponentCapability("DOPC", (), "DOPC", True, True, True, headgroup_beads=("PO4",)),
    "DOPE": ComponentCapability("DOPE", (), "DOPE", True, True, True, headgroup_beads=("PO4",)),
    "DOPS": ComponentCapability("DOPS", (), "DOPS", True, True, True, headgroup_beads=("PO4",)),
    "DOPG": ComponentCapability("DOPG", (), "DOPG", True, True, True, headgroup_beads=("PO4",)),
    "CHOL": ComponentCapability("CHOL", ("CHL1",), "CHOL", True, True, True, headgroup_beads=("ROH",)),
    "DPG3": ComponentCapability("DPG3", ("GM3",), "DPG3", True, True, False, "glycolipid_template_finalize", ("GLC", "GAL", "NMC")),
}


ALIASES = {alias: name for name, item in COMPONENTS.items() for alias in item.aliases}


def component_for(name: str) -> ComponentCapability | None:
    key = str(name).upper()
    return COMPONENTS.get(ALIASES.get(key, key))


def validate_membrane_composition(config: dict[str, Any]) -> dict[str, Any]:
    membrane = config.get("membrane", {})
    composition = membrane.get("composition", {}) if isinstance(membrane, dict) else {}
    if not isinstance(composition, dict):
        raise ChemistryError("membrane.composition must be a mapping")
    units = composition.get("units", "count")
    if units not in {"count", "ratio", "fraction", "percent"}:
        raise ChemistryError("membrane.composition.units must be count, ratio, fraction, or percent")
    plugin_enabled = bool(
        config.get("stage3", {})
        .get("plugins", {})
        .get("glycolipid_template_finalize", {})
        .get("enabled", False)
    )
    require_backmap = bool(config.get("backmapping", {}).get("enabled", True))
    require_aa = bool(config.get("all_atom", {}).get("enabled", True))
    normalized: dict[str, dict[str, float]] = {"upper": {}, "lower": {}}
    for leaflet in ("upper", "lower"):
        values = composition.get(leaflet, {})
        if not isinstance(values, dict):
            raise ChemistryError(f"membrane.composition.{leaflet} must be a mapping")
        for raw_name, raw_value in values.items():
            try:
                value = float(raw_value)
            except (TypeError, ValueError) as exc:
                raise ChemistryError(f"membrane.composition.{leaflet}.{raw_name} must be numeric") from exc
            if value <= 0.0:
                raise ChemistryError(f"membrane.composition.{leaflet}.{raw_name} must be positive")
            component = component_for(str(raw_name))
            if component is None or component.insane_name is None:
                raise ChemistryError(f"unsupported membrane component for INSANE build: {raw_name}")
            if not component.martini_topology:
                raise ChemistryError(f"membrane component lacks Martini topology support: {raw_name}")
            if require_backmap and not component.mstool_mapping:
                raise ChemistryError(f"membrane component lacks mstool mapping support: {raw_name}")
            if require_aa and not component.charmm_topology:
                if component.optional_plugin and plugin_enabled:
                    pass
                else:
                    raise ChemistryError(
                        f"membrane component {raw_name} lacks all-atom topology support without "
                        f"{component.optional_plugin or 'an optional adapter'}"
                    )
            normalized[leaflet][component.insane_name] = normalized[leaflet].get(component.insane_name, 0.0) + value
    if not normalized["upper"] and not normalized["lower"]:
        raise ChemistryError("membrane.composition.upper/lower must contain at least one positive component")
    return {"status": "PASS", "units": units, "normalized_insane_composition": normalized}
