"""
availability.py - Query the ionchannels data tables and report what
channel configurations, model types, species, and compartment maps
are defined for a given cell type or model name.

Usage:
    from cnmodel.data import report_available_configurations
    print(report_available_configurations(celltype='tstellate', species='mouse'))
"""

from . import _db


def report_available_configurations(celltype=None, species=None):
    """
    Scan the ionchannels data tables and report what is available.

    Parameters
    ----------
    celltype : str or None
        Optional filter string.  Any table whose name contains this string
        (case-insensitive) is included.  Pass None to list everything.
    species : str or None
        Optional filter: only show entries for this species.

    Returns
    -------
    str
        A human-readable multi-line summary.
    """
    # Collect {tablename: {dimension_name: set(values)}} from DATA
    tables = {}
    for k in _db.DATA.keys():
        # The table name is the sole string element in the key tuple.
        # Kwarg pairs are always 2-tuples; the table name is a bare string.
        # Two or more kwarg keys can sort before the table name, so we
        # cannot assume a fixed position — find the first string element.
        tablename = next((elem for elem in k if isinstance(elem, str)), None)
        if tablename is None:
            continue

        if tablename not in tables:
            tables[tablename] = {}
        for elem in k:
            if isinstance(elem, tuple):
                dim_name, dim_val = elem
                tables[tablename].setdefault(dim_name, set()).add(dim_val)

    # Separate channel tables from compartment tables


    channel_tables = {
        t: v for t, v in tables.items()
        if t.endswith('_channels') and not t.endswith('_compartments')
    }
    compartment_tables = {
        t: v for t, v in tables.items()
        if t.endswith('_channels_compartments')
    }

    # Apply celltype filter (matches table name substring, case-insensitive)
    if celltype is not None:
        filt = celltype.lower()
        channel_tables = {t: v for t, v in channel_tables.items()
                         if filt in t.lower()}
        compartment_tables = {t: v for t, v in compartment_tables.items()
                             if filt in t.lower()}

    # Apply species filter
    if species is not None:
        channel_tables = {
            t: v for t, v in channel_tables.items()
            if 'species' not in v or species in v.get('species', set())
        }
        compartment_tables = {
            t: v for t, v in compartment_tables.items()
            if 'species' not in v or species in v.get('species', set())
        }

    lines = []
    header = "Available ion channel configurations"
    if celltype:
        header += f" matching '{celltype}'"
    if species:
        header += f" for species='{species}'"
    lines.append(header + ":")

    if not channel_tables:
        lines.append("  (none found — check that ionchannels.py defines a table for this cell)")
        return "\n".join(lines)

    for tname in sorted(channel_tables):
        tinfo = channel_tables[tname]
        lines.append(f"\n  Channel table: {tname}")
        if 'species' in tinfo:
            lines.append(f"    species      : {sorted(tinfo['species'])}")
        if 'model_type' in tinfo:
            lines.append(f"    model_type   : {sorted(tinfo['model_type'])}")
        if 'field' in tinfo:
            # Show only gbar parameters to keep output concise
            gbars = sorted(f for f in tinfo['field'] if '_gbar' in f or f == 'units')
            if gbars:
                lines.append(f"    conductances : {gbars}")

        # Look for matching compartments table
        comp_name = tname + '_compartments'
        if comp_name in compartment_tables:
            ct = compartment_tables[comp_name]
            lines.append(f"    compartments table: {comp_name}")
            if 'model_type' in ct:
                lines.append(f"      model_type  : {sorted(ct['model_type'])}")
            if 'compartment' in ct:
                lines.append(f"      compartments: {sorted(ct['compartment'])}")
            if 'parameter' in ct:
                gbars = sorted(p for p in ct['parameter'] if '_gbar' in p)
                if gbars:
                    lines.append(f"      parameters  : {gbars}")
        else:
            lines.append(f"    compartments table: (none — uniform soma density used as fallback)")

    return "\n".join(lines)
