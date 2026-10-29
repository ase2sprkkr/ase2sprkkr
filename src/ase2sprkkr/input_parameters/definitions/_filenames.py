"""Native SPR-KKR dataset prefixes shared by defaults and file manifests."""


def dataset_adsi(parameters):
    """Return native DATSET: DATASET with a nonempty ADSI joined by '_'."""
    control = parameters.CONTROL
    dataset = control.DATASET()
    if dataset is None:
        return None
    adsi = control.ADSI()
    return dataset + (f"_{adsi}" if adsi else "")
