"""Version of the remote helper's control commands and JSON payloads."""

CONTROL_PROTOCOL_VERSION = 1


def check_control_protocol(version):
    """Reject an incompatible client before performing a remote operation."""
    if type(version) is not int or version != CONTROL_PROTOCOL_VERSION:
        raise ValueError(
            f"Unsupported remote control protocol {version!r}; "
            f"expected {CONTROL_PROTOCOL_VERSION}. "
            "Use compatible ASE2SPRKKR versions on the client and remote host."
        )
