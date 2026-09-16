# pytest sorting stuff to run the fast tests first

def pytest_collection_modifyitems(items):
    items.sort(
        key=lambda item: item.get_closest_marker("slow") is not None
    )

def pytest_configure(config):
    config.addinivalue_line(
        "markers",
        "slow: marks tests as slow",
    )
