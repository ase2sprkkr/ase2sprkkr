import numpy as np

from ..unique_values import UniqueValuesMapping


def test_merge_empty_mapping_with_sized_iterator():
    mapping = UniqueValuesMapping.from_values([])

    merged = mapping.merge([])

    assert isinstance(merged.mapping, np.ndarray)
    assert merged.mapping.dtype == np.int32
    assert merged.mapping.size == 0
    assert merged.value_to_class_id == {}
